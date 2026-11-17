#!/usr/bin/env nextflow
//
// Neoantigen candidate generation across MHC class I predictors.
//
//   MC3 MAF ──subset──> per-sample MAF ──maf2vcf──> VCF ──VEP──> annotated VCF
//        └──(optional expression annotation)──> pVACseq × {predictors}
//
// The fan-out over predictors is the point of the study: every predictor sees
// an identical annotated VCF, so differences downstream are attributable to
// the predictor alone.

nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
    Usage:
      nextflow run . -profile docker \\
        --input assets/samplesheet.csv \\
        --maf data/mc3.v0.2.8.PUBLIC.maf.gz \\
        --ref_fasta ref/GRCh37.fa \\
        --vep_cache ref/vep_cache \\
        --vep_plugins ref/VEP_plugins

    Run bin/smoke_test.sh first, and set --predictors from its output.
    """.stripIndent()
}

// ---------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------

process SUBSET_MAF {
    // MC3 is one ~2 GB MAF covering all of TCGA. Pull out one sample's rows and
    // keep only PASS calls, which is the documented convention for this file.
    tag "$sample"
    label 'process_low'
    publishDir "${params.outdir}/maf", mode: 'copy'

    input:
    tuple val(sample), path(maf)

    output:
    tuple val(sample), path("${sample}.maf"), emit: maf

    script:
    """
    set -euo pipefail

    reader() { if [[ "$maf" == *.gz ]]; then zcat "$maf"; else cat "$maf"; fi; }

    # Header is the first line that is not a '#' comment.
    reader | awk 'NR<=10 && /^#/ {next} {print; exit}' > header.tsv

    # Resolve the column indices by name rather than position: MAF column order
    # is not guaranteed stable between releases.
    bc=\$(awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="Tumor_Sample_Barcode") print i; exit}' header.tsv)
    fc=\$(awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="FILTER") print i; exit}' header.tsv)
    if [[ -z "\$bc" ]]; then echo "no Tumor_Sample_Barcode column" >&2; exit 1; fi

    cp header.tsv "${sample}.maf"
    reader | awk -F'\\t' -v OFS='\\t' -v bc="\$bc" -v fc="\$fc" -v s="${sample}" '
        \$0 ~ /^#/ {next}
        \$bc == s {
            if (fc != "" && \$fc != "PASS" && \$fc != "") next
            print
        }' >> "${sample}.maf"

    n=\$(( \$(wc -l < "${sample}.maf") - 1 ))
    echo "${sample}: \$n PASS variants" >&2
    if [[ "\$n" -eq 0 ]]; then
        echo "ERROR: no variants for ${sample}. Check the barcode matches the MAF." >&2
        exit 1
    fi
    """

    stub:
    """
    printf 'Hugo_Symbol\\tTumor_Sample_Barcode\\tFILTER\\n' > "${sample}.maf"
    printf 'TP53\\t%s\\tPASS\\n' "${sample}" >> "${sample}.maf"
    """
}

process MAF2VCF {
    tag "$sample"
    label 'process_low'
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    tuple val(sample), path(maf)
    path ref_fasta

    output:
    tuple val(sample), path("${sample}.vcf"), emit: vcf

    script:
    """
    set -euo pipefail
    maf2vcf.pl \\
        --input-maf "$maf" \\
        --output-dir . \\
        --output-vcf "${sample}.vcf" \\
        --ref-fasta "$ref_fasta"
    """

    stub:
    """
    touch "${sample}.vcf"
    """
}

process VEP_ANNOTATE {
    // pVACseq requires the Wildtype and Frameshift plugins from
    // griffithlab/VEP_plugins. Without them the annotated VCF lacks the
    // protein sequence fields pVACseq reads, and it fails later with an
    // unhelpful message. This is the step most likely to cost an evening.
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/vep", mode: 'copy'

    input:
    tuple val(sample), path(vcf)
    path ref_fasta
    path vep_cache
    path vep_plugins

    output:
    tuple val(sample), path("${sample}.vep.vcf"), emit: vcf
    path "${sample}.vep.summary.html", optional: true

    script:
    """
    set -euo pipefail
    vep \\
        --input_file "$vcf" \\
        --output_file "${sample}.vep.vcf" \\
        --stats_file "${sample}.vep.summary.html" \\
        --format vcf --vcf --offline --cache \\
        --dir_cache "$vep_cache" \\
        --cache_version ${params.vep_cache_version} \\
        --assembly ${params.vep_assembly} \\
        --species homo_sapiens \\
        --fasta "$ref_fasta" \\
        --dir_plugins "$vep_plugins" \\
        --plugin Wildtype \\
        --plugin Frameshift \\
        --symbol --terms SO --tsl --biotype --hgvs \\
        --transcript_version \\
        --pick --flag_pick \\
        --fork ${task.cpus} \\
        --force_overwrite
    """

    stub:
    """
    touch "${sample}.vep.vcf" "${sample}.vep.summary.html"
    """
}

process ANNOTATE_EXPRESSION {
    // Adds gene-level expression to the VCF so pVACseq's --expn-val filter has
    // something to act on. Only runs for samples with an expression file;
    // samples without one go straight through, which is also the control arm
    // for the filter-sensitivity question.
    tag "$sample"
    label 'process_low'
    publishDir "${params.outdir}/vep", mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(expression)

    output:
    tuple val(sample), path("${sample}.expr.vcf"), emit: vcf

    script:
    // NOTE: verify the id/expression column names against the GDC STAR-Counts
    // TSV you actually download before trusting this.
    """
    set -euo pipefail
    vcf-expression-annotator \\
        -i gene_id \\
        -e tpm_unstranded \\
        -s "${sample}" \\
        -o "${sample}.expr.vcf" \\
        "$vcf" "$expression" custom gene
    """

    stub:
    """
    touch "${sample}.expr.vcf"
    """
}

process PVACSEQ {
    // One task per (sample, predictor). Identical input VCF across predictors.
    tag "${sample}:${predictor}"
    label 'process_high'
    // A dynamic publish path must be a closure so it is resolved per task.
    publishDir { "${params.outdir}/pvacseq/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(vcf), val(hla), val(predictor)

    output:
    tuple val(sample), val(predictor), path("${predictor}/**"), emit: results

    script:
    def expn = params.expn_val ? "--expn-val ${params.expn_val}" : ''
    """
    set -euo pipefail
    pvacseq run \\
        "$vcf" \\
        "${sample}" \\
        "${hla}" \\
        ${predictor} \\
        "${predictor}" \\
        --iedb-install-directory /opt/iedb \\
        -e1 ${params.epitope_lengths} \\
        -b ${params.binding_threshold} \\
        ${expn} \\
        --pass-only \\
        -t ${task.cpus}
    """

    stub:
    """
    mkdir -p "${predictor}/MHC_Class_I"
    touch "${predictor}/MHC_Class_I/${sample}.filtered.tsv"
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    if (params.help) { helpMessage(); exit 0 }

    // Required inputs, checked up front rather than three processes deep.
    def required = ['input', 'maf', 'ref_fasta', 'vep_cache', 'vep_plugins']
    def missing  = required.findAll { !params[it] }
    if (missing) {
        exit 1, "Missing required parameter(s): ${missing.collect{ "--$it" }.join(', ')}\n" +
                "Run with --help for usage."
    }

    ch_ref     = Channel.fromPath(params.ref_fasta,  checkIfExists: true).first()
    ch_cache   = Channel.fromPath(params.vep_cache,  checkIfExists: true).first()
    ch_plugins = Channel.fromPath(params.vep_plugins, checkIfExists: true).first()
    ch_maf     = Channel.fromPath(params.maf,        checkIfExists: true).first()

    // Samplesheet: sample,hla,expression  (expression may be blank)
    ch_samples = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample?.trim()) exit 1, "Samplesheet row missing 'sample': $row"
            if (!row.hla?.trim())    exit 1, "Samplesheet row missing 'hla': $row"
            // pVACseq wants a comma-separated allele list; accept ';' too.
            def hla = row.hla.trim().replaceAll(';', ',').replaceAll(/\s+/, '')
            def expr = row.expression?.trim() ? file(row.expression.trim(), checkIfExists: true) : null
            tuple(row.sample.trim(), hla, expr)
        }

    // Pair every sample name with the one shared MC3 MAF.
    SUBSET_MAF( ch_samples.map { s, hla, e -> s }.combine(ch_maf) )

    MAF2VCF( SUBSET_MAF.out.maf, ch_ref )

    VEP_ANNOTATE( MAF2VCF.out.vcf, ch_ref, ch_cache, ch_plugins )

    // Split by whether the sample has an expression file.
    ch_branched = VEP_ANNOTATE.out.vcf
        .join( ch_samples.map { s, hla, e -> tuple(s, e) } )
        .branch { s, vcf, expr ->
            with_expr:    expr != null
            without_expr: true
        }

    ANNOTATE_EXPRESSION( ch_branched.with_expr )

    ch_ready = ANNOTATE_EXPRESSION.out.vcf
        .mix( ch_branched.without_expr.map { s, vcf, expr -> tuple(s, vcf) } )

    // Fan out: every predictor against every sample's identical VCF.
    ch_predictors = Channel.fromList(
        params.predictors.toString().split(',').collect { it.trim() }.findAll { it }
    )

    PVACSEQ(
        ch_ready
            .join( ch_samples.map { s, hla, e -> tuple(s, hla) } )
            .combine( ch_predictors )
    )
}
// The completion handler lives in nextflow.config: strict syntax (Nextflow 25+)
// forbids top-level statements here, and inside a workflow body the name
// `workflow` is shadowed by the workflow scope, so `workflow.success` is null.
