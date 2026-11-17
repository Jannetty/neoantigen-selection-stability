# Step one: which predictors do we actually have?

Everything in this repo rests on one unverified assumption: that the
`griffithlab/pvactools` image lets you run the NetMHC-family predictors,
PickPocket, SMM and SMMPMBEC without applying for a license. That claim comes
from the pVACtools documentation, which states that the image bundles a
standalone IEDB installation at `/opt/iedb` alongside MHCflurry, MHCnuggets,
BigMHC and DeepImmuno, and that using pVACtools implicitly accepts the
underlying tools' license terms.

It has not been checked against the image. Check it before writing any analysis
code, because the answer determines whether this study is worth doing.

## Run it

```bash
bin/smoke_test.sh
```

Roughly 9.4 GB of image download on the first run. On Apple Silicon the script
adds `--platform linux/amd64` automatically.

## Reading the result

The script writes the predictors that responded to
`results/smoke_test/available_predictors.txt`.

| Outcome | What it means | What to do |
|---|---|---|
| ~10 predictors OK | The assumption holds. | Copy the list into `params.predictors`. Proceed as planned. |
| 4 predictors OK (MHCflurry, MHCnuggets, BigMHC, DeepImmuno) | `/opt/iedb` is absent or gated. | Still a viable study, and still double what nf-core/epitopeprediction offers license-free. Drop the per-allele question if allele coverage is thin, and say so in the write-up. |
| 2 or fewer | No comparison to make. | Stop. Reconsider the project before spending compute. |

If the NetMHC family is missing, note that the remaining predictors are *more*
independent of one another than the NetMHC family is internally, which makes the
concordance numbers more informative per predictor even though there are fewer
of them. Worth a sentence in the findings either way.

## Then: a one-sample end-to-end run

Before the full matrix, prove the chain works on a single sample with a single
predictor:

```bash
nextflow run . -profile docker,test \
  --input assets/samplesheet.csv \
  --maf data/mc3.v0.2.8.PUBLIC.maf.gz \
  --ref_fasta ref/GRCh37.fa \
  --vep_cache ref/vep_cache \
  --vep_plugins ref/VEP_plugins
```

Expected failure points, in the order you will hit them:

1. **Barcode mismatch.** `SUBSET_MAF` exits if a sample has zero variants. MC3
   `Tumor_Sample_Barcode` values are full barcodes; the HLA table keys on a
   shorter patient barcode. Match them deliberately.
2. **Genome build.** MC3 v0.2.8 is GRCh37. A GRCh38 cache or FASTA produces
   silent nonsense rather than an error. Confirm with
   `zcat mc3.v0.2.8.PUBLIC.maf.gz | head -2 | cut -f4`.
3. **Missing VEP plugins.** Without `Wildtype` and `Frameshift` from
   [griffithlab/VEP_plugins](https://github.com/griffithlab/VEP_plugins),
   pVACseq fails later with an unhelpful message rather than at the VEP step.
4. **Expression column names.** `ANNOTATE_EXPRESSION` assumes `gene_id` and
   `tpm_unstranded`. Check the header of the GDC STAR-Counts TSV you actually
   downloaded.

## Reference setup

```bash
mkdir -p ref data

# VEP cache (GRCh37, matching MC3) — large, do this once
curl -O https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh37.tar.gz
tar -xzf homo_sapiens_vep_116_GRCh37.tar.gz -C ref/vep_cache

# pVACtools VEP plugins
git clone https://github.com/griffithlab/VEP_plugins.git ref/VEP_plugins
```

Verify the release-116 GRCh37 cache path above before relying on it; Ensembl
moves GRCh37 assets between the main and `grch37` FTP trees between releases.
