#!/usr/bin/env bash
#
# Step one of this project. Reports which MHC class I prediction algorithms the
# pvactools container can actually answer for, without any licensed external
# tool. Everything else in the repo depends on the answer.
#
# Usage:  bin/smoke_test.sh [image]
#
# Writes a summary to stdout and a machine-readable list to
# results/smoke_test/available_predictors.txt
#
set -uo pipefail

IMAGE="${1:-griffithlab/pvactools:7.1.6}"
OUTDIR="results/smoke_test"
TEST_ALLELE="HLA-A*02:01"

# Class I predictors that the docs say are bundled in the image. MixMHCpred and
# PRIME are deliberately absent: they need a separate academic license.
PREDICTORS=(
  MHCflurry
  MHCflurryEL
  MHCnuggetsI
  BigMHC_EL
  BigMHC_IM
  DeepImmuno
  NetMHC
  NetMHCpan
  NetMHCpanEL
  NetMHCcons
  PickPocket
  SMM
  SMMPMBEC
)

PLATFORM=""
if [[ "$(uname -m)" == "arm64" || "$(uname -m)" == "aarch64" ]]; then
  PLATFORM="--platform linux/amd64"
  echo "note: arm64 host detected, running image under linux/amd64 emulation"
fi

mkdir -p "$OUTDIR"
AVAILABLE="$OUTDIR/available_predictors.txt"
: > "$AVAILABLE"

run_in_image() {
  # shellcheck disable=SC2086
  docker run --rm $PLATFORM "$IMAGE" "$@" 2>&1
}

echo "=== image: $IMAGE ==="
if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
  echo "pulling image (~9.4 GB, this takes a while)..."
  # shellcheck disable=SC2086
  docker pull $PLATFORM "$IMAGE" || { echo "FAIL: could not pull $IMAGE"; exit 1; }
fi

echo
echo "=== pvactools version ==="
run_in_image pvacseq --version || echo "WARN: version call failed"

echo
echo "=== bundled IEDB install at /opt/iedb ==="
# The NetMHC-family, PickPocket and SMM predictors come from here. If this
# directory is missing or empty, those predictors will not be usable.
run_in_image ls /opt/iedb || echo "WARN: /opt/iedb not listable"

echo
echo "=== predictor availability (probe: valid_alleles for $TEST_ALLELE) ==="
printf '%-16s %s\n' "PREDICTOR" "RESULT"
for p in "${PREDICTORS[@]}"; do
  out="$(run_in_image pvacseq valid_alleles -p "$p")"
  rc=$?
  if [[ $rc -eq 0 ]] && grep -qF "$TEST_ALLELE" <<<"$out"; then
    n="$(grep -c . <<<"$out")"
    printf '%-16s %s\n' "$p" "OK (${n} alleles)"
    echo "$p" >> "$AVAILABLE"
  elif [[ $rc -eq 0 ]]; then
    printf '%-16s %s\n' "$p" "responds, but test allele absent"
    echo "$p" >> "$AVAILABLE"
  else
    printf '%-16s %s\n' "$p" "UNAVAILABLE ($(head -1 <<<"$out" | cut -c1-60))"
  fi
done

echo
n_ok="$(grep -c . "$AVAILABLE" || true)"
echo "=== $n_ok predictor(s) available ==="
echo "written to: $AVAILABLE"
echo
if [[ "$n_ok" -lt 2 ]]; then
  echo "STOP: fewer than two predictors available. There is no comparison to make."
  echo "Re-read docs/smoke-test.md before continuing."
  exit 1
fi
if [[ "$n_ok" -lt 5 ]]; then
  echo "NOTE: a reduced predictor set. The study still works, but say so in the"
  echo "      write-up and drop the per-allele question if coverage is thin."
fi
echo "Next: set params.predictors in nextflow.config from the list above."
