#!/usr/bin/env bash
set -euo pipefail

# Create directory for GIAB HG002 reads
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$SCRIPT_DIR/giab_hg002_reads"
cd "$SCRIPT_DIR/giab_hg002_reads"

# Base URL for GIAB HG002 reads
BASE_URL="https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads/"

download_if_missing() {
    local url="$1"
    local destination="$2"
    local partial="${destination}.partial"

    if [[ -s "$destination" ]]; then
        echo "Already exists: $destination"
        return
    fi

    if command -v wget >/dev/null 2>&1; then
        wget -c -O "$partial" "$url"
    elif command -v curl >/dev/null 2>&1; then
        curl --fail --location --retry 3 --continue-at - \
            --output "$partial" "$url"
    else
        echo "ERROR: neither wget nor curl is available for downloading GIAB reads." >&2
        exit 1
    fi

    mv "$partial" "$destination"
}

# Download L001 lane files (001-017)
for i in {1..17}; do
    chunk=$(printf "%03d" $i)
    echo "Downloading L001 chunk ${chunk}..."
    download_if_missing "${BASE_URL}D1_S1_L001_R1_${chunk}.fastq.gz" "D1_S1_L001_R1_${chunk}.fastq.gz"
    download_if_missing "${BASE_URL}D1_S1_L001_R2_${chunk}.fastq.gz" "D1_S1_L001_R2_${chunk}.fastq.gz"
done

# Download L002 lane files (001-017)
for i in {1..17}; do
    chunk=$(printf "%03d" $i)
    echo "Downloading L002 chunk ${chunk}..."
    download_if_missing "${BASE_URL}D1_S1_L002_R1_${chunk}.fastq.gz" "D1_S1_L002_R1_${chunk}.fastq.gz"
    download_if_missing "${BASE_URL}D1_S1_L002_R2_${chunk}.fastq.gz" "D1_S1_L002_R2_${chunk}.fastq.gz"
done

merge_if_missing() {
    local destination="$1"
    shift

    if [[ -s "$destination" ]]; then
        echo "Already exists: $destination"
        return
    fi

    local partial="${destination}.partial"
    cat "$@" > "$partial"
    mv "$partial" "$destination"
}

echo "Creating merged FASTQs expected by run_giab.sh..."
merge_if_missing HG002_R1.fq.gz D1_S1_L00*_R1_*.fastq.gz
merge_if_missing HG002_R2.fq.gz D1_S1_L00*_R2_*.fastq.gz

echo "Done."
echo "  $SCRIPT_DIR/giab_hg002_reads/HG002_R1.fq.gz"
echo "  $SCRIPT_DIR/giab_hg002_reads/HG002_R2.fq.gz"
