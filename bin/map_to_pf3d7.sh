#!/usr/bin/env bash
set -euo pipefail

DEMUX_DIR="$1"
REF="$2"
THREADS="$3"

if [ -z "${DEMUX_DIR}" ]; then
    echo "[error] demux dir not provided" >&2
    exit 1
fi

mkdir -p mapped_pf

# Find barcode directories (one level deep)
find "${DEMUX_DIR}" -maxdepth 2 -type d -name 'barcode*' | sort | while IFS= read -r bc; do
    bcname=$(basename "${bc}")

    # Collect fastq files (gzipped or not)
    fq_list=$(find "${bc}" -type f | grep -Ei '\.fq$|\.fq\.gz$|\.fastq$|\.fastq\.gz$' | sort)
    if [ -z "${fq_list}" ]; then
        echo "[warn] No fastq under ${bc}" >&2
        continue
    fi

    # If any file is gzipped, use zcat; otherwise cat
    if echo "${fq_list}" | grep -qE '\.gz$'; then
        cat_cmd="zcat"
    else
        cat_cmd="cat"
    fi

    # Stream reads into minimap2 and pipe to samtools to produce a sorted BAM
    ${cat_cmd} ${fq_list} | minimap2 -t ${THREADS} -ax map-ont "${REF}" - \
        | samtools sort -@ ${THREADS} -o mapped_pf/${bcname}.bam -

    samtools index -@ ${THREADS} mapped_pf/${bcname}.bam
done
