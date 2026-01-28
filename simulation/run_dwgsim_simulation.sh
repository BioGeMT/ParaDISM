#!/usr/bin/env bash
set -euo pipefail

# Small helper script to:
#  1) Ensure DWGSIM is available (use system dwgsim if present, otherwise
#     clone and build a local copy from GitHub).
#  2) Run a PKD1/PKD1P* read simulation against ParaDISM/ref.fa.
#
# You can override key parameters via environment variables:
#   REF        : reference FASTA (default: ParaDISM/ref.fa)
#   OUT_DIR    : output directory for DWGSIM files (default: ParaDISM/pkd_reads_dwgsim)
#   OUT_PREFIX : DWGSIM output prefix (default: ${OUT_DIR}/pkd_reads_dwgsim)
#   NUM_READS  : number of read pairs (default: 100000)
#   READ_LEN   : read length for R1/R2 (default: 150)
#   FRAG_MEAN  : mean insert size (default: 350)
#   FRAG_SD    : insert size std dev (default: 35)
#   ERROR_RATE : per-base error rate for R1/R2 (default: 0.01)
#   SNP_RATE   : per-base SNP rate (default: 0.005)
#   INDEL_RATE : per-base indel rate (default: 0.0005)
#   INDEL_EXT  : indel extension probability (default: 0.5)
#   DWGSIM_DIR : location to clone/build DWGSIM (default: ../dwgsim relative to this script)
#   ZLIB_PREFIX: prefix for a locally built zlib (default: $HOME/local/zlib)
#
# Example:
#   cd ParaDISM
#   bash simulation/run_dwgsim_simulation.sh
#
# After running, FASTQs will be in ${OUT_PREFIX}.bwa.read1.fastq.gz and
# ${OUT_PREFIX}.bwa.read2.fastq.gz

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

REF="${REF:-${PROJECT_ROOT}/ref.fa}"
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/pkd_reads_dwgsim}"
OUT_PREFIX="${OUT_PREFIX:-${OUT_DIR}/pkd_reads_dwgsim}"
DWGSIM_DIR="${DWGSIM_DIR:-${SCRIPT_DIR}/dwgsim}"
ZLIB_PREFIX="${ZLIB_PREFIX:-${HOME}/local/zlib}"

NUM_READS="${NUM_READS:-100000}"
READ_LEN="${READ_LEN:-150}"
FRAG_MEAN="${FRAG_MEAN:-350}"
FRAG_SD="${FRAG_SD:-35}"
ERROR_RATE="${ERROR_RATE:-0.01}"
SNP_RATE="${SNP_RATE:-0.005}"
INDEL_RATE="${INDEL_RATE:-0.0005}"
INDEL_EXT="${INDEL_EXT:-0.5}"
# Space-separated list of seeds to simulate. Example: SEEDS="1 2 3".
SEEDS="${SEEDS:-1}"

ZLIB_CFLAGS=""
ZLIB_LDFLAGS="-lz"

ensure_tool() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Error: required tool '$1' not found in PATH." >&2
    exit 1
  fi
}

ensure_zlib() {
  # If system zlib headers are available, use them.
  if echo '#include <zlib.h>' | gcc -E - >/dev/null 2>&1; then
    ZLIB_CFLAGS=""
    ZLIB_LDFLAGS="-lz"
    echo "[run_dwgsim] Using system zlib."
    return
  fi

  # If a local zlib install already exists, use it.
  if [ -f "${ZLIB_PREFIX}/include/zlib.h" ]; then
    ZLIB_CFLAGS="-I${ZLIB_PREFIX}/include"
    ZLIB_LDFLAGS="-L${ZLIB_PREFIX}/lib -lz"
    export LD_LIBRARY_PATH="${ZLIB_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
    echo "[run_dwgsim] Using existing local zlib at ${ZLIB_PREFIX}."
    return
  fi

  echo "[run_dwgsim] No zlib headers found; building local zlib under ${ZLIB_PREFIX}."

  mkdir -p "${ZLIB_PREFIX}"
  local deps_root="${PROJECT_ROOT}/../deps"
  mkdir -p "${deps_root}"
  cd "${deps_root}"

  local zver="1.3.1"
  local ztar="zlib-${zver}.tar.gz"

  if [ ! -d "zlib-${zver}" ]; then
    if [ ! -f "${ztar}" ]; then
      if command -v curl >/dev/null 2>&1; then
        curl -L "https://zlib.net/${ztar}" -o "${ztar}"
      elif command -v wget >/dev/null 2>&1; then
        wget "https://zlib.net/${ztar}"
      else
        echo "Error: neither curl nor wget available to download zlib." >&2
        exit 1
      fi
    fi
    tar xzf "${ztar}"
  fi

  cd "zlib-${zver}"
  ./configure --prefix="${ZLIB_PREFIX}"
  make -j"$(command -v nproc >/dev/null 2>&1 && nproc || echo 2)"
  make install

  ZLIB_CFLAGS="-I${ZLIB_PREFIX}/include"
  ZLIB_LDFLAGS="-L${ZLIB_PREFIX}/lib -lz"
  export LD_LIBRARY_PATH="${ZLIB_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

  echo "[run_dwgsim] Built local zlib at ${ZLIB_PREFIX}."
}

DWGSIM_BIN=""

ensure_dwgsim() {
  # Prefer a system-installed dwgsim if present.
  if command -v dwgsim >/dev/null 2>&1; then
    DWGSIM_BIN="$(command -v dwgsim)"
    echo "[run_dwgsim] Using DWGSIM from PATH: ${DWGSIM_BIN}"
    return
  fi

  mkdir -p "${DWGSIM_DIR}"

  if [ ! -f "${DWGSIM_DIR}/src/dwgsim.c" ]; then
    echo "[run_dwgsim] Cloning DWGSIM into ${DWGSIM_DIR}."
    git clone https://github.com/nh13/DWGSIM.git "${DWGSIM_DIR}"
  else
    echo "[run_dwgsim] DWGSIM source already present at ${DWGSIM_DIR}."
  fi

  cd "${DWGSIM_DIR}"

  if [ ! -x "./dwgsim" ]; then
    echo "[run_dwgsim] Compiling local DWGSIM binary."
    rm -f src/*.o dwgsim
    (
      cd src
      gcc -c -g -Wall -O3 -I.. ${ZLIB_CFLAGS} -DPACKAGE_VERSION=\"0.1.16-dev\" \
        dwgsim_opt.c mut.c contigs.c regions_bed.c \
        mut_txt.c mut_bed.c mut_vcf.c mut_input.c dwgsim.c
    )
    gcc -g -Wall -O3 -o dwgsim \
      src/dwgsim_opt.o src/mut.o src/contigs.o src/regions_bed.o \
      src/mut_txt.o src/mut_bed.o src/mut_vcf.o src/mut_input.o src/dwgsim.o \
      ${ZLIB_LDFLAGS} -lm -lpthread
  fi

  DWGSIM_BIN="${DWGSIM_DIR}/dwgsim"
  echo "[run_dwgsim] Using local DWGSIM at ${DWGSIM_BIN}."
}

main() {
  ensure_tool gcc
  ensure_tool git

  if [ ! -f "${REF}" ]; then
    echo "Error: reference FASTA not found at ${REF}" >&2
    exit 1
  fi

  mkdir -p "${OUT_DIR}"

  ensure_zlib
  ensure_dwgsim
  
  # Expand SEEDS into a list. Supports either:
  #   SEEDS="1 2 3"      (explicit list)
  #   SEEDS="1-1000"     (numeric range)
  local seeds_expanded=()
  if [[ "${SEEDS}" =~ ^[0-9]+-[0-9]+$ ]]; then
    local start="${SEEDS%-*}"
    local end="${SEEDS#*-}"
    for s in $(seq "${start}" "${end}"); do
      seeds_expanded+=("${s}")
    done
  else
    # Treat SEEDS as a plain space‑separated list
    for s in ${SEEDS}; do
      seeds_expanded+=("${s}")
    done
  fi

  for seed in "${seeds_expanded[@]}"; do
    local prefix_seed="${OUT_PREFIX}_seed${seed}"

    echo "[run_dwgsim] Running DWGSIM simulation (seed=${seed})."
    echo "  Reference : ${REF}"
    echo "  Output    : ${prefix_seed}"
    echo "  Pairs     : ${NUM_READS}"
    echo "  Read len  : ${READ_LEN}"
    echo "  Insert    : mean=${FRAG_MEAN}, sd=${FRAG_SD}"
    echo "  Error     : ${ERROR_RATE}"
    echo "  SNP       : ${SNP_RATE}"
    echo "  Indel     : rate=${INDEL_RATE}, ext=${INDEL_EXT}"

    "${DWGSIM_BIN}" \
      -z "${seed}" \
      -N "${NUM_READS}" \
      -1 "${READ_LEN}" -2 "${READ_LEN}" \
      -d "${FRAG_MEAN}" -s "${FRAG_SD}" \
      -y 0 \
      -e "${ERROR_RATE}" -E "${ERROR_RATE}" \
      -r "${SNP_RATE}" -R "${INDEL_RATE}" -X "${INDEL_EXT}" \
      "${REF}" "${prefix_seed}"

    # Rename DWGSIM's default BWA-style FASTQ names to simpler R1/R2 names.
    local orig_r1="${prefix_seed}.bwa.read1.fastq.gz"
    local orig_r2="${prefix_seed}.bwa.read2.fastq.gz"
    local new_r1="${prefix_seed}_R1.fastq.gz"
    local new_r2="${prefix_seed}_R2.fastq.gz"

    if [ -f "${orig_r1}" ]; then
      mv -f "${orig_r1}" "${new_r1}"
    fi
    if [ -f "${orig_r2}" ]; then
      mv -f "${orig_r2}" "${new_r2}"
    fi

    echo
    echo "[run_dwgsim] Simulation complete (seed=${seed}). Key outputs:"
    echo "  ${new_r1}"
    echo "  ${new_r2}"
    echo "  ${prefix_seed}.mutations.txt"
    echo "  ${prefix_seed}.mutations.vcf"
  done
}

main "$@"
