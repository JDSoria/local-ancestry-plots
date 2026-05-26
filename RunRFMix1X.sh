#!/bin/bash
#SBATCH --job-name=RFMix_1X_submit
#SBATCH --partition=mono
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=logRFMix1X/submit_%j.out
#SBATCH --error=logRFMix1X/submit_%j.err

set -euo pipefail

# This launcher submits one RFMix job per chromosome for lcWGS/1X data.
# It is intended to be run on the cluster with: sbatch RunRFMix1X.sh

BASE=${BASE:-/home/dsoria/PoblAR/LocalAncestryMuestrasPiloto}
INPUT=${INPUT:-${BASE}/Script/Input}
OUTDIR=${OUTDIR:-${BASE}/1X}
SINGULARITY_IMAGE=${SINGULARITY_IMAGE:-/home/dsoria/miniconda3/rfmix.sif}
BIND_PATH=${BIND_PATH:-/home/pluisi/}

QUERY_VCF=${QUERY_VCF:-${INPUT}/LocalAncestry.33IndYControl.vcf.gz}
REFERENCE_VCF=${REFERENCE_VCF:-/home/pluisi/PoblAR/Analyse100Genomas_1X/WithRef_Downsampled/LocalAncestry.GENO0.001.NoMonomorphic.1X.vcf.gz}
SAMPLE_MAP=${SAMPLE_MAP:-${INPUT}/SampleMapF.txt}
GENETIC_MAP_DIR=${GENETIC_MAP_DIR:-${INPUT}/MapasGeneticos}

CRF_SPACING=${CRF_SPACING:-0.5}
RF_WINDOW=${RF_WINDOW:-0.5}
GENERATIONS=${GENERATIONS:-11}
EM_ITERATIONS=${EM_ITERATIONS:-1}
NUM_TREES=${NUM_TREES:-5}

mkdir -p "${OUTDIR}" logRFMix1X

[[ -f "${QUERY_VCF}" ]] || { echo "Missing query VCF: ${QUERY_VCF}" >&2; exit 1; }
[[ -f "${REFERENCE_VCF}" ]] || { echo "Missing reference VCF: ${REFERENCE_VCF}" >&2; exit 1; }
[[ -f "${SAMPLE_MAP}" ]] || { echo "Missing sample map: ${SAMPLE_MAP}" >&2; exit 1; }
[[ -f "${SINGULARITY_IMAGE}" ]] || { echo "Missing Singularity image: ${SINGULARITY_IMAGE}" >&2; exit 1; }

for chr in $(seq 1 22); do
  chr_name="chr${chr}"
  map_file="${GENETIC_MAP_DIR}/${chr_name}.b38.gmap.final.corrected.modified"
  out_prefix="${OUTDIR}/${chr_name}/cromosomaPPieFinal${chr}"

  [[ -f "${map_file}" ]] || { echo "Missing genetic map: ${map_file}" >&2; exit 1; }
  mkdir -p "$(dirname "${out_prefix}")"

  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=RFMix1X_${chr_name}
#SBATCH --partition=mono
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=${out_prefix}.out
#SBATCH --error=${out_prefix}.err

set -euo pipefail

/usr/bin/singularity exec --bind ${BIND_PATH} ${SINGULARITY_IMAGE} \
  rfmix \
  -f ${QUERY_VCF} \
  -r ${REFERENCE_VCF} \
  -m ${SAMPLE_MAP} \
  -g ${map_file} \
  -o ${out_prefix} \
  --chromosome=${chr_name} \
  -n ${NUM_TREES} \
  -s ${RF_WINDOW} \
  -c ${CRF_SPACING} \
  -G ${GENERATIONS} \
  -e ${EM_ITERATIONS} \
  --reanalyze-reference \
  > ${out_prefix}.log 2>&1
EOF
done
