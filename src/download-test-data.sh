#!/usr/bin/env sh

# OUTDIR="put here your preferred EXISTING download directory"
OUTDIR="/gpfs/projects/bsc40/project/pipelines/genepred/test_data"

module load sratoolkit-precomp/2.9.6

cd ${OUTDIR}

unset X
declare -a X=(SRR18362107 SRR13978640 SRR19090837 SRR15178417 SRR14687189 SRR14687233 SRR5989374 SRR6352887 SRR6352891)

for i in ${X[@]}
do
  fasterq-dump \
  --outdir ${OUTDIR} \
  --threads 4 \
  --progress \
  --split-files \
  ${i}
done
