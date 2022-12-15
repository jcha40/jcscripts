#!/bin/bash
#SBATCH --job-name=bulk_pileup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=6-822:6
#SBATCH --mem=16GB
#SBATCH --time=2:00:00

set -exo pipefail

module load anaconda3
source activate julia

SCIDXLIST=/storage/home/jsc6015/work/deep_motif_discovery/scidx_filepaths.txt
OUTDIR=/scratch/jsc6015/yep_heatmap
SCRIPTPATH=/storage/home/jsc6015/work/deep_motif_discovery/pileup.jl
BEDLIST=/storage/home/jsc6015/work/deep_motif_discovery/bed_filepaths.txt
CHROMSIZES=/storage/home/jsc6015/work/sacCer3/sacCer3_corrected.chrom.sizes

SCIDXS=$(head -n $SLURM_ARRAY_TASK_ID $SCIDXLIST | tail -n 6)
for SCIDX in $SCIDXS; do
  OUT=$OUTDIR/$(basename $SCIDX .tab).pileup.txt
  julia $SCRIPTPATH -s $SCIDX -b $BEDLIST -c $CHROMSIZES -o $OUT
done
