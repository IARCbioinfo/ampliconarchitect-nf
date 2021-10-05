#!/bin/bash
#SBATCH --job-name="CNVkit"
#SBATCH --partition=low_p
#SBATCH -c 1
#SBATCH --mem-per-cpu=5000

g=`awk NR==${SLURM_ARRAY_TASK_ID} tumors.txt | awk '{print "SAMPLE="$1" CRAM="$2}'`
CONTAINER=/data/scratch/digenovaa/mesomics/AMPARCHITECT/aarchitect_v3.0.sif
#echo $g
singularity exec ${CONTAINER} make -f CNVKIT.mk  $g all
