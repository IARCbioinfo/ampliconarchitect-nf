.DELETE_ON_ERROR:

#programs
CNVKIT=/home/digenovaa/.local/bin/cnvkit.py
SSEED=/data/scratch/digenovaa/mesomics/AMPARCHITECT/CNVKIT_MESO/select_seeds_cnvkit.py
AA_CODE=/home/programs/AmpliconArchitect-master/src
MOSEKLM_LICENSE_FILE=/data/scratch/digenovaa/mesomics/AMPARCHITECT/mosek/mosek.lic
AA_DATA_REPO=/home/digenovaa/scratch/mesomics/AMPARCHITECT/CNVKIT_MESO/
export AA_DATA_REPO
export MOSEKLM_LICENSE_FILE



#compute the copy segments for tumor-only mode
bname=$(notdir $(basename ${CRAM}))
#${SAMPLE}/$(basename ${CRAM}).cnr:
${SAMPLE}/${bname}.cnr:
	${CNVKIT} batch -m wgs -y -r GRCh38/GRCh38_cnvkit_filtered_ref.cnn -p 1 -d ${SAMPLE} ${CRAM}

#compute copy number calls for tumor-only mode
#${SAMPLE}/$(basename ${CRAM}).cns: ${SAMPLE}/$(basename ${CRAM}).cnr
${SAMPLE}/${bname}.segment.cns: ${SAMPLE}/${bname}.cnr
	${CNVKIT} segment  -p 1 -o $@ $<
# we select the seeds for AA
${SAMPLE}/${bname}_CNV_GAIN.bed:${SAMPLE}/${bname}.segment.cns
	python2 ${SSEED} -o ${SAMPLE} --sorted_bam ${CRAM}

#we select the seeds with amplidied intervals
${SAMPLE}.cnvkit.aa.bed:${SAMPLE}/${bname}_CNV_GAIN.bed
	python2 ${AA_CODE}/amplified_intervals.py --gain 4.5 --cnsize_min 50000 --ref GRCh38 --bed $<  --out $(subst .bed,,$@)

#we run AA
${SAMPLE}.cnvkit.aa.done:${SAMPLE}.cnvkit.aa.bed
	perl run-aa.pl $<  $(subst .cnvkit.aa.bed,,$<) > $(subst .cnvkit.aa.bed,,$<).aa.log 2> $(subst .cnvkit.aa.bed,,$<).aa.err

all : ${SAMPLE}/${bname}.segment.cns ${SAMPLE}/${bname}_CNV_GAIN.bed ${SAMPLE}.cnvkit.aa.bed ${SAMPLE}.cnvkit.aa.done

