#!/usr/bin/env nextflow


//help function for the tool
def show_help (){
  log.info IARC_Header()
  log.info tool_header()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run iarcbioinfo/purple-nf -singularity [OPTIONS]

    Mandatory arguments:
      --tumor_file		    [file] File containing list of tumor CRAM files to be processed [sample_id CRAM_file]
      --cram_dir         [dir]  directory where the CRAM file are stored (including indexed)
      --moseklm_license [file] License for mosek program
      --aa_repo_dir               [directory] directory containing the AA_REPO [See Readme.md]
    Optional arguments:
      --output_folder       [string] name of output folder [def:amplicon_results]
      --cpu                 [Integer]  Number of CPUs[def:2]
      --mem 		            [Integer] Max memory [def:8Gb]
      --debug               [boolean] enable debuging of the pipeline
      """.stripIndent()
}





//we display help information
if (params.help){ show_help(); exit 0;}
//we display the header of the tool
log.info IARC_Header()
log.info tool_header()
//Check mandatory parameters
assert (params.aa_repo_dir != null) : "please specify --aa_repo_dir  directory"
assert (params.tumor_file != null ) : "please specify --tumor_file"
assert (params.cram_dir != null ) : "please specify --cram_dir"
assert (params.moseklm_license != null ) : "please specify --moseklm_license"

//function that read the tumors to process from a tn_file
if(params.tumor_file){
  def cram = params.bam ? false:true
 tumors_crams = parse_tumorfile(params.tumor_file,params.cram_dir,cram)
 //we duplicate the tn_pairs channel
 tumors_crams.into {cnvkit_input; aa_input}
}

// //chanel for reference genome
// ref_fasta = Channel.value(file(params.ref)).ifEmpty{exit 1, "reference file not found: ${params.ref}"}
// ref_fai = Channel.value(file(params.ref+'.fai')).ifEmpty{exit 1, "index file not found: ${params.ref}.fai"}
// //BWA indexes for re-mapping MHC reads
// ref_sa  = file(params.ref+'.sa')
// ref_bwt =  file(params.ref+'.bwt')
// ref_ann =  file(params.ref+'.ann')
// ref_amb =  file(params.ref+'.amb')
// ref_pac =  file(params.ref+'.pac')

aa_repo_dir_path = file(params.aa_repo_dir)
mosek_license=file(params.moseklm_license)

print_params()


process cnvkit {
 cpus params.cpu
 memory params.mem+'G'

  publishDir params.output_folder+'/cnvkit/', mode: 'copy'
  input:
  set val(tumor_id), file(cram), file(cram_index) from cnvkit_input
  file (aa_repo_dir_path)

  output:
  set val(tumor_id), file("${tumor_id}/") into cnvkit_out
  set val(tumor_id), file("${tumor_id}.cnvkit.bed") into cnvkit_seeds

  script:
    if(params.debug == false)
       """
      #we run the cnvkit program
      cnvkit.py batch -m wgs -y -r ${aa_repo_dir_path}/GRCh38_cnvkit_filtered_ref.cnn -p 1 -d ${tumor_id} ${cram}
      #step 2
      cnvkit.py segment -p 1 -o ${tumor_id}/${cram.baseName}.segment.cns ${tumor_id}/${cram.baseName}.cnr
      #step 3 creates ${tumor_id}/${cram.baseName}_CNV_GAIN.bed
      python2 ${baseDir}/aux_scripts/select_seeds_cnvkit.py -o ${tumor_id} --sorted_bam ${cram}
      #we copy the cnv selected seeds
      cp ${tumor_id}/${cram.baseName}_CNV_GAIN.bed ${tumor_id}.cnvkit.bed
      """
    //debug mode
    else
      """
      echo cnvkit.py batch -m wgs -y -r ${aa_repo_dir_path}/GRCh38_cnvkit_filtered_ref.cnn -p 1 -d ${tumor_id} ${cram}
      #step 2
      echo cnvkit.py segment -p 1 -o ${tumor_id}/${cram.baseName}.segment.cns ${tumor_id}/${cram.baseName}.cnr
      #step 3 creates ${tumor_id}/${cram.baseName}_CNV_GAIN.bed
      echo python2 ${baseDir}/aux_scripts/select_seeds_cnvkit.py -o ${tumor_id} --sorted_bam ${cram}
      #we copy the cnv selected seeds
      echo cp ${tumor_id}/${cram.baseName}_CNV_GAIN.bed ${tumor_id}.cnvkit.bed

      mkdir ${tumor_id}
      touch ${tumor_id}.cnvkit.bed ${tumor_id}/${cram.baseName}_CNV_GAIN.bed
      touch ${tumor_id}/${cram.baseName}.segment.cns ${tumor_id}/${cram.baseName}.cnr
      """
}

//we select the ammplicon architect seeds
process amplified_intervals {
 cpus params.cpu
 memory params.mem+'G'

  publishDir params.output_folder+'/amplified_intervals/', mode: 'copy'
  input:
  set val(tumor_id), file(cnvkit_bed) from cnvkit_seeds
  file (aa_repo_dir_path)
  file (mosek_license)

  output:
  set val(tumor_id), file("${tumor_id}.cnvkit.aa.bed") into ampinter_seeds

  script:
    if(params.debug == false)
       """
      #we select the seeds with amplidied intervals using the ref genome
        export MOSEKLM_LICENSE_FILE=${mosek_license}
        export AA_DATA_REPO=\$PWD
      	python2 /home/programs/AmpliconArchitect-master/src/amplified_intervals.py \\
                  --gain 4.5 --cnsize_min 50000 --ref GRCh38 \\
                  --bed ${cnvkit_bed}  --out ${cnvkit_bed.baseName}.aa
      """
    //debug mode
    else
      """
      echo python2 /home/programs/AmpliconArchitect-master/src/amplified_intervals.py \\
                --gain 4.5 --cnsize_min 50000 --ref GRCh38 \\
                --bed ${cnvkit_bed}  --out ${cnvkit_bed.baseName}.aa.bed
      export MOSEKLM_LICENSE_FILE=${mosek_license}
      export AA_DATA_REPO=${aa_repo_dir_path}
      echo \$MOSEKLM_LICENSE_FILE > ${cnvkit_bed.baseName}.aa.bed
      echo \$AA_DATA_REPO >> ${cnvkit_bed.baseName}.aa.bed
      """
}

//we sync both crams+selected seeds
seed_crams=ampinter_seeds.join(aa_input)
//we run the amplicon architect tool
process amplicon_architect{
  cpus params.cpu
  memory params.mem+'G'

   publishDir params.output_folder+'/amplicon_predictions/', mode: 'copy'
   input:
   set val(tumor_id), file(cnvseeds),file(cram),file(cram_index) from seed_crams
   file (aa_repo_dir_path)
   file (mosek_license)

   output:
   file("${tumor_id}.aa") into aa_results
   set val(tumor_id), file("${tumor_id}_aa_summary.txt") into aa_results_log

   script:
     if(params.debug == false)
        """
         export MOSEKLM_LICENSE_FILE=\$PWD/mosek.lic
         export AA_DATA_REPO=\$PWD
         touch coverage.stats
       #we select the seeds with amplidied intervals using the ref genome
         mkdir ${tumor_id}.aa
         cd ${tumor_id}.aa
         #ln  -s ${aa_repo_dir_path} .
         n_s=`wc -l ../${cnvseeds}  | awk '{print \$1}'`
         if [[ \$n_s -gt 0 && \$n_s -lt 71 ]]
          then
            python2 /home/programs/AmpliconArchitect-master/src/AmpliconArchitect.py \\
                            --bed ../${cnvseeds} --bam ../${cram} --out ${tumor_id} --ref GRCh38
	    cp ${tumor_id}_summary.txt ../${tumor_id}_aa_summary.txt
          else
            echo "canceling amplicon architect run, ${tumor_id} with no (0) or many (>70) seeds (\$n_s)" > ../${tumor_id}_aa_summary.txt
          fi
         """
     //debug mode
     else
       """
       #we select the seeds with amplidied intervals using the ref genome
         export MOSEKLM_LICENSE_FILE=\$PWD/mosek.lic
         export AA_DATA_REPO=\$PWD
         touch coverage.stats
         mkdir ${tumor_id}.aa
         cd ${tumor_id}.aa
         n_s=`wc -l ../${cnvseeds}  | awk '{print \$1}'`
         if [[ \$n_s -gt 0 && \$n_s -lt 71 ]]
          then
            echo python2 /home/programs/AmpliconArchitect-master/src/AmpliconArchitect.py \\
                            --bed ../${cnvseeds} --bam ../${cram} --out ${tumor_id} --ref GRCh38
            echo "running amplicon architect, ${tumor_id} with the right number of seeds (\$n_s)" > ../${tumor_id}_aa_summary.txt
          else
            echo "canceling amplicon architect run, ${tumor_id} with no (0) or many (>70) seeds (\$n_s)" > ../${tumor_id}_aa_summary.txt
          fi
        """
}


//we gatther all the amplicon architect results to compute amplicon classes
process amplicon_classifier{
  cpus params.cpu
  memory params.mem+'G'

   publishDir params.output_folder+'/amplicon_classes/', mode: 'copy'
   input:
   file("*.aa") from aa_results.collect()
   file (aa_repo_dir_path)
   file (mosek_license)

   output:
   file("all_amplicons*") into aa_clases

   script:
     if(params.debug == false)
        """
         export MOSEKLM_LICENSE_FILE=\$PWD/mosek.lic
         export AA_DATA_REPO=\$PWD
         #we find all the amplicon results
         find \$PWD -name "*_cycles.txt" | sed 's/_cycles.txt//' | \\
         awk '{split(\$0,a,"/"); print a[length(a)]" "\$0"_cycles.txt "\$0"_graph.txt"}' > all_amplicons.txt
         #we run amplicon classifier
         python2 /home/programs/AmpliconClassifier/amplicon_classifier.py --ref GRCh38 --input all_amplicons.txt
       """
     //debug mode
     else
       """
       export MOSEKLM_LICENSE_FILE=\$PWD/mosek.lic
       export AA_DATA_REPO=\$PWD
       #we find all the amplicon results
       find \$PWD -name "*_cycles.txt" | sed 's/_cycles.txt//' | \\
       awk '{split(\$0,a,"/"); print a[length(a)]" "\$0"_cycles.txt "\$0"_graph.txt"}' > all_amplicons.txt
       #we run amplicon classifier
       echo python2 /home/programs/AmpliconClassifier/amplicon_classifier.py --ref GRCh38 --input all_amplicons.txt
       touch all_amplicons.txt all_amplicons_gene_list.tsv
       mkdir all_amplicons_classification_bed_files
       touch all_amplicons_amplicon_classification_profiles.tsv
       """
}

/*
*
* Functions to create channels from TSV or directories containing BAM/CRAM
*
*/

//we read the pairs from tn_file
def parse_tumorfile (tumor_file,path_cram,cram){
	    // FOR INPUT AS A TAB DELIMITED FILE
			def file_ext = cram ? '.crai':'.bai'
			//[sample cram]
      //id  cram
    def tumors_crams=Channel.fromPath(tumor_file)
      .splitCsv(header: true, sep: '\t', strip: true)
      .map{row -> [ row.id,
               file(path_cram + "/" + row.cram),
               file(path_cram + "/" + row.cram+file_ext)]}
      .ifEmpty{exit 1, "${tumor_file} was empty - no tumors supplied" }
	//we return the channel
  return tumors_crams
}

// print the calling parameter to the log and a log file
def print_params () {
  //software versions for v2.0
  def software_versions = ['cnvkit' : '0.9.8',
                           'ampliconarchitect' : '1.2']
  //we print the parameters
  log.info "\n"
  log.info "-\033[2m------------------Calling PARAMETERS--------------------\033[0m-"
  log.info params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
  log.info "-\033[2m--------------------------------------------------------\033[0m-"
  log.info "\n"
  log.info "-\033[2m------------------Software versions--------------------\033[0m-"
  log.info software_versions.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
  log.info "-\033[2m--------------------------------------------------------\033[0m-"
  log.info "\n"


  //we print the parameters to a log file
   def output_d = new File("${params.output_folder}/nf-pipeline_info/")
   if (!output_d.exists()) {
       output_d.mkdirs()
   }
   def output_tf = new File(output_d, "run_parameters_report.txt")
   def  report_params="------------------Calling PARAMETERS--------------------\n"
        report_params+= params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
        report_params+="\n--------------------------------------------------------\n"
        report_params+="\n------------------NEXTFLOW Metadata--------------------\n"
        report_params+="nextflow version : "+nextflow.version+"\n"
        report_params+="nextflow build   : "+nextflow.build+"\n"
        report_params+="Command line     : \n"+workflow.commandLine.split(" ").join(" \\\n")
        report_params+="\n--------------------------------------------------------\n"
        report_params+="-----------------Software versions--------------------\n"
        report_params+=software_versions.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
        report_params+="\n--------------------------------------------------------\n"

   output_tf.withWriter { w -> w << report_params}
}


//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        AmpliconArchitect-nf: Pipeline to predict ecDNA from WGS data (${workflow.manifest.version})
        """
}

//header for the IARC tools
// the logo was generated using the following page
// http://patorjk.com/software/taag  (ANSI logo generator)
def IARC_Header (){
     return  """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pipelines for cancer genomics.########################################
"""
}
