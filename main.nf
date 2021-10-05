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
      cnvkit.py batch segment -p 1 -o ${tumor_id}/${cram.baseName}.segment.cns ${tumor_id}/${cram.baseName}.cnr
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
      echo cnvkit.py batch segment -p 1 -o ${tumor_id}/${cram.baseName}.segment.cns ${tumor_id}/${cram.baseName}.cnr
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
        export AA_DATA_REPO=${aa_repo_dir_path}
      	python2 /home/programs/AmpliconArchitect-master/src/amplified_intervals.py \\
                  --gain 4.5 --cnsize_min 50000 --ref GRCh38 \\
                  --bed ${cnvkit_bed}  --out ${cnvkit_bed.baseName}.aa.bed


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


//we run the amplicon architect tool




// // Annnot the VCF with VEP tools
// //create a local VEP database (gencode 33) ~ 16Gb size
// //vep_install -a cf -s homo_sapiens -y GRCh38 -c vep-db-99 --CONVERT
// //https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/vep.html
// process VEP {
//  cpus params.cpu
//  memory params.mem+'G'
//
//   publishDir params.output_folder+'/VEP/', mode: 'copy'
//   input:
//   set val(tumor_id), file(vcf), file(normal), file(normal_index), val(normal_id), val(tumor_id_name) from xVEP_input
//   file (vep_dir_path)
//   output:
//     set val(tumor_id), file("${tumor_id}.vep.vcf") into xVEP_out
//   script:
//        """
//         vep -i ${vcf} \\
//         -o ${tumor_id}.vep.vcf \\
//         --cache --offline \\
//         --dir_cache ${vep_dir_path} \\
//         --format vcf \\
//         --vcf \\
//         --symbol  \\
//         --terms SO \\
//         --tsl \\
//         --hgvs \\
//         --fasta ${vep_dir_path}/homo_sapiens/99_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \\
//         --plugin Frameshift \\
//         --plugin Wildtype \\
//         --dir_plugins ${baseDir}/VEP_plugins \\
//         --pick  --transcript_version
// 	# we remove the VAF from the VCF
// 	bcftools annotate -x  FORMAT/AF ${tumor_id}.vep.vcf > ${tumor_id}.vep.noAF.vcf
// 	mv ${tumor_id}.vep.noAF.vcf ${tumor_id}.vep.vcf
//        """
// }
//
// //we have to sync the xHLA, VEP ouputs and pvac input
// xHLA_xVEP=xHLA_out.join(xVEP_out, remainder: true)
// pvac_hla_vep=pvactools_input.join(xHLA_xVEP,remainder: true)
//
// process pVactools {
//  cpus params.cpu
//  memory params.mem+'G'
//
//   publishDir params.output_folder+'/pVACTOOLS/', mode: 'copy'
//   input:
//   set val(tumor_id), file(vcf), file(normal), file(normal_index), val(normal_id), val(tumor_id_name), file(hla_dir_out),file(vcf_vep) from pvac_hla_vep
//
//   output:
//     set val(tumor_id), path("${tumor_id}*_pvactools") into pVACTOOLS_out
//     file("${tumor_id}.pvactools.log")
//   script:
//        """
//          perl ${baseDir}/scripts/pbactools_wrapper.pl -a ${hla_dir_out} \\
//             -b ${baseDir}/db/xHLA2PVAC_alleles.txt -c ${normal_id}   -d ${vcf_vep} -t ${tumor_id_name} -p ${tumor_id} \\
//             -e ${params.pvactools_predictors} > ${tumor_id}.pvactools.log
//        """
// }


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
  def software_versions = ['cnvkit' : '0.8.9',
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
