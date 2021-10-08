# ampliconarchitect-nf
Nextflow pipeline to discover ecDNA in cancer genomes.


## Usage
```
  #using a tn_pairs file
  nextflow run iarcbioinfo/ampliconarchitect-nf -r v1.0 \ 
  -profile singularity  --tumor_file tumor_file.txt \
  --cram_dir $PWD/CRAM \
  --moseklm_license mosek.lic --aa_repo_dir GRCh38 \
  --output_folder ecDNA_public1
  
  #testing a small fake run
   nextflow run iarcbioinfo/ampliconarchitect-nf -r v1.0 \ 
   --tumor_file test/tumor_file.txt \
   --cram_dir test/cram \
   --moseklm_license test/mosic_license/mosic.license \
   --aa_repo_dir test/amplicon_repo --debug true
   ```
 
 ## Dependencies

###  License for Mosek optimization tool

Obtain license file mosek.lic (https://www.mosek.com/products/academic-licenses/ or https://www.mosek.com/try/).

### Amplicon architect repository
To run amplicon architect is necessary to download the Amplicon Architect repository from [here](https://drive.google.com/drive/folders/0ByYcg0axX7udeGFNVWtaUmxrOFk), currently we have tested this nextflow pipeline with GRCh38 genome built.

You can avoid installing all the external software by only installing Docker or singularity.
See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.

## Input 
### Mandatory
  | Parameter      | Description   |
  |-----------|---------------|
  | --cram_dir   | Directory containing all BAM/CRAM files |  
  | --tumor_file    | File containing list of tumor CRAM files to be processed |
  | --moseklm_license       |  License for mosek program [mosic.license]|
  | --aa_repo_dir | directory containing the Amplicon Architect repository [See above]|
  
### Optional
  | Parameter      | Description   |
  |-----------|---------------|
  | --output_folder   | name of output folder [def:amplicon_results]|  
  | --cpu  | Number of CPUs[def:1] |
  |  --mem     |  Max memory [def:8Gb]|
  | --debug | enable debuging of the pipeline|
  
### Example of Tumor file (--tumor_file)
A text file tabular separated, with the following header:
```
id	cram
sample1_T1	sample1_T.cram
sample2_T1	sample2_T.cram
sample3_T1	sample3_T.cram
``` 


## Output

```
├── amplicon_classes                # amplicon classified
│   └── all_amplicons_classification_bed_files # amplicon in genomic coordiantes 
├── amplicon_predictions           # Amplicon architect predictions
│   ├── S00016_T.aa
│   ├── S01493_T.aa
│   ├── S01501_T.aa
├── amplified_intervals          #Selected seeds per sample
├── cnvkit							#CNVkit copyNumber segments
│   ├── S00016_T
│   ├── S01493_T
│   ├── S01501_T
│   ├── S01502_T
└── nf-pipeline_info           # Nextflow logs
```

## Common errors

### Singularity
The first time that the container is built from the docker image, the TMPDIR  should be defined in a non parallel file-system, you can set this like:

```
export TMPDIR=/tmp
```

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Matthieu Foll*    |            follm@iarc.fr | Developer to contact for support (link to specific gitter chatroom) |
  | Alex Di Genova | digenovaa@fellows.iarc.fr| Developer |