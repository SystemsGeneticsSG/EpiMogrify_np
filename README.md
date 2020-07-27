## Scripts used in this study

i) Proprietary codes: If you are interested in using EpiMogrify proprietary codes to prioritize cell identity, cell maintenance and cell conversion factors, please contact us at owen.rackham@duke-nus.edu.sg .  
ii) Non proprietary codes: Please see below for the non-proprietary codes which contains scripts for generating cell-type H3K4me3 profile from ENCODE data. 

### I. generateH3K4me3prolifes_ENCODE.py:  

Contains the script to generate H3K4me3 peak profiles from ENCODE bam files and the steps are:  
1. Create analysis directory for each ENCODE experiment ID
2. Check the presence of control ChIP-seq for the H3K4me3 ChIP-seq files. Set parameters for MACS2 peak calling software according to the control ChIP-seq information.
3. Check if the control sample is already mapped to hg38 genome orelse remap the control bam files to hg38 genome.
4. Check if H3K4me3 treatment sample bam file is derived from a single-end or paired-end sequencing and remap the bam file to hg38 genome.
5. Use MACS2 peak caller with broad parameter to generate broad H3K4me3 profiles.
6. Annotate the ChIP-seq peak to the nearest Ensemble transcripts.

#### List of softwares used
samtools  
bwa  
run_spp.R (https://github.com/kundajelab/phantompeakqualtools)  
macs2  
bedtools  


```
## Assign paths to directories
path_to_runSPP=/path_to_runSPP.R    ## Path to phantompeakqualtools directory that contains runSPP.R 
export PATH=$PATH:$path_to_runSPP

maindir=/path_to_H3K4me3Analysis
download_meta=$maindir/metadata_ChipSeq_27Feb2017.csv          # metadata downloaded from ENCODE on 27Feb2017
     
# For each ENCODE experiment ID, the metadata file is parsed and the treatment H3K4me3 sample is assigned to the control sample 
ChipSeqHistonMetafile=$maindir/metadata_ChipSeq_27Feb2017all-matchingcontrols.csv 

reference=/path_to_ref/GRCh38.primary_assembly.genome.fa       # Hg38 version v25 reference genome
annotation=/path_to_ref/gencode.v25.annotation.gtf


PeakType="Broad"                                                # Select Narrow/Broad. In our study we used "Broad" parameter in MACS. 
wrkdir=$maindir/Analysis_H3K4me3_Encode                         # directory where the results containing H3K4me3 ChIP-seq peak files will be stored 
exp_dir=$processdir/ChipSeq_Encode/H3K4me3Encode                # directory containing the downloaded bam files of H3K4me3 ChIP-seq from ENCODE
control_dir=$processdir/ChipSeq_Encode/ControlEncode            # directory containing the downloaded bam files of the control ChIP-seq from ENCODE


## H3K4me3_ENCODE_IDs_2017.txt contains the list of ENCODE experiments ids used in this study  
## For a list of ENCODE experiment IDs as an example. 
for experiment_id in ["ENCSR662PLB", "ENCSR000DTK", "ENCSR000DQH", "ENCSR000DUJ", "ENCSR000DWS" ]:
do  
   generateH3K4me3prolifes_ENCODE.py $ChipSeqHistonMetafile $experiment_id $reference $annotation $wrkdir $exp_dir $control_dir $PeakType
done
```


### II. celltype_H3K4me3profile.R  
Combines H3K4me3 profiles from different samples to generate cell-type H3K4me3 profile at the specified ChIP-seq q-value filter.  

R script:
```
## set directories
hmdir="path"

peakdir = paste0(hmdir, "/path_to_H3K4me3peakfiles/")
resdir = paste0(hmdir, "/path_to_Resultdir")

celltypes = c("primary_cell__astrocytes", "primary_cell__cardiac_muscle_cell", "stem_cell__H9")
qValue_cutoff = 5 ## -log10 pvalue (10^-5)

dir.create(resdir)


## Merge samples in each celltype. The union of broad peaks across samples are generated.
celltype_unionprofile(peakdir, resdir, celltypes, qValue_cutoff)
```


