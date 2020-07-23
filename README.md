# Scripts to generate H3K4me3 peak profiles from ENCODE bam files  


1. Create analysis directory
2. Check the number of Control files and set parameters for MACS2
3. Check if the control sample is already mapped to hg38 genome orelse remap the bam files to hg38
4. Check if Experiment sample Bam derived is from Single or paired and Remap to hg38 v25
5. MACS2 peak caller: Broad or narrow
6. Annotate the peak to the nearest Ensemble transcripts


```
## Assign paths to directories
path_to_runSPP=/path_to_runSPP.R    ## Path to phantompeakqualtools directory that contains runSPP.R 
export PATH=$PATH:$path_to_runSPP

maindir=/path_to_H3K4me3Analysis
download_meta=$maindir/metadata_ChipSeq_27Feb2017.csv
ChipSeqHistonMetafile=$maindir/metadata_ChipSeq_27Feb2017all-matchingcontrols.csv

reference=/path_to_ref/GRCh38.primary_assembly.genome.fa
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


