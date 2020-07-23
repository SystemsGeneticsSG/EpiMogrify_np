#!/usr/bin/python

"""Pipeline to realign all bam files """

import os, sys
import csv
import glob
import subprocess
from operator import itemgetter
import time
from datetime import datetime
import logging
import numpy as np


def celltypeBamSelect(ChipSeqHistonMetafile):
  # Select the non-redundant bam files to process along controls used
  exp_control = {}
  H3K4me3 = {}
  Control = {}
  with open(ChipSeqHistonMetafile, "r") as f:
    X = f.readlines()
    for x in X:
      tmp = csv.reader([x], delimiter=',').next()
      if tmp[2] == "H3K4me3-human":
        H3K4me3[tmp[0]] = tmp 
        bamfile_set =  tmp[4].replace("'", "").replace('[[', '').replace(']]', '').split("], [")
        if "GRCh38" in [i.split(", ")[2] for i in bamfile_set]:
          align_index = [i.split(", ") for i in bamfile_set if (i.split(", ")[0] == "alignments") and (i.split(", ")[2] == "GRCh38")]
        else: 
          align_index = [i.split(", ") for i in bamfile_set if (i.split(", ")[0] == "alignments") and (i.split(", ")[2] == "hg19")]
        exp_control[tmp[0]] = [[a[1] for a in align_index],  tmp[-1].replace("[", "").replace("]", "").replace("'", "").split(", ")] 
      elif tmp[2] == "Control-human":
        Control[tmp[0]] = tmp
  for e in exp_control:
    if exp_control[e][-1][0] not in Control: print exp_control[e][-1][0] , " Not found"
    else:
      #print  exp_control[e][-1][0]
      con_set = Control[exp_control[e][-1][0]][4].replace("'", "").replace('[[', '').replace(']]', '').split("], [")  # Assuming the control is 1 Experiment
      if "GRCh38" in [i.split(", ")[2] for i in con_set]:
        Calign_index = [i.split(", ") for i in con_set if (i.split(", ")[0] == "alignments") and (i.split(", ")[2] == "GRCh38")]
      else: 
        Calign_index = [i.split(", ") for i in con_set if (i.split(", ")[0] == "alignments") and (i.split(", ")[2] == "hg19")]
      exp_control[e][-1] = [a[1] for a in Calign_index]
  return(H3K4me3, exp_control)



def peakCallingPipeline(exp_details, exp_contol, reference, ref_annot, wrkdir, exp_dir, control_dir, PeakType):
  import time
  from datetime import datetime
  import logging
  start_time = datetime.now()
  
  if not os.path.exists(wrkdir + "/" + exp_details[1].replace(" ", "_")):
    os.makedirs(wrkdir + "/" + exp_details[1].replace(" ", "_"))
  
  os.chdir(wrkdir + "/" +  exp_details[1].replace(" ", "_"))
  
  # Check whether the experiment bam files have already been analysed:
  flag = "Analysed"
  print "flag: ", flag

  for expbamfilestocheck in exp_contol[exp_details[0]][0]:
    print expbamfilestocheck
    if not os.path.exists(expbamfilestocheck + "_peaks.broadPeak"):  
      flag = "NOT analysed"
      print flag       

  if flag != "NOT analysed":
    print "The Experiment has already been processed: ", exp_details[0]
  else:
  # Initiate logging"
    logger = logging.getLogger(exp_details[1].replace(" ", "_"))
    logger.setLevel(logging.DEBUG)

    fh = logging.FileHandler(wrkdir + "/" + exp_details[1].replace(" ", "_") +'/Run.log')
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s ... %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    logger.info("\n### 1. Create analysis directory")
    logger.info("Analysis directory created: " + wrkdir + "/" + exp_details[1].replace(" ", "_"))
    process_log = open( wrkdir + "/" + exp_details[1].replace(" ", "_") + "/Process_outerr.log" , "a")
  
    logger.info("\n\nPeak calling Pipeline starts")

    logger.info("\n\n### 2. Check the number of Control files and set parameters for MACS2" + str(exp_contol[exp_details[0]][1]))
    if len(exp_contol[exp_details[0]][1]) == 1: 
      controlfile = exp_contol[exp_details[0]][1][0] + ".bam"    # for Bwa
      macs_controlargs =  " -c " + control_dir + "/" + controlfile.replace(".bam", "_hg38.bam")  # For MACS

    elif len(exp_contol[exp_details[0]][1]) > 1: 
      all_bamfiles = ""
      for c in exp_contol[exp_details[0]][1]:
        all_bamfiles = all_bamfiles + " " + control_dir + "/" + c + ".bam"
      controlfile =  "mergeControl_" + exp_details[0] + ".bam"
      macs_controlargs = " -c " + control_dir + "/" + controlfile.replace(".bam", "_hg38.bam") 
      logger.info("Multiple control files, Merging : " + str(exp_contol[exp_details[0]][1]))
      mergecmd = "samtools merge " + control_dir + "/" + controlfile + " " + all_bamfiles
      logger.info("Merge : " + mergecmd)
      subprocess.call(mergecmd, shell=True, stdout=process_log, stderr=subprocess.STDOUT)
    elif len(exp_contol[exp_details[0]][1]) == 0:
      macs_controlargs = " --nolambda "
    logger.info("MACS2 control parameters: " + macs_controlargs)

    logger.info("\n\n### 3. Check if the control sample is already mapped to hg38 genome orelse remap the bam files to hg38")
    if not os.path.exists(control_dir + "/" + controlfile.replace(".bam", "_hg38.bam")) and macs_controlargs != " --nolambda ":
      flag = subprocess.check_output("samtools view -c -f 1 " + control_dir + "/" + controlfile, shell=True).strip()
      if flag == "0": # Single end
        logger.info("Bam file : Single-end sequencing")
        controlMapcmd = "samtools bam2fq " + control_dir + "/" + controlfile + \
      " | bwa mem " + reference + " - | samtools view -Sb - > " + control_dir + "/" + controlfile.replace(".bam", "_hg38.bam")
      else: #Paired-end
        logger.info("Bam file : Paired-end sequencing")
        controlMapcmd = "samtools bam2fq " + control_dir + "/" + controlfile + \
      " | bwa mem -p " + reference + " -  | samtools view -Sb - > " + control_dir + "/" + controlfile.replace(".bam", "_hg38.bam")
      logger.info(controlMapcmd)
      subprocess.call(controlMapcmd, shell=True, stdout=process_log, stderr=subprocess.STDOUT)  
    else: 
      logger.info("Control file present: " + controlfile)

    logger.info("\n\n### 4. Check if Experiment sample Bam derived is from Single or paired and Remap to hg38 v25")

    for b in exp_contol[exp_details[0]][0]:
      logger.info("\n\n *** Processing: " + b + " out of " + str(exp_contol[exp_details[0]][0]))
      bamfile = b + ".bam"
      flag = subprocess.check_output("samtools view -c -f 1 " + exp_dir + "/" + bamfile, shell=True).strip()
    # Convert Experiment BAM to fastq
      if flag == "0": # Single end
        logger.info("Bam file : Single-end sequencing")
        mappBamcmd = "samtools bam2fq " + exp_dir + "/" + bamfile + \
        " | bwa mem " + reference + " - | samtools view -Sb - > " +  bamfile.replace(".bam", "_hg38.bam")
        macs_seq = " -f BAM  "
      else: #Paired-end
        logger.info("Bam file : Paired-end sequencing")
        mappBamcmd = "samtools bam2fq " + exp_dir + "/" + bamfile + \
          " | bwa mem -p " + reference + " - | samtools view -Sb - > " + bamfile.replace(".bam", "_hg38.bam")
        macs_seq = " -f BAMPE  "
      logger.info(mappBamcmd)
      if not os.path.isfile(bamfile.replace(".bam", "_hg38.bam")):
        subprocess.call(mappBamcmd, shell=True, stdout=process_log, stderr=subprocess.STDOUT)
      else: logger.info("Bam file present: " + bamfile )
  
      logger.info("\n\n### 5. MACS2 peak caller: Broad or narrow")
      # Find the fragment size of the bam file....  
      qcBampeak =  "R CMD BATCH  --no-save --no-restore '--args -c=" + bamfile.replace(".bam", "_hg38.bam") + \
	" -savp=" + bamfile.replace(".bam", "_peak.pdf") +" -out=" + bamfile.replace(".bam", "_peak.stats") + " ' run_spp.R "

      logger.info( qcBampeak)
      subprocess.call(qcBampeak, shell=True, stdout=process_log, stderr=subprocess.STDOUT) 
  
      if PeakType == "Broad":
        with open( bamfile.replace(".bam", "_peak.stats"), 'r') as f:
          stat = f.readline()
	  if len(stat.split("\t")[2].split(",")) == 1: fragmentsiz = stat.split("\t")[2]  
          else: 
	    fragmentsiz =  stat.split("\t")[2].split(",")[0] # If multiple: the majority has first frag size
        print stat.split("\t")[2]
        peakCallcmd = "macs2 callpeak -B --broad --broad-cutoff 0.1 " + macs_seq + " -t " + bamfile.replace(".bam", "_hg38.bam") + " "  \
          + macs_controlargs + " -n " + bamfile.replace(".bam", "") + " --nomodel --extsize "  + fragmentsiz
      elif PeakType == "Narrow":
        peakCallcmd ="macs2 callpeak -B -q 0.01 " + macs_seq + " -t " + bamfile.replace(".bam", "_hg38.bam") + " " + macs_controlargs \
          + " -n " + bamfile.replace(".bam", "")
      logger.info( peakCallcmd)
      subprocess.call(peakCallcmd, shell=True, stdout=process_log, stderr=subprocess.STDOUT)

      logger.info("\n\n### 6. Annotate the peak to the nearest Ensemble transcripts")
      if PeakType == "Broad":
        annotPeak = "bedtools closest -d -b  " + ref_annot + " -a " + bamfile.replace(".bam", "_peaks.broadPeak")  + \
          ' -D  "b" >  ' + bamfile.replace(".bam", "_peaks.broadPeakAnnot")
      elif PeakType == "Narrow":
        annotPeak = "bedtools closest -d -b " + ref_annot + " -a " + bamfile.replace(".bam", "_peaks.narrowPeak")  + \
          ' -D  "b" >  ' + bamfile.replace(".bam", "_peaks.narrowPeakAnnot")
      logger.info( annotPeak)
      subprocess.call(annotPeak, shell=True, stdout=process_log, stderr=subprocess.STDOUT)
      covCmd = "bedtools coverage -abam " + bamfile.replace(".bam", "_hg38.bam") + " -b " + \
  	  bamfile.replace(".bam", "_peaks.broadPeak") + " -counts > " + bamfile.replace(".bam", ".depth") 
      logger.info(covCmd)
      subprocess.call(covCmd, shell=True, stdout=process_log, stderr=subprocess.STDOUT)

    end_time = datetime.now()
    logger.info('\n\nCompleted Reprocessing..... Duration: {}'.format(end_time - start_time)) 
  

if len(sys.argv) == 9:
  [H3K4me3, exp_control] = celltypeBamSelect(sys.argv[1])
  peakCallingPipeline(H3K4me3[sys.argv[2]], exp_control, sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
else:
  print sys.argv  
