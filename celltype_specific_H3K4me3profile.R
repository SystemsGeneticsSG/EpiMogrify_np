
## load libraries
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(grid)
require(gridExtra)
library(data.table)
options(stringsAsFactors = FALSE)


## Merge samples in each celltype. The union of broad peaks across samples are generated.
celltype_unionprofile <- function(peakdir, resdir, celltypes, qValue_cutoff ){
  extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric")
  for (ct in celltypes){
    files = list.files(peakdir, ct)
    i=0
    if (length(files) != 0){
      x = GRanges()
      cat ("Processing cell type ", ct, "\n" )
      for (f in files){
        cat ("Sample: ", f[1])
        pks = import.bed(paste0(peakdir,"/", f[1]), extraCols = extraCols_broadPeak)
        pks_cutoff = pks[elementMetadata(pks)$qValue  >= qValue_cutoff]
        cat("\t Number of peaks: ", length(pks_cutoff), "\n" )
        if (length(pks_cutoff) >= 10000){
          x = c(x, pks_cutoff) 
        }
      }
      if (length(x) !=0){
        x = reduce(x)
        export.bed(x, paste0(resdir, "/",  ct, "-union.bed"))
      }
    }
  }
}

## set directories
hmdir="path"

peakdir = paste0(hmdir, "/path_to_H3K4me3peakfiles/")
resdir = paste0(hmdir, "/path_to_Resultdir")

celltypes = c("primary_cell__astrocytes", "primary_cell__cardiac_muscle_cell", "stem_cell__H9")
qValue_cutoff = 5 ## -log10 pvalue = 10^-5

dir.create(resdir)


## Merge samples in each celltype. The union of broad peaks across samples are generated.
celltype_unionprofile(peakdir, resdir, celltypes, qValue_cutoff)

