This is version 3:
- use new adduct information and retention times for *standards* defined by Mar.
- use `integrate = 2` and `peakwidth = c(2, 20)` for the peak detection.)

# VAMS vs intraveinous

Comparison of capillary whole blood Volumetric Absorbtive Microsampling to
intravenous whole blood sampling to intravenous plasma to intravenous red blood
cells (RBC)

All data files as well as supporting information files (at the moment only
injection sequence) are on the biobank nas under the folder 'data/vams2018'

In the same folder there is also files that contain the injection order in a
sub-folder called 'csv'.

The file names follow a naming format convention although the pool samples break
the format to a lesser degree.  The format is:
date-of-analysis_sampleid_matrixid_replicate_polarity.  SampleID is just a
running number, matrixID is one of a few values: '1', '2', '3', 'RBC' and
'POOLxx'.  The xx in the 'POOLxx' indicates the injection order of the pooled
sample, while 1, 2, and 3 correspond to 'plasma', 'intravenous whole blood', and
'capillary whole blood' respectively. RBC and POOL should be self explanatory.

Note that most of the files other than the peack picking and feature
identification rely on R objects that are generated in the feature
identification script. If you haven't run that then chances are the other ones
will not work.


## Files and structure of the analysis

Data preprocessing and normalization:

- [vams_preprocessing.Rmd](vams_preprocessing.Rmd): Defines the pre-processing
  of the positive polarity data (chromatographic peak detection, alignment and
  correspondence).
- [vams_preprocessing_neg.Rmd](vams_preprocessing_neg.Rmd): pre-processing of
  the negative polarity data.
- [vams_normalization.Rmd](vams_normalization.Rmd): performs the normalization
  of the feature abundances (positive polarity).
- [vams_normalziation_neg.Rmd](vams_normalization_neg.Rmd): normalization of the
  negative polarity data.

Differences between matrices:

- [vams_matrix_pos.Rmd](vams_matrix_pos.Rmd): analysis of the matrix specific
  effects (positive polarity).
- [vams_matrix_neg.Rmd](vams_matrix_neg.Rmd): matrix specific effects for
  negative polarity data.
- [vams_matrix_pos_semtar.Rmd](vams_matrix_pos_semtar.Rmd): *semi-targeted*
  analysis for the matrix specific effects (positive polarity). The analysis is
  restricted to features potentially measuring ~ 120 *standards* that were
  measured previously and for which approximate retention times are available.
- [vams_matrix_neg_semtar.Rmd](vams_matrix_neg_semtar.Rmd): *semi-targeted*
  analysis for the matrix specific effects (negative polarity).

Differences between male and female individuals:

- [vams_sex_pos.Rmd](vams_sex_pos.Rmd): analysis to identify features with
  significant abundances between male and female participants.
- [vams_sex_neg.Rmd](vams_sex_neg.Rmd): analysis to identify features with
  significant abundances between male and female participants, negative polarity
  data.
- [vams_sex_pos_semtar.Rmd](vams_sex_pos_semtar.Rmd): *semi-targeted* analysis
  for the sex differences, positive polarity.
- [vams_sex_neg_semtar.Rmd](vams_sex_neg_semtar.Rmd): *semi-targeted* analysis
  for the sex differences, negative polarity.

## Required packages and setup

The analysis requires a recent version of R (>= 3.6.0) and the following R
packages. All required R packages can be installed with the code below.

```r
install.packages("BiocManager")
BiocManager::install(c("BiocStyle",
                       "xcms",
                       "RColorBrewer",
                       "pander",
                       "doParallel",
                       "magrittr",
                       "pheatmap",
                       "DESeq2",
                       "edgeR",
                       "NormalizeMets",
                       "MetNorm",
                       "ruv",
                       "SummarizedExperiment",
                       "UpSetR",
                       "limma",
                       "writexl",
                       "devtools"))
#' Install from github if versions are "too old"
if (packageVersion("MSnbase") < "2.9.3")
    devtools::install_github("lgatto/MSnbase")
if (packageVersion("xcms") < "3.5.2")
    devtools::install_github("sneumann/xcms")
```

## Raw mzML data location

mzML files with the raw spectrum data are available from the calculation
clusters in */data/massspec/mzML/*.

