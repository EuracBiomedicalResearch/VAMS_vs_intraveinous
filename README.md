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

- [vams_preprocessing.Rmd](vams_preprocessing.Rmd): Defines the pre-processing
  of the data (chromatographic peak detection, alignment and correspondence).
- [vams_normalization.Rmd](vams_normalization.Rmd): performs the normalization
  of the feature abundances.


