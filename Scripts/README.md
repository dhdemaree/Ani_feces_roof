# Scripts

This file contains the R scripts that were used for the meta analysis and synthesis of the outcomes drawn from the systematic review.

The scripts in this folder and the folders in the folder are set up in the way that they would be used and created in the Rscripts.  So if the whole Scripts file is downloaded you should be able to run N_Am_Meta_DHD_J-EPC0033278-QP-1-0_20260311_v04.R and HighDef_Figures_DHD_J-EPC0033278-QP-1-0_20260226_v01.R. and it should use the data from the Rawdata and it should output the files into Data_out, HighDef, and Forest.


# R Scripts

N_Am_Meta_DHD_J-EPC0033278-QP-1-0_20260311_v04.R
This file is the main R script that conducts the meta analysis for the journal article.  The input files are from the Rawdata and are the outcomes from the journal articles.  This script processes the outcomes from Rawdata to create Data_out which is the text based data that was created during the meta analysis and Forest which contains the images of the figures that were created while processing the data. This script requires Rstudio to set the working directory.   

HighDef_Figures_DHD_J-EPC0033278-QP-1-0_20260226_v01.R
This R scipt makes Figures in the resolution for the target journal Water Reasearch.  It reqires that N_Am_Meta_DHD_J-EPC0033278-QP-1-0_20260311_v04.R has been run first or that the workspace that was created in that script is loaded.  This script will then make figured based on the data in the workspace.   Though this script is mostly from N_Am_Meta_DHD_J-EPC0033278-QP-1-0_20260311_v04.R it is being kept seperate for the purpose of creating figures. The ouput of this script is placed in HighDef.

# Input

Rawdata
This file contains cvs files that have all the outcomes that were drawn from the reports.  Each file is seperated by pathogen.

# Output

Data_Out
This file contains all of the text based data that was generated from the script N_Am_Meta_DHD_J-EPC0033278-QP-1-0_20260311_v04.R.

Forest
This file contains all the images of the figures that were created by the script N_Am_Meta_DHD_J-EPC0033278-QP-1-0_20260311_v04.R.

HighDef
This file contains the same images as Forest but the images are higher resolution and were created by the script HighDef_Figures_DHD_J-EPC0033278-QP-1-0_20260226_v01.R.
