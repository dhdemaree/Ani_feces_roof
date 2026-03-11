# Ani_feces_roof

Repository for the Journal article: systematic review of animal fecal pathogens that have the potential to affect quality of roof harvested rain water.  
Work conducted by David Demaree
Affiliation ORISE and US EPA
3/11/2026

This file contains all of the data that was used to create the figures for the journal article "Systematic review of animal fecal pathogens that have the potential to affect quality of roof harvested rain water".

# Figures
This folder contains all of the final figures that are going to be used for the journal article.

# Scripts
This file contains the scripts that were used to create the tables and figures for the metaanalysis

    Rawdata
      Tables of the outcomes that were drawn from the journal articles and are the inputs for the scripts
    
    Data_out
      Tables of data processed during the meta analysis
      
    Forest
      Images of the figures created by the scripts.
      
    HighDef
      High resolution images of the figures created by the scripts.

The Main Branch is where all the files that were directly held to creat the journal article.
RawData Branch unprocessed data.
Animal Taxonomy Branch processes for adding animal taxonomic information to unprocessed data.
GRADE Branch is where information on where GRADE and how GRADE was process and added to the unprocessed data.

N_Am_Meta_J-EPC0033278-QP-0_20230627_v03.csv  
This is the R script that was used to conduct the metaanalysis.

Inputs  
Campylobacter_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv  
Coliform_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv  
Cryptosporidium_STUDYN_DHD_J-EPC0033278-QP-1-0_220230627_v03.csv  
EColi_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv  
Enteriococci_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv  
Giardia_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv  
Salmonella_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv  

These are the input files that were used used as part of the meta analysis.  The files are CSV tables that contain all the information that is needed to conduct the meta-analysis.  The files are seperated by pathogen.  The file naming convention is pathogen, STUDYN, initial of file creator (DHD PI of project), QAPP#, date, and version.

Outputs  

Enumerated_box_display.csv  
This is the table that is displayed in the main article and gives a numeric representation of the data in Journal_EnumBoxplot.png.  

Enumerated_boxplot.csv  
This table is the same as the information displayed in Enumerated_box_display.csv, but also includes data that could not fit in the journal article and includes the standard error and the names of the references for the outcomes that were part of the sysnesized data.  

Total_Enumerated_Dataset.csv  
This table is the same as Enumerated_boxplot.csv except for is also includes data for reports that were from the US and from groups that were from 1 study so were excluded from the table in the journal article but were included in the total.  

Total_Prevalence_Dataset1.csv  
This table contains all the information that was in Fig. XX Journal_forestplot.png, and also included the fixed effect model, random effect model confidence limits, and the references for all the reports that the data came from.  

MetaanalysisReference.docx
This file contains the references for all the review articles that were part of the meta-analysis for this study.  It does not include references from the raw data.

Eligibility Criteria.docx  
This word document contains a table of all the elegibility criteria for reports for the systematic review.  In the table are a detailed explaination of each of the elegibility criteria as well as a the justification for why each of the elegibility criteria was used.  This is also in the Supplimentary Materials for the paper.

Grouping Pathogens.docx  
This word document contains a table with how the pathogens were grouped.  The first column Pathogen is the pathogen category that the pathogen falls under.  The thrid column is the name of the pathogen strain, species, serotype, etc.. from the report.  The second column is meta-analysis group which is the group that the retrieved pathogen names are grouped under.  The fourth column is the justification.    This is also in the Supplimentary Documents for the Paper.

ExampleMeta.docx  
This is a word document that contains a walk through of the processing of the data for birds and Escherichia coli from Fricker et al. 2018.  
