# Vaccination_VoC_Publication

R code related to the publication in Nature Communications: "Immune response to SARS-CoV-2 variants of concern in vaccinated individuals"

The code in the file "Analysis_submission.R" used to generate the graphs used in the paper. In general graphs are exported to an output folder in .svg format, from where they are further processed into the finalized figures using a vector graphics editing software such as InkScape.

R version used was version 3.6.1 (2019-07-05) -- "Action of the Toes" within RStudio 1.2.5001

Required are an input folder at "/input" with the following files containing annotated measurement data for the different assays employed:
"ACE2_Results.csv"
"ELISA_Saliva_Results.csv"
"MULTICOV-AB_Results.csv"
"MULTICOV-AB_Saliva_Results.csv"
"NBP_Results.csv"
"VNT_Results.csv"
"DilSeries.csv"
