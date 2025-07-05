This folder contains all the script used to analyse the Masoko dataset.
Analyses carried out:
1. GWAS
  a. egg_spot_Area
  b. egg_spot_size
  c. egg_spot_number
  d. egg_spot_number_Theis
  e. egg_spot_Hue
2. SELECTIVE SCANS
  a. FST
  b. DXY
  c. Tajima's D
  d. Pi
  e. xpEHH
3. Expresssion analyses using:
   a. RNAseq data from 58 wildcaught individuals
   b. RNAseq data from Aaron's experiment
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
1. GWAS  _Metdata [[ADD A LINK]] were subsetted based on the analyses to be carried out - NA's dropped and minimum sequence depth set to atleast 10X or above._
All GWAS were carried out using GEMMA and ANGSD: [see: https://github.com/Santos-cichlids/Methods-and-tools/tree/main/Bioinformatics/GWAS]
   a. egg_spot_Area
This analysis includes 288 male individuals, using the same input files and filtering criteria as the masoko_eggspot_size GWAS.

Location: /rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_eggspot_Area

Phenotype: #Aligned_Total_Eggspot_area_PropOfFin (from the metafile)

Filtering criteria:
Samples were filtered using the following conditions:

R
Copy
Edit
subset(masoko, Sex == "Male" & seq_depth > 10 & Aligned_SpotSize_PropOfFin != "NA" & SL != "NA")
Covariates used: Only SL (Standard Length) was included as a covariate in the GWAS.

b. egg_spot_size

Location: /rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_eggspot_size/

Shares the same filtering pipeline and input data as the egg_spot_Area analysis.
   c.
   d.
   e.
   

2. Selective scans
   a.
   b.
   c.
   d.
   e.
   f.
   
3. Expresssion analyses using:
   a. RNAseq data from 58 wildcaught individuals
   b. RNAseq data from Aaron's experiment
