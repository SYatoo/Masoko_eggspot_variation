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
### 1. GWAS  
_Metdata (https://docs.google.com/spreadsheets/d/1rJda7cw5H4N52hMqD5NrAQfIXbaVSL6f/edit?gid=1587235546#gid=1587235546) were subsetted based on the analyses to be carried out - NA's dropped and minimum sequence depth set to atleast 10X or above._
All GWAS were carried out using GEMMA and ANGSD: [see: https://github.com/Santos-cichlids/Methods-and-tools/tree/main/Bioinformatics/GWAS]
#### a. `egg_spot_Area`

This analysis includes **288 male individuals**.

- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_eggspot_Area`

- **Phenotype column** (from `metafile`):  
  `Aligned_Total_Eggspot_area_PropOfFin`

- **Filtering criteria**:  
  ```r
  subset(masoko, 
         Sex == "Male" & 
         seq_depth > 10 & 
         Aligned_SpotSize_PropOfFin != "NA" & 
         SL != "NA")
   
#### b. `egg_spot_size`
This analysis includes **288 male individuals**.

- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_eggspot_size`

- **Phenotype column** (from `metafile`):  
  `Aligned_SpotSize_PropOfFin`

- **Filtering criteria**:  
  ```r
  subset(masoko, 
         Sex == "Male" & 
         seq_depth > 10 & 
         Aligned_SpotSize_PropOfFin != "NA" & 
         SL != "NA")

#### c. `egg_spot_number`
This analysis includes **260 male individuals**.

- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_GWAS_eggsot_number`

- **Phenotype column** (from `metafile`):  
  `EggspotNumber`

- **Filtering criteria**:  
  ```r
  subset(masoko, 
         Sex == "Male" & 
         seq_depth > 10 & 
         EggspotNumber != "NA" & 
         SL != "NA")
#### d. `egg_spot_number_Theis`
This analysis includes **243 male individuals**.

- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_GWAS_eggSpot_theis`

- **Phenotype column** (from `metafile`):  
  `EggspotNumberTheis`

- **Filtering criteria**:  
  ```r
  subset(masoko, 
         Sex == "Male" & seq_depth > 10 &
  EggspotNumberTheis != "NA" & 
         SL != "NA")
 
#### e. `egg_spot_number_Hue`
This analysis includes **160 male individuals**.

- **Directory**:  
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/GWAS/masoko_HUE`

- **Phenotype column** (from `metafile`):  
  `Aligned_SpotSize_PropOfFin`

- **Filtering criteria**:  
  ```r
  subset(masoko, 
         Sex == "Male" & 
         seq_depth > 10 & 
         spot.hue != "NA" & 
         SL != "NA")

2. Selective scans
   Filtering:...............
#### a. Fst 
#### b. dxy
#### c. Tajima's D 
#### d. Pi 
#### e. xp-EHH:  
_(requires phased data - data were phased to fAstCal1.2 reference)_
This scan was carried out on 187+98 (littoral+benthic respectively).
`Calc_IHH_xpEHH.R` was used which outputs IHH for each population, in addition to xp-EHH.
- **Directory**:
  `/rds/project/rds-8b3VcZwY7rY/projects/cichlid/suhaib/masoko/genome_scans_ALL/BENTHIC_LITTORAL/EHH_IHH_Selection/biallelics/HIGH_SEQ_depth/`

3. Expresssion analyses using:
   a. RNAseq data from 58 wildcaught individuals
#### b. RNAseq data from Aaron's experiment:
These data were analysed and plotted using `Domain_Sex_deseq_Aaron.R`
   
