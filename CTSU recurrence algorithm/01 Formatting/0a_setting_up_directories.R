# ********************* CTSU Recurrence Algorithm  ************************


# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)

# Date: 15/05/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# ------------------------------------------------------------------------------
# Cancer summary preparation - read in frozen version May 2023
library(dplyr)
library(tidyverse)
library(magrittr)
library(stats)
library(lubridate)
library(ggplot2)
library(janitor)
library(stringr)
library(summarytools)
library(scales)
library(gmodels)
library(flextable)
library(officer)
library(writexl)
library(gtsummary)
library(gt)
library(flextable)

# ------------------------------------------------------------------------------
# Set working directory and load in treatment datasets

# read in cohort code

# cohort in dataframe BC_dat - with pseudonymised id (BGS_ID) and primary 
# breast cancer data

setwd()

# read in following datasets
hes_APC <- read.csv()
hes_APC

hes_OP <- read.csv()
hes_OP

# dataset with ICD code and corresponding ICD description
lines <- readLines("ICD10codes3char.txt")
icd10 <- data.frame(text = lines, stringsAsFactors = FALSE)
icd10

# dataset with OPCS-4 code and corresponding OPCS-4 description
lines2 <- readLines("OPCS410codes.txt")
opcs4 <- data.frame(text = lines2, stringsAsFactors = FALSE)
opcs4

sact <- read.csv()
sact

rtds <- read.csv()
rtds

# vital status dataset
end_status <- read.csv()
end_status

# cancer registry patient + tumour characteristic dataset
tumour_phe <- read.csv()
tumour_phe

# used when considering treatment extended in her2 positive
wks <- 64
her2wks <- 26

# cancer registry treatment dataset
treatment <- read.csv()
treatment

############################
#### END OF CODE BLOCK #####
############################