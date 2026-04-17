# ********************* CTSU Recurrence Algorithm  ************************

# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 04/06/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# ------------------------------------------------------------------------------

# merge cancer pathway mainspef and tretspef codes to the HES OP and HES APC

# ------------------------------------------------------------------------------
# STATA Script 6a. HES main trt speciality all cancers
# Purpose: Coding HES speciality codes
# ------------------------------------------------------------------------------

hes_OP$appdate_hesop <- hes_OP$apptdate

hes_OP <- hes_OP %>%
  arrange(bgs_id, appdate_hesop) %>%
  group_by(bgs_id) %>% mutate(OPid_n = row_number()) %>% ungroup()

hes_OP$mainspef_hesop <- hes_OP$mainspef
hes_OP$tretspef_hesop <- hes_OP$tretspef

#Create medical oncology flags and radiology/radiotherapy flags
hes_OP <- hes_OP %>%
 mutate(HEStrtoncol_flag = case_when(
   tretspef == "800" ~ 1, # Clinical Oncology (previously known as radiotherapy)
   tretspef == "370" ~ 1, # Medical Oncology
   tretspef == "503" ~ 1, # Gynaecological Oncology
   tretspef == "811" ~ 1, # Interventional Radiology
   tretspef == "103" ~ 2, # Breast Surgery (Includes suspected neoplasms, cysts etc, does not include cosmetic)
   TRUE ~ NA),
   HESmainoncol_flag = case_when(
     mainspef == "800" ~ 1,
     mainspef == "370" ~ 1,
     mainspef == "810" ~ 1, # Radiology
     mainspef == "103" ~ 2,
     TRUE ~ NA
     ))
#NOTE: paedatric medical oncology not included as not relevant for this data

# ------------------------------------------------------------------------------
# HES APC

hes_APC$admidate_hesapc <- hes_APC$admidate

hes_APC <- hes_APC %>%
  arrange(bgs_id, admidate_hesapc) %>%
  group_by(bgs_id) %>% mutate(APCid_n = row_number()) %>% ungroup()

hes_APC$mainspef_hesapc <- hes_APC$mainspef
hes_APC$tretspef_hesapc <- hes_APC$tretspef

# Create medical oncology flags (oncology and RT)
hes_APC <- hes_APC %>%
  mutate(HEStrtoncol_flag = case_when(
    tretspef == "800" ~ 1, # Clinical Oncology (previously known as radiotherapy)
    tretspef == "370" ~ 1, # Medical Oncology
    tretspef == "503" ~ 1, # Gynaecological Oncology
    tretspef == "811" ~ 1, # Interventional Radiology
    tretspef == "103" ~ 2, # Breast Surgery (Includes suspected neoplasms, cysts etc, does not include cosmetic)
    TRUE ~ NA),
    HESmainoncol_flag = case_when(
      mainspef == "800" ~ 1,
      mainspef == "370" ~ 1,
      mainspef == "810" ~ 1, # Radiology
      mainspef == "103" ~ 2,
      TRUE ~ NA
      ))

# ------------------------------------------------------------------------------
# NEW STATA SCRIPT
# Script: 6e. POETIC HES long format
# Purpose: Put HES data into long format to code
# ------------------------------------------------------------------------------
# 1 Long format diagnosis

hes_OP$temp <- "_OP"
hes_OP <- transform(hes_OP, HESid_n = paste0(OPid_n, temp))
hes_OP <- hes_OP %>% select(-temp)

head(hes_OP)

hes_OP_diag <- hes_OP

# check diag is character not numeric and only use columns with data in
summary_diagvars <- hes_OP_diag %>%
  select(starts_with("diag"))
summary(summary_diagvars)

# list of variables to process
varlist <- c('diag_01', 'diag_02', 'diag_03', 'diag_04', 'diag_05',
             'diag_06', 'diag_07', 'diag_08', 'diag_09', 'diag_10',
             'diag_11', 'diag_12')

#for (v in varlist) {
  #stopifnot(all(is.na(hes_OP_diag[[v]]))) # will give error if false
#}


# only keep diagnosis data with some data included
hes_OP_diag <- hes_OP_diag %>% select(bgs_id, HESid_n, diag_01:diag_12)

# rename the diagnosis codes so can put it into long format
#initialize the dummy variable
varlist <- c('diag_01', 'diag_02', 'diag_03', 'diag_04', 'diag_05',
             'diag_06', 'diag_07', 'diag_08', 'diag_09', 'diag_10',
             'diag_11', 'diag_12')

hes_OP_diag$dum <- 0

for (i in seq_along(varlist)) {
  
  hes_OP_diag$dum <- hes_OP_diag$dum + 1 # increment the dummy variable
  
  d <- hes_OP_diag$dum[i] # display the current value of dum
  n <- substr(varlist[i], 1, 5) # extract and modify the variable name
  
  new_name <- paste0(n, d)
  
  hes_OP_diag <- hes_OP_diag %>% rename(!!new_name := !!sym(varlist[i]))
  
}

head(hes_OP_diag)

# Reshape the diagnosis data in HES OP into a long format
hes_OP_diag_long <- hes_OP_diag %>%
  pivot_longer(
    cols = starts_with("diag_"),
    names_to = "no_n",
    names_prefix = "diag_",
    values_to = "diag_")
hes_OP_diag_long

hes_OP_diag_long <- hes_OP_diag_long %>%
  rename(HESdiag_org = diag_) %>% # rename the variable diag_to HESdiag_org
  filter(HESdiag_org !="") %>%
  select(-dum)
hes_OP_diag_long

# ------------------------------------------------------------------------------
# 2 Sorting diagnosis codes

hes_OP_diag_long <- hes_OP_diag_long %>%
  mutate(HESdiag_nc4 = HESdiag_org,
         HESdiag_nc3 = substr(HESdiag_org, 1, 3),
         HESdiag_nc1 = substr(HESdiag_org, 1, 1),
         HESdiag_nn2 = substr(HESdiag_org, 2, 3),
         HESdiag_nn3 = substr(HESdiag_org, 2, 4))
hes_OP_diag_long

# some 3 letter codes don't exist (most detail is 2 letter) replace with 2 letter
hes_OP_diag_long <- hes_OP_diag_long %>%
  mutate(templen = nchar(HESdiag_org),
         tempchar = substr(HESdiag_org, 4, 4),
         tempchar2 = substr(HESdiag_org, 5, 5))
hes_OP_diag_long

hes_OP_diag_long$HESdiag_nc4 <- ifelse(hes_OP_diag_long$tempchar2 == "-", substr(hes_OP_diag_long$HESdiag_org, 1, 4), hes_OP_diag_long$HESdiag_nc4)
hes_OP_diag_long$HESdiag_nc4 <- ifelse(hes_OP_diag_long$tempchar2 == "X", substr(hes_OP_diag_long$HESdiag_org, 1, 3), hes_OP_diag_long$HESdiag_nc4)

hes_OP_diag_long$HESdiag_nn3 <- ifelse(hes_OP_diag_long$tempchar == "X", hes_OP_diag_long$HESdiag_nn2, hes_OP_diag_long$HESdiag_nn3)

hes_OP_diag_long$HESdiag_nn2 <- as.numeric(hes_OP_diag_long$HESdiag_nn2)
hes_OP_diag_long$HESdiag_nn3 <- as.numeric(hes_OP_diag_long$HESdiag_nn3)

# create ICD description
icd10HES <- icd10 %>%
  rename(
    HESdiag_nc3 = first,
    HESdiag_d4 = text
  )

hes_OP_diag_long <- left_join(hes_OP_diag_long, icd10HES, by = "HESdiag_nc3")

hes_OP_diag_long$HESdiag_d3 <- hes_OP_diag_long$HESdiag_nc3

hes_OP_diag_long$HESdiag_d4 <- ifelse(hes_OP_diag_long$HESdiag_d4=="" & hes_OP_diag_long$HESdiag_org!="", "no code", hes_OP_diag_long$HESdiag_d4)
hes_OP_diag_long$HESdiag_d3 <- ifelse(hes_OP_diag_long$HESdiag_d3=="" & hes_OP_diag_long$HESdiag_org!="", "no code", hes_OP_diag_long$HESdiag_d3)

hes_OP_diag_long <- hes_OP_diag_long %>%
  select(-templen, -tempchar, -tempchar2)
hes_OP_diag_long

# recode gen data = "OP"

# ------------------------------------------------------------------------------
# 3 Long format operation

hes_OP_oper <- hes_OP

hes_OP_oper <- hes_OP_oper %>%
  mutate_at(vars(starts_with("opertn_")), ~na_if(., "-")) %>%
  mutate_at(vars(starts_with("opertn_")), ~na_if(., "&"))

# check diag is character not numeric and only use columns with data in
summary_opertnvars <- hes_OP_oper %>%
  select(starts_with("opertn"))
summary(summary_opertnvars)

# list of variables to process
varlist <- c('opertn_01', 'opertn_02', 'opertn_03', 'opertn_04', 'opertn_05',
             'opertn_06', 'opertn_07', 'opertn_08', 'opertn_09', 'opertn_10',
             'opertn_11', 'opertn_12', 'opertn_13', 'opertn_14', 'opertn_15',
             'opertn_16', 'opertn_17', 'opertn_18', 'opertn_19', 'opertn_20',
             'opertn_21', 'opertn_22', 'opertn_23', 'opertn_24')

#for (v in varlist) {
 # stopifnot(all(is.na(hes_OP_oper[[v]]))) # will give error if false
#}


# only keep diagnosis data with some data included
hes_OP_oper <- hes_OP_oper %>% select(bgs_id, HESid_n, opertn_01:opertn_24)

# rename the diagnosis codes so can put it into long format
#initialize the dummy variable

hes_OP_oper$dum <- 0

for (i in seq_along(varlist)) {
  
  hes_OP_oper$dum <- hes_OP_oper$dum + 1 # increment the dummy variable
  
  d <- hes_OP_oper$dum[i] # display the current value of dum
  n <- substr(varlist[i], 1, 7) # extract and modify the variable name
  
  new_name <- paste0(n, d)
  
  hes_OP_oper <- hes_OP_oper %>% rename(!!new_name := !!sym(varlist[i]))
  
}

head(hes_OP_oper)

# Reshape the operations code data in HES OP into a long format
hes_OP_oper_long <- hes_OP_oper %>%
  pivot_longer(
    cols = starts_with("opertn_"),
    names_to = "no_n",
    names_prefix = "opertn_",
    values_to = "opertn_")
hes_OP_oper_long

hes_OP_oper_long <- hes_OP_oper_long %>%
  rename(HESoper_org = opertn_) %>% # rename the variable opertn_ to HESoper_org
  filter(HESoper_org !="") %>%
  select(-dum)
hes_OP_oper_long

# ------------------------------------------------------------------------------
# 4 Sorting out operation codes

hes_OP_oper_long <- hes_OP_oper_long %>%
  mutate(HESoper_nc4 = HESoper_org,
         HESoper_nc3 = substr(HESoper_org, 1, 3),
         HESoper_nc1 = substr(HESoper_org, 1, 1),
         HESoper_nn2 = substr(HESoper_org, 2, 3),
         HESoper_nn3 = substr(HESoper_org, 2, 4))
hes_OP_oper_long

# some 3 letter codes don't exist (most detail is 2 letter) replace with 2 letter
hes_OP_oper_long <- hes_OP_oper_long %>%
  mutate(templen = nchar(HESoper_org))
hes_OP_oper_long

# checking all numeric
hes_OP_oper_long$non_numeric <- grepl("[^0-9.-]", hes_OP_oper_long$HESoper_nn3)
hes_OP_oper_long

freq(hes_OP_oper_long$non_numeric) # 152 instances of non numeric in HESoper_nn3

hes_OP_oper_long$HESoper_nn2[hes_OP_oper_long$non_numeric] <- ""
hes_OP_oper_long$HESoper_nn3[hes_OP_oper_long$non_numeric] <- ""

hes_OP_oper_long$HESoper_nn2 <- as.numeric(hes_OP_oper_long$HESoper_nn2)
hes_OP_oper_long$HESoper_nn3 <- as.numeric(hes_OP_oper_long$HESoper_nn3)

hes_OP_oper_long <- hes_OP_oper_long %>% select(-non_numeric, -templen)
hes_OP_oper_long

# ------------------------------------------------------------------------------
# 5 Linking operation codes

# CTSU then link this data set with OPCS49.dta to get OPCS-4 text

opcs4 <- opcs4 %>%
  rename(
    HESoper_nc3 = opcs4code,
    HESoper_d4 = opcs4text
  )

hes_OP_oper_long <- left_join(hes_OP_oper_long, opcs4, by = "HESoper_nc3")

hes_OP_oper_long <- hes_OP_oper_long %>%
  mutate(HESoper_d4 = ifelse((is.na(HESoper_d4) | HESoper_d4 == "") & 
                               (!is.na(HESoper_org) | HESoper_org !=""), "no code", HESoper_d4))

# ------------------------------------------------------------------------------
# HES APC
# ------------------------------------------------------------------------------

hes_APC$temp <- "_APC"
hes_APC <- transform(hes_APC, HESid_n = paste0(APCid_n, temp))
hes_APC <- hes_APC %>% select(-temp)
hes_APC

# 6 Diagnosis

hes_diag_APC <- hes_APC

# check diag is character not numeric and only use columns with data in
summary_diagvars <- hes_diag_APC %>%
  select(starts_with("DIAG"))
summary(summary_diagvars)

hes_diag_APC <- hes_diag_APC %>% select(bgs_id, HESid_n, diag_01:diag_20)

varlist <- c('diag_01', 'diag_02', 'diag_03', 'diag_04', 'diag_05', 'diag_06',
             'diag_07', 'diag_08', 'diag_09', 'diag_10', 'diag_11', 'diag_12',
             'diag_13', 'diag_14', 'diag_15', 'diag_16', 'diag_17', 'diag_18',
             'diag_19', 'diag_20')

# rename the diagnosis codes so can put it into long format
#initialize the dummy variable

hes_diag_APC$dum <- 0

for (i in seq_along(varlist)) {
  
  hes_diag_APC$dum <- hes_diag_APC$dum + 1 # increment the dummy variable
  
  d <- hes_diag_APC$dum[i] # display the current value of dum
  n <- substr(varlist[i], 1, 5) # extract and modify the variable name
  
  new_name <- paste0(n, d)
  
  hes_diag_APC <- hes_diag_APC %>% rename(!!new_name := !!sym(varlist[i]))
  
}

head(hes_diag_APC)

# Reshape the operations code data in HES OP into a long format
hes_diag_APC_long <- hes_diag_APC %>%
  pivot_longer(
    cols = starts_with("diag_"),
    names_to = "no_n",
    names_prefix = "diag_",
    values_to = "diag_")
hes_diag_APC_long

hes_diag_APC_long <- hes_diag_APC_long %>%
  rename(HESdiag_org = diag_) %>% # rename the variable DIAG_ to HESdiag_org
  filter(HESdiag_org !="") %>%
  select(-dum)
hes_diag_APC_long

# ------------------------------------------------------------------------------
# 7 Sorting diagnosis codes

hes_diag_APC_long <- hes_diag_APC_long %>%
  mutate(HESdiag_nc4 = HESdiag_org,
         HESdiag_nc3 = substr(HESdiag_org, 1, 3),
         HESdiag_nc1 = substr(HESdiag_org, 1, 1),
         HESdiag_nn2 = substr(HESdiag_org, 2, 3),
         HESdiag_nn3 = substr(HESdiag_org, 2, 4))
hes_diag_APC_long

# some 3 letter codes don't exist (most detail is 2 letter) replace with 2 letter
hes_diag_APC_long <- hes_diag_APC_long %>%
  mutate(templen = nchar(HESdiag_org),
         tempchar = substr(HESdiag_org, 4, 4))
hes_diag_APC_long

hes_diag_APC_long$HESdiag_nc4 <- ifelse(hes_diag_APC_long$tempchar == "X", substr(hes_diag_APC_long$HESdiag_nc4, 1, 3), hes_diag_APC_long$HESdiag_nc4)
hes_diag_APC_long$HESdiag_nn3 <- ifelse(hes_diag_APC_long$tempchar == "X", hes_diag_APC_long$HESdiag_nn2, hes_diag_APC_long$HESdiag_nn3)

hes_diag_APC_long$HESdiag_nn2 <- as.numeric(hes_diag_APC_long$HESdiag_nn2)
hes_diag_APC_long$HESdiag_nn3 <- as.numeric(hes_diag_APC_long$HESdiag_nn3)

# create ICD description
hes_diag_APC_long <- left_join(hes_diag_APC_long, icd10HES, by = "HESdiag_nc3")

hes_diag_APC_long$HESdiag_d3 <- hes_diag_APC_long$HESdiag_nc3

hes_diag_APC_long$HESdiag_d4 <- ifelse(hes_diag_APC_long$HESdiag_d4=="" & hes_diag_APC_long$HESdiag_org!="", "no code", hes_diag_APC_long$HESdiag_d4)
hes_diag_APC_long$HESdiag_d3 <- ifelse(hes_diag_APC_long$HESdiag_d3=="" & hes_diag_APC_long$HESdiag_org!="", "no code", hes_diag_APC_long$HESdiag_d3)

hes_diag_APC_long <- hes_diag_APC_long %>%
  select(-templen, -tempchar)
hes_diag_APC_long

# ------------------------------------------------------------------------------
# 8 Operations APC

hes_oper_APC <- hes_APC

hes_oper_APC <- hes_oper_APC %>%
  mutate_at(vars(starts_with("opertn_")), ~na_if(., "-")) %>%
  mutate_at(vars(starts_with("opertn_")), ~na_if(., "&"))

# check diag is character not numeric and only use columns with data in
summary_opertnvars <- hes_oper_APC %>%
  select(starts_with("opertn"))
summary(summary_opertnvars)

# list of variables to process
varlist <- c('opertn_01', 'opertn_02', 'opertn_03', 'opertn_04', 'opertn_05',
             'opertn_06', 'opertn_07', 'opertn_08', 'opertn_09', 'opertn_10',
             'opertn_11', 'opertn_12', 'opertn_13', 'opertn_14', 'opertn_15',
             'opertn_16', 'opertn_17', 'opertn_18', 'opertn_19', 'opertn_20',
             'opertn_21', 'opertn_22', 'opertn_23', 'opertn_24')

#for (v in varlist) {
 # stopifnot(all(is.na(hes_oper_APC[[v]]))) # will give error if false
#}


# only keep diagnosis data with some data included
hes_oper_APC <- hes_oper_APC %>% select(bgs_id, HESid_n, starts_with("opertn"))

# rename the diagnosis codes so can put it into long format
#initialize the dummy variable

hes_oper_APC$dum <- 0

for (i in seq_along(varlist)) {
  
  hes_oper_APC$dum <- hes_oper_APC$dum + 1 # increment the dummy variable
  
  d <- hes_oper_APC$dum[i] # display the current value of dum
  n <- substr(varlist[i], 1, 7) # extract and modify the variable name
  
  new_name <- paste0(n, d)
  
  hes_oper_APC <- hes_oper_APC %>% rename(!!new_name := !!sym(varlist[i]))
  
}

head(hes_oper_APC)

# Reshape the operations code data in HES OP into a long format
hes_oper_APC_long <- hes_oper_APC %>%
  pivot_longer(
    cols = starts_with("opertn_"),
    names_to = "no_n",
    names_prefix = "opertn_",
    values_to = "opertn_")
hes_oper_APC_long

hes_oper_APC_long <- hes_oper_APC_long %>%
  rename(HESoper_org = opertn_) %>% # rename the variable opertn_ to HESoper_org
  filter(HESoper_org !="") %>%
  select(-dum)
hes_oper_APC_long

# ------------------------------------------------------------------------------
# 9 Sorting operation codes

hes_oper_APC_long <- hes_oper_APC_long %>%
  mutate(HESoper_nc4 = HESoper_org,
         HESoper_nc3 = substr(HESoper_org, 1, 3),
         HESoper_nc1 = substr(HESoper_org, 1, 1),
         HESoper_nn2 = substr(HESoper_org, 2, 3),
         HESoper_nn3 = substr(HESoper_org, 2, 4))
hes_oper_APC_long

# some 3 letter codes don't exist (most detail is 2 letter) replace with 2 letter
hes_oper_APC_long <- hes_oper_APC_long %>%
  mutate(templen = nchar(HESoper_org))
hes_oper_APC_long

# checking all numeric
hes_oper_APC_long$non_numeric <- grepl("[^0-9.-]", hes_oper_APC_long$HESoper_nn3)
hes_oper_APC_long

freq(hes_oper_APC_long$non_numeric) 

hes_oper_APC_long$HESoper_nn2[hes_oper_APC_long$non_numeric] <- ""
hes_oper_APC_long$HESoper_nn3[hes_oper_APC_long$non_numeric] <- ""

hes_oper_APC_long$HESoper_nn2 <- as.numeric(hes_oper_APC_long$HESoper_nn2)
hes_oper_APC_long$HESoper_nn3 <- as.numeric(hes_oper_APC_long$HESoper_nn3)

hes_oper_APC_long <- hes_oper_APC_long %>% select(-non_numeric, -templen)
hes_oper_APC_long

# ------------------------------------------------------------------------------
# 10 Linking operation codes

# CTSU then link this data set with OPCS49.dta to get OPCS-4 text

hes_oper_APC_long <- left_join(hes_oper_APC_long, opcs4, by = "HESoper_nc3")

hes_oper_APC_long <- hes_oper_APC_long %>%
  mutate(HESoper_d4 = ifelse((is.na(HESoper_d4) | HESoper_d4 == "") & 
                               (!is.na(HESoper_org) | HESoper_org !=""), "no code", HESoper_d4))

# ------------------------------------------------------------------------------
# 11 Combine OP and APC data for diagnosis

# OP
head(hes_OP_diag_long)
head(hes_OP_oper_long)

hes_OP_long <- merge(hes_OP_diag_long, hes_OP_oper_long, by = c('bgs_id', 'HESid_n', 'no_n'),
                      all.x = TRUE)
head(hes_OP_long)

head(hes_OP)
hes_OP_2 <- hes_OP %>%
  arrange(bgs_id, appdate_hesop) %>%
  select(bgs_id, appdate_hesop, OPid_n, mainspef_nc, tretspef_nc, HEStrtcontrol_flag,
         HESmainoncol_flag, mainspef, tretspef)

hes_OP_2$HESdat_n <- hes_OP_2$appdate_hesop
hes_OP_2$OPid_n <- as.character(hes_OP_2$OPid_n)

hes_OP_2$temp <- "_OP"
hes_OP_2 <- transform(hes_OP_2, HESid_n = paste0(OPid_n, temp))
hes_OP_2 <- hes_OP_2 %>% select(-temp, -OPid_n)

head(hes_OP_2)
head(hes_OP_long)

hes_OP_merge <- merge(hes_OP_2, hes_OP_long, by = c('bgs_id', 'HESid_n'),
                      all.x = TRUE)
head(hes_OP_merge)
nrow(hes_OP_merge)

# APC

head(hes_diag_APC_long)
head(hes_oper_APC_long)

hes_APC_long <- merge(hes_diag_APC_long, hes_oper_APC_long, by = c('bgs_id', 'HESid_n', 'no_n'),
                     all.x = TRUE)
head(hes_APC_long)

head(hes_APC)

hes_APC_2 <- hes_APC %>%
  arrange(bgs_id, admidate_hesapc) %>%
  select(bgs_id, admidate_hesapc, APCid_n, mainspef_nc, tretspef_nc, HEStrtoncol_flag,
         mainspef, tretspef)

hes_APC_2$HESdat_n <- hes_APC_2$admidate_hesapc
hes_APC_2$APCid_n <- as.character(hes_APC_2$APCid_n)

hes_APC_2$temp <- "_APC"
hes_APC_2 <- transform(hes_APC_2, HESid_n = paste0(APCid_n, temp))
hes_APC_2 <- hes_APC_2 %>% select(-temp, -APCid_n)

head(hes_APC_2)
head(hes_diag_APC_long)

hes_APC_merge <- merge(hes_APC_2, hes_APC_long, by = c('bgs_id', 'HESid_n'),
                      all.x = TRUE)
head(hes_APC_merge)
nrow(hes_APC_merge)

hes_APC_merge[hes_APC_merge == ""] <- "APC" # some data has no diag or oper codes
n_distinct(hes_APC_merge$bgs_id)

HESlongcoded <- bind_rows(hes_APC_merge, hes_OP_merge)
head(HESlongcoded)
nrow(HESlongcoded)

n_distinct(HESlongcoded$bgs_id)

HESlongcoded$HESdiag_d4 <- tolower(HESlongcoded$HESdiag_d4)
HESlongcoded$HESdiag_d3 <- tolower(HESlongcoded$HESdiag_d3)
HESlongcoded$HESoper_d4 <- tolower(HESlongcoded$HESoper_d4)


##########################
#### END OF CODE BLOCK ###
##########################