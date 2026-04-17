# ********************* CTSU Recurrence Algorithm  ************************

# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 05/06/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# ------------------------------------------------------------------------------

# Code converted from the stata script 6f. HES flags
# Purpose: Create flags to indicate cancer diagnosis and treatment related to Ca

# ------------------------------------------------------------------------------
# Use HESlongcoded data

class(HESlongcoded$HESdiag_nn2)

can <- HESlongcoded %>%
  arrange(bgs_id, HESdat_n, HESid_n, no_n) %>%
  select(bgs_id, HESdat_n, HESid_n, no_n, everything())

freq(can$HESdiag_nc3)
freq(can$HESoper_nc3)

# 1 Cancer diagnosis flags

can <- can %>%
  mutate(HEScan_flag = ifelse(HESdiag_nc1 == "C", 1, NA)) %>%
  # distant lymph nodes
  mutate(HEScan_flag = ifelse((HESdiag_nn2 >= 77 & HESdiag_nn2 < 78) & HESdiag_nc1 == "C", 3, HEScan_flag),
         # local lymph nodes
         HEScan_flag = ifelse(
           HESdiag_nn3 == 773 & HESdiag_nc1 == "C", 2, HEScan_flag),
         # metastatic disease
         HEScan_flag = ifelse((HESdiag_nn2 >= 78 & HESdiag_nn2 < 80) & 
                                HESdiag_nc1 == "C", 4, HEScan_flag),
         # benign/non invasive
         HEScan_flag = ifelse(
           HESdiag_nc1 == "D" & (HESdiag_nn2 >= 0 & HESdiag_nn2 < 10), 11, HEScan_flag),
         HEScan_flag = ifelse(
           HESdiag_nc1 == "D" & (HESdiag_nn2 >= 10 & HESdiag_nn2 <= 37), 12, HEScan_flag),
         # update should D46 be included - some categories within could be removed
         HEScan_flag = ifelse(
           HESdiag_nc1 == "D" & (HESdiag_nn2 >= 37 & HESdiag_nn2 <= 48), 13, HEScan_flag))
can

freq(can$HEScan_flag)

# add a corresponding description variables with following categories
can <- can %>%
  mutate(HEScan_flag_c = case_when(
    HEScan_flag == 1 ~ "Primary cancer",
    HEScan_flag == 2 ~ "Local lymph nodes",
    HEScan_flag == 3 ~ "Metatstatic lymph nodes",
    HEScan_flag == 4 ~ "Secondary cancer",
    HEScan_flag == 11 ~ "In situ neoplasms",
    HEScan_flag == 12 ~ "Benign neoplasms",
    HEScan_flag == 13 ~ "Neoplasms of uncertain or unknown behaviour"
  ))
can

freq(can$HEScan_flag_c)

# only keep data with a cancer diagnosis
can <- can %>%
  filter(!is.na(HEScan_flag))

n_distinct(can$bgs_id)

# breast cancer flag
#HESlongflag <- HESlongflag %>%
#  mutate(HESbrcan_flag = ifelse(HESdiag_nc1 == "C" & HESdiag_nn2 == 50, 1, 0))

# re-arrange/format the data
can <- can %>%
  arrange(bgs_id, HESdat_n, HESid_n, HEScan_flag) %>%
  group_by(bgs_id, HESdat_n, HESid_n) %>%
  mutate(tempno = row_number()) %>% 
  ungroup()

# Drop HESanycan_site
can <- can %>%
  mutate(HESanycan_site = ifelse(tempno == 1, HESdiag_nc3, NA)) %>%
  mutate(
    HESanycan_site = ifelse(
    HESid_n == lead(HESid_n) & !is.na(lead(HEScan_flag)) & tempno == 1,
    paste(HESanycan_site, lead(HESdiag_nc3), sep = "|"), HESanycan_site),
    HESanycan_site = ifelse(
      HESid_n == lead(HESid_n, 2) & !is.na(lead(HEScan_flag, 2)) & tempno == 1,
      paste(HESanycan_site, lead(HESdiag_nc3, 2), sep = "|"), HESanycan_site),
    HESanycan_site = ifelse(
      HESid_n == lead(HESid_n, 3) & !is.na(lead(HEScan_flag, 3)) & tempno == 1,
      paste(HESanycan_site, lead(HESdiag_nc3, 3), sep = "|"), HESanycan_site),
    HESanycan_site = ifelse(
      HESid_n == lead(HESid_n, 4) & !is.na(lead(HEScan_flag, 4)) & tempno == 1,
      paste(HESanycan_site, lead(HESdiag_nc3, 4), sep = "|"), HESanycan_site),
    HESanycan_site = ifelse(
      HESid_n == lead(HESid_n, 5) & !is.na(lead(HEScan_flag, 5)) & tempno == 1,
      paste(HESanycan_site, lead(HESdiag_nc3, 5), sep = "|"), HESanycan_site),
    HESanycan_site = ifelse(
      HESid_n == lead(HESid_n, 6) & !is.na(lead(HEScan_flag, 6)) & tempno == 1,
      paste(HESanycan_site, lead(HESdiag_nc3, 6), sep = "|"), HESanycan_site),
    HESanycan_site = ifelse(
      HESid_n == lead(HESid_n, 7) & !is.na(lead(HEScan_flag, 7)) & tempno == 1,
      paste(HESanycan_site, lead(HESdiag_nc3, 7), sep = "|"), HESanycan_site))
can

freq(can$HESanycan_site)

# Other Cancer type
can <- can %>%
  mutate(HEScantype_flagc = "") %>%
  mutate(
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C07", "Partoid gland", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C16", "Stomach", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C17", "Colorectal", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C18", "Colorectal", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C19", "Colorectal", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C20", "Colorectal", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C22", "Liver", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C23", "Gall bladder", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C25", "Pancreas", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C34", "Lung", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C38", "Heart, mediastinum and pleura", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C43", "Melanoma", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C44", "Skin", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C45", "Mesothelioma", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C48", "Retroperitoneum and peritoneum", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C50", "Breast", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C51", "Vulva", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C52", "Vagina", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C53", "Cervix uteri", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C54", "Corpus uteri", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C55", "Uterus", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C56", "Ovary", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C57", "Unspecified female genital organs", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C64", "Kidney", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C67", "Bladder", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C71", "Brain", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C72", "Spinal cord", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C77", "Secondary lymph nodes", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C78", "Secondary respiratory and digestive organs", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C79", "Secondary malignant neoplasm of other and unspecified organs", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C80", "Without specification of site", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C81" | HESdiag_nc3 == "C82" | HESdiag_nc3 == "C83" |
        HESdiag_nc3 == "C84" | HESdiag_nc3 == "C85" | HESdiag_nc3 == "C86", "Lymphoma", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C88", "Malignant immunoproliferative diseases", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C90", "Multiple myeloma and malignant plasma neoplasms", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C91" | HESdiag_nc3 == "C92" | HESdiag_nc3 == "C93" | 
        HESdiag_nc3 == "C94" | HESdiag_nc3 == "C95", "Leukaemia", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C96", "Lymphoid, haematopoietic and related tissue", HEScantype_flagc),
    HEScantype_flagc = ifelse(
      HESdiag_nc3 == "C97", "Independent (primary) multiple sites", HEScantype_flagc))

freq(can$HEScantype_flagc)

# drop HESanycan_sitec
can <- can %>%
  mutate(HESanycan_sitec = ifelse(tempno == 1, HEScantype_flagc, NA)) %>%
  mutate(
    HESanycan_sitec = ifelse(
      HESid_n == lead(HESid_n) & !is.na(lead(HEScan_flag)) & tempno == 1,
      paste(HESanycan_sitec, lead(HEScantype_flagc), sep = "|"), HESanycan_sitec),
    HESanycan_sitec = ifelse(
      HESid_n == lead(HESid_n, 2) & !is.na(lead(HEScan_flag, 2)) & tempno == 1,
      paste(HESanycan_sitec, lead(HEScantype_flagc, 2), sep = "|"), HESanycan_sitec),
    HESanycan_sitec = ifelse(
      HESid_n == lead(HESid_n, 3) & !is.na(lead(HEScan_flag, 3)) & tempno == 1,
      paste(HESanycan_sitec, lead(HEScantype_flagc, 3), sep = "|"), HESanycan_sitec),
    HESanycan_sitec = ifelse(
      HESid_n == lead(HESid_n, 4) & !is.na(lead(HEScan_flag, 4)) & tempno == 1,
      paste(HESanycan_sitec, lead(HEScantype_flagc, 4), sep = "|"), HESanycan_sitec),
    HESanycan_sitec = ifelse(
      HESid_n == lead(HESid_n, 5) & !is.na(lead(HEScan_flag, 5)) & tempno == 1,
      paste(HESanycan_sitec, lead(HEScantype_flagc, 5), sep = "|"), HESanycan_sitec),
    HESanycan_sitec = ifelse(
      HESid_n == lead(HESid_n, 6) & !is.na(lead(HEScan_flag, 6)) & tempno == 1,
      paste(HESanycan_sitec, lead(HEScantype_flagc, 6), sep = "|"), HESanycan_sitec),
    HESanycan_sitec = ifelse(
      HESid_n == lead(HESid_n, 7) & !is.na(lead(HEScan_flag, 7)) & tempno == 1,
      paste(HESanycan_sitec, lead(HEScantype_flagc, 7), sep = "|"), HESanycan_sitec)
  )

freq(can$HESanycan_sitec)

# creating overall cancer site
can <- can %>%
  mutate(HEScan_flagtemp = ifelse(HEScan_flag == 13, 1000000, NA)) %>%
  mutate(HEScan_flagtemp = ifelse(HEScan_flag == 12, 100000, HEScan_flagtemp),
         HEScan_flagtemp = ifelse(HEScan_flag == 11, 10000, HEScan_flagtemp),
         HEScan_flagtemp = ifelse(HEScan_flag == 4, 1000, HEScan_flagtemp),
         HEScan_flagtemp = ifelse(HEScan_flag == 3, 100, HEScan_flagtemp),
         HEScan_flagtemp = ifelse(HEScan_flag == 2, 10, HEScan_flagtemp),
         HEScan_flagtemp = ifelse(HEScan_flag == 1, 1, HEScan_flagtemp))

# create unique Id for combinations
can <- can %>%
  group_by(bgs_id, HESid_n) %>%
  mutate(HESanycan_flag = sum(HEScan_flagtemp)) %>%
  ungroup()

can <- can %>%
  distinct(bgs_id, HESid_n, .keep_all = TRUE)   

can <- can %>%
  select(bgs_id, HESid_n, HESanycan_flag, HESanycan_sitec, HESanycan_site) %>%
  arrange(bgs_id, HESid_n, HESanycan_flag, HESanycan_site, HESanycan_sitec)
head(can)

# ------------------------------------------------------------------------------

# 2 Treatment Flags

HESlongflag <- left_join(HESlongcoded, can, by = c("bgs_id", "HESid_n"))

HESlongflag <- HESlongflag %>%
  arrange(bgs_id, HESdat_n)

n_distinct(HESlongflag$bgs_id)

# Identify chemotherapy
#HESlongflag <- HESlongflag %>%
#  mutate(temp = ifelse(grepl("chemo", HESdiag_d4), 1, 0)) %>% 
#  select(-temp)

HESlongflag <- HESlongflag %>%
  # Z511 - chemo for neoplasms
  mutate(HEStrtchemo_flag = ifelse(HESdiag_nc1 == "Z" & HESdiag_nn3 == 511, 1, NA)) %>%
  # Z512 - Other chemo
  mutate(HEStrtchemo_flag = ifelse(HESdiag_nc1 == "Z" & HESdiag_nn3 == 512, 1, HEStrtchemo_flag))

freq(HESlongflag$HEStrtchemo_flag)

# Identify radiotherapy
#HESlongflag <- HESlongflag %>%
#  mutate(temp = ifelse(grepl("radio", HESdiag_d4), 1, 0)) %>% 
#  select(-temp)

#HESlongflag <- HESlongflag %>%
#  mutate(temp = ifelse(grepl("malignant", HESdiag_d4), 1, NA)) %>% 
#  select(-temp)

HESlongflag <- HESlongflag %>%
  # Z510 - RT session
  mutate(HEStrtrt_flag = ifelse(HESdiag_nc1 == "Z" & HESdiag_nn3 == 510, 1, NA))

freq(HESlongflag$HEStrtrt_flag)

# Operation flags

HESlongflag <- HESlongflag %>%
  mutate(HESop_br = ifelse(HESoper_nc1 == "B" & 
                             (HESoper_nn2 >= 27 & HESoper_nn2 <=41), 1, NA),
         HESop_lymp = ifelse(HESoper_nc1 == "T" &
                               (HESoper_nn2 >= 85 & HESoper_nn2 <= 92), 1, NA)) %>%
  mutate(HESsite_flag = ifelse(HESoper_nc3 == "Z15", 1, NA),
         # Lymph nodes
         HESsite_flag = ifelse(HESoper_nc3 == "Z61", 2, HESsite_flag),
         # Liver
         HESsite_flag = ifelse(HESoper_org == "Z301", 3, HESsite_flag),
         # Lung
         HESsite_flag = ifelse(HESoper_org == "Z246" | HESoper_nc3 == "Z52" |
                                 HESoper_org == "Z245", 4, HESsite_flag),
         # Bladder
         HESsite_flag = ifelse(HESoper_org == "Z421", 5, HESsite_flag),
         # Brain
         HESsite_flag = ifelse(HESoper_nc3 == "Z01" | HESoper_nc3 == "Z02" |
                                 HESoper_nc3 == "Z03" | HESoper_nc3 == "Z04" |
                                 HESoper_nc3 == "Z05", 6, HESsite_flag),
         # Bone
         HESsite_flag = ifelse(HESoper_nc1 == "Z" & 
                                 (HESoper_nn2 >= 63 & HESoper_nn2 <= 80), 7, HESsite_flag),
         HESsite_flag = ifelse(HESoper_nc1 == "Z" & is.na(HESsite_flag), 20, HESsite_flag))

freq(HESlongflag$HESop_br)
freq(HESlongflag$HESop_lymp)
freq(HESlongflag$HESsite_flag)

# add a corresponding description variables with following categories
HESlongflag <- HESlongflag %>%
  mutate(HESsite_flag_c = case_when(
    HESsite_flag == 1 ~ "Breast",
    HESsite_flag == 2 ~ "Lymph node",
    HESsite_flag == 3 ~ "Liver",
    HESsite_flag == 4 ~ "Lung",
    HESsite_flag == 5 ~ "Bladder",
    HESsite_flag == 6 ~ "Brain",
    HESsite_flag == 7 ~ "Bone",
    HESsite_flag == 20 ~ "Other"
  ))

freq(HESlongflag$HESsite_flag)

# Treatments
# Y98 Radiology procedures, X65 Radiotherapy delivery, X69 Other radiotherapy
HESlongflag <- HESlongflag %>%
  mutate(HESoprt_flag = ifelse(
    (HESoper_nc1 == "Y" & (HESoper_nn2 == 91 | HESoper_nn2 == 98)) |
      (HESoper_nc1 == "X" & (HESoper_nn2 == 65 | HESoper_nn2 == 69)), 1, NA),
    # X70 Procurement of drugs for chemotherapy neoplasm in Bands 1-5
    # X71 Procurement of drugs for chemotherapy neoplasm in Bands 6-10
    # X72 Delivery of chemotherapy for neoplasm
    # X73 Delivery of oral chemotherapy for neoplasm
    # X74 Other chemotherapy drugs
    HESopchemo_flag = ifelse(
      HESoper_nc1 == "X" & (HESoper_nn2 == 70 | HESoper_nn2 == 71 | 
                              HESoper_nn2 == 72 | HESoper_nn2 == 73 |
                              HESoper_nn2 == 74), 1, NA))

freq(HESlongflag$HESoprt_flag)
freq(HESlongflag$HESopchemo_flag)

freq(HESlongflag$HESanycan_flag)

head(HESlongflag)
nrow(HESlongflag)

##########################
#### END OF CODE BLOCK ###
##########################