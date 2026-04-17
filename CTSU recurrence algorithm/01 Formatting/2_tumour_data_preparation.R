# ********************* CTSU Recurrence Algorithm  ************************


# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 17/06/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# # Code converted from the stata script 5. Preparing tumour data set
# Purpose: Preparing tumour data set to be used to identify cancer diagnosis

# ------------------------------------------------------------------------------

# 1 Identify breast/non breast cancer diagnosis

# Based on site_icd10_o2_3char and site_coded_desc_fields

n_distinct(tumour_phe$bgs_id)

# how up to date is the diagnosis date?
tumour_phe <- tumour_phe %>%
  mutate(maxdiagdat_n = max(diagnosisdatebest))

class(tumour_phe$maxdiagdat_n)

# coding ICD codes
icd10_tumourphe <- icd10 %>%
  rename(
    site_coded_3char = first,
    icd10text = text
  )

tumour_phe <- left_join(tumour_phe, icd10_tumourphe, by = "site_coded_3char")

# identifying breast in the site_coded_desc field
tumour_phe$site_coded_desc <- tolower(tumour_phe$site_coded_desc)
tumour_phe$icd10text <- tolower(tumour_phe$icd10text)

tumour_phe <- tumour_phe %>%
  mutate(breastdiag_n = grepl("breast", site_coded_desc))

# those with 'breast' in site_coded_desc have ICD code C50 or D05
# and D48 (unspecified code)
freq(tumour_phe$breastdiag_n) # 8110
xtabs(~breastdiag_n + site_coded_3char, data = tumour_phe)

# two type of breast diagnosis under ICD-10 3 char
# malignant neoplasm of breast - c50
# carcinoma in situ of breast - d05
tumour_phe <- tumour_phe %>%
  mutate(breastICD3_flag = ifelse(site_coded_3char == "C50", 1, NA),
         breastICD3_flag = ifelse(site_coded_3char == "D05", 2, breastICD3_flag),
         breastICD3_flag = ifelse(is.na(breastICD3_flag), 0, breastICD3_flag))

freq(tumour_phe$site_coded_3char)
freq(tumour_phe$breastICD3_flag)
freq(tumour_phe$breastdiag_n)

# create any breast cancer flag
tumour_phe <- tumour_phe %>%
  mutate(br_flag = ifelse(
    breastICD3_flag == 1 | breastICD3_flag == 2 | breastdiag_n == "TRUE", 1, NA),
    br_flag = ifelse(is.na(br_flag), 0, br_flag))

freq(tumour_phe$br_flag)

# ------------------------------------------------------------------------------
# Not all cancers included - identifying cancers that need to be removed
# squamous in situ and basal
# benign
# check to see if these have duplicates post the initial breast diagnosis
# ------------------------------------------------------------------------------

table(tumour_phe$behaviour_coded_desc)
tumour_phe$behaviour_coded_desc <- tolower(tumour_phe$behaviour_coded_desc)
table(tumour_phe$histology_coded_desc)
tumour_phe$histology_coded_desc <- tolower(tumour_phe$histology_coded_desc)

tumour_phe$squamous_n <- as.integer(grepl("squamous", tumour_phe$histology_coded_desc))
tumour_phe$basal_n <- as.integer(grepl("basal", tumour_phe$histology_coded_desc))

freq(tumour_phe$squamous_n)
freq(tumour_phe$basal_n)

# basal cell carcinoma based on morph_icd10_o2
freq(tumour_phe$morph_icd10_o2)
class(tumour_phe$morph_icd10_o2)
tumour_phe$morph_icd10_o2 <- as.numeric(tumour_phe$morph_icd10_o2)

tumour_phe <- tumour_phe %>%
  mutate(bassqua_n = ifelse(
    morph_icd10_o2 == 8090 | morph_icd10_o2 == 8091 | morph_icd10_o2 == 8092 |
      morph_icd10_o2 == 8093 | morph_icd10_o2 == 8094 | morph_icd10_o2 == 8095, 1, NA),
  # squamous cell carcinoma based on morph_icd10_o2
  bassqua_n = ifelse(
    morph_icd10_o2 == 8050 | morph_icd10_o2 == 8051 | morph_icd10_o2 == 8052 |
      morph_icd10_o2 == 8070 | morph_icd10_o2 == 8071 | morph_icd10_o2 == 8072 |
      morph_icd10_o2 == 8073 | morph_icd10_o2 == 8074 | morph_icd10_o2 == 8075 |
      morph_icd10_o2 == 8076 | morph_icd10_o2 == 8082, 2, bassqua_n
  ),
  cervix_n = ifelse(site_icd10_o2 == "D069", 1, NA))

freq(tumour_phe$bassqua_n)
freq(tumour_phe$basal_n)
freq(tumour_phe$cervix_n)

# skin cancers based on site_icd10_3_n
tumour_phe$skin_n <- as.integer(grepl("skin", tumour_phe$icd10text))
tumour_phe$melanoma_n <- as.integer(grepl("melanoma", tumour_phe$icd10text))

freq(tumour_phe$skin_n)
freq(tumour_phe$melanoma_n)

# when melanoma and skin just keep the melanoma
tumour_phe <- tumour_phe %>%
  mutate(skin_n = ifelse(melanoma_n == 1, 0, skin_n))

tumour_phe$temp <- as.integer(grepl("malignant neoplasm of", tumour_phe$icd10text))
tumour_phe$temp1 <- as.integer(grepl("carcinoma in situ of", tumour_phe$icd10text))
tumour_phe$temp2 <- as.integer(grepl("other malignant neoplasms of", tumour_phe$icd10text))
tumour_phe$temp3 <- as.integer(grepl("neoplasm of uncertain or unknown behaviour of", tumour_phe$icd10text))
tumour_phe$temp4 <- as.integer(grepl("benign", tumour_phe$icd10text))

freq(tumour_phe$temp)
freq(tumour_phe$temp1)
freq(tumour_phe$temp2)
freq(tumour_phe$temp3)
freq(tumour_phe$temp4)

tumour_phe <- tumour_phe %>%
  mutate(shortdiag_n = ifelse(temp == 1, substr(icd10text, 23, nchar(icd10text)), NA),
         shortdiag_n = ifelse(temp1 == 1, substr(icd10text, 22, nchar(icd10text)), shortdiag_n),
         shortdiag_n = ifelse(temp2 == 1, substr(icd10text, 30, nchar(icd10text)), shortdiag_n),
         shortdiag_n = ifelse(temp3 == 1, substr(icd10text, 47, nchar(icd10text)), shortdiag_n),
         shortdiag_n = ifelse(temp4 == 1, "benign", shortdiag_n)) %>%
  select(-temp, -temp1, -temp2, -temp3, -temp4) %>%
  mutate(shortdiag_n = ifelse(is.na(shortdiag_n), icd10text, shortdiag_n))

freq(tumour_phe$shortdiag_n)

# create field that highlights cancers to be removed
tumour_phe <- tumour_phe %>%
  mutate(benign_n = ifelse(shortdiag_n == "benign", 1, NA),
         benign_n = ifelse(is.na(benign_n), 0, benign_n),
         insitu_n = ifelse(behaviour_coded_desc == "in situ", 1, NA),
         insitu_n = ifelse(is.na(insitu_n), 0, insitu_n))

freq(tumour_phe$benign_n)
freq(tumour_phe$insitu_n)

# Create one field to highlight diagnosis that are invasive or not 
# The following rules were applied:
# - Squamous cell carcinoma has been removed where it was also classified as skin
#   cancer for the site of the disease otherwise, it has been left in.
# - All cancers classified as being 'in situ' have been removed.
# - Cancers classified as being 'in situ' have been removed.
# - Cancers classified as benign have been removed

tumour_phe <- tumour_phe %>%
  mutate(tumdiag_inv = ifelse(
    benign_n == 1 | insitu_n == 1 | basal_n == 1 | (squamous_n == 1 & skin_n == 1), 0, NA
  ),
  tumdiag_inv = ifelse(is.na(tumdiag_inv), 1, tumdiag_inv))

freq(tumour_phe$tumdiag_inv)

# Possible consider removing the squamous cell carcinoma events if sites are 
# unknown there is an odd one in IH as also says site_icd10_4_n = secondary
# malignant neoplasm of bone and bone marrow

# ------------------------------------------------------------------------------
# checking staging data
table(tumour_phe$stage_best)
table(tumour_phe$m_best)
table(tumour_phe$n_best)

# ------------------------------------------------------------------------------

# code NCRAS data using same convention as the study 2nd primary data
# Study data 2nd primaries have been coded

# in future when study data is not available this would be
# ------------------------------------------------------------------------------

# Leukemia
tumour_phe$leuk_n <- as.integer(grepl("leukaemia", tumour_phe$icd10text))
tumour_phe$leuk_n2 <- as.integer(grepl("leukaemia", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(leuk_n == 1, "AML/CLL", NA),
         secprim_n2 = ifelse(leuk_n2 == 1, "AML/CLL", NA))

freq(tumour_phe$leuk_n)
freq(tumour_phe$leuk_n2)
freq(tumour_phe$secprim_n)
freq(tumour_phe$secprim_n2)

# Multiple myeloma
tumour_phe$myel_n <- as.integer(grepl("multiple myeloma", tumour_phe$icd10text))
tumour_phe$myel_n2 <- as.integer(grepl("multiple myeloma", tumour_phe$icd10text))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(myel_n == 1, "multiple myeloma", secprim_n),
         secprim_n2 = ifelse(myel_n2 == 1, "multiple myeloma", secprim_n2))

# Colon
tumour_phe$colo_n <- as.integer(grepl("colon", tumour_phe$icd10text) |
                                  grepl("rectum", tumour_phe$icd10text),
                                grepl("rectosigmoid", tumour_phe$icd10text))

tumour_phe$colo_n2 <- as.integer(grepl("colon", tumour_phe$site_coded_desc) |
                                   grepl("rectum", tumour_phe$site_coded_desc),
                                 grepl("rectosigmoid", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(colo_n == 1, "colon", secprim_n),
         secprim_n2 = ifelse(colo_n2 == 1, "colon", secprim_n2))

# Renal
tumour_phe$rena_n <- as.integer(grepl("renal", tumour_phe$icd10text))
tumour_phe$rena_n2 <- as.integer(grepl("renal", tumour_phe$site_coded_desc) |
                                   grepl("kidney", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(rena_n == 1, "Renal", secprim_n),
         secprim_n2 = ifelse(rena_n2 == 1, "Renal", secprim_n2))

# GI
tumour_phe$gast_n <- as.integer(grepl("gastrointestinal", tumour_phe$icd10text) |
                                  grepl("stomach", tumour_phe$icd10text) |
                                  grepl("intestina", tumour_phe$icd10text) |
                                  grepl("oesophagus", tumour_phe$icd10text))

tumour_phe$gast_n2 <- as.integer(grepl("gastrointestinal", tumour_phe$site_coded_desc) |
                                   grepl("stomach", tumour_phe$site_coded_desc) |
                                   grepl("intestine", tumour_phe$site_coded_desc) |
                                   grepl("oesophagus", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(gast_n == 1, "GI", secprim_n),
         secprim_n2 = ifelse(gast_n2 == 1, "GI", secprim_n2))

# Lung
tumour_phe$lung_n <- as.integer(grepl("lung", tumour_phe$icd10text) |
                                  grepl("mesothelioma", tumour_phe$icd10text))

tumour_phe$lung_n2 <- as.integer(grepl("lung", tumour_phe$site_coded_desc) |
                                  grepl("mesothelioma", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(lung_n == 1, "Lung", secprim_n),
         secprim_n2 = ifelse(lung_n2 == 1, "Lung", secprim_n2))

# Bladder
tumour_phe$blad_n <- as.integer(grepl("bladder", tumour_phe$icd10text) |
                                  grepl("ureter", tumour_phe$icd10text))

tumour_phe$blad_n2 <- as.integer(grepl("bladder", tumour_phe$site_coded_desc))

freq(tumour_phe$blad_n)

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(blad_n == 1, "Bladder", secprim_n),
         secprim_n2 = ifelse(blad_n2 == 1, "Bladder", secprim_n2))

# Lymphoma
tumour_phe$lymp_n <- as.integer(grepl("lymphoma", tumour_phe$icd10text))
tumour_phe$lymp_n2 <- as.integer(grepl("lymphoma", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(lymp_n == 1, "Lymphoma", secprim_n),
         secprim_n2 = ifelse(lymp_n2 == 1, "Lymphoma", secprim_n2))

# Melanoma
tumour_phe$mela_n <- as.integer(grepl("melanoma", tumour_phe$icd10text))
tumour_phe$mela_n2 <- as.integer(grepl("melanoma", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(mela_n == 1, "Melanoma", secprim_n),
         secprim_n2 = ifelse(mela_n2 == 1, "Melanoma", secprim_n2))

# Meningioma
tumour_phe$meni_n <- as.integer(grepl("meningioma", tumour_phe$icd10text) |
                                  grepl("meninges", tumour_phe$icd10text))
tumour_phe$meni_n2 <- as.integer(grepl("meningioma", tumour_phe$histology_coded_desc) |
                                   grepl("meninges", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(meni_n == 1, "Meningioma", secprim_n),
         secprim_n2 = ifelse(meni_n2 == 1, "Meningioma", secprim_n2))

# Ovary
tumour_phe$ovar_n <- as.integer(grepl("ovary", tumour_phe$icd10text))
tumour_phe$ovar_n2 <- as.integer(grepl("ovary", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(ovar_n == 1, "Ovary", secprim_n),
         secprim_n2 = ifelse(ovar_n2 == 1, "Ovary", secprim_n2))

# Peritoneum or retroperitoneum
tumour_phe$retper_n <- as.integer(grepl("peritoneum", tumour_phe$icd10text) |
                                    grepl("retroperitoneum", tumour_phe$icd10text))
tumour_phe$retper_n2 <- as.integer(grepl("peritoneum", tumour_phe$site_coded_desc) |
                                    grepl("retroperitoneum", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(retper_n == 1, "retroperitoneum or peritoneum", secprim_n),
         secprim_n2 = ifelse(retper_n2 == 1, "retroperitoneum or peritoneum", secprim_n2))

# Pancreas
tumour_phe$panc_n <- as.integer(grepl("pancreas", tumour_phe$icd10text))
tumour_phe$panc_n2 <- as.integer(grepl("pancreas", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(panc_n == 1, "Pancreas", secprim_n),
         secprim_n2 = ifelse(panc_n2 == 1, "Pancreas", secprim_n2))

# Sarcoma
tumour_phe$sarc_n <- as.integer(grepl("sarcoma", tumour_phe$icd10text))
tumour_phe$sarc_n2 <- as.integer(grepl("sarcoma", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(sarc_n == 1, "Sarcoma", secprim_n),
         secprim_n2 = ifelse(sarc_n2 == 1, "Sarcoma", secprim_n2))

# Code histology of sarcoma
tumour_phe$sarchist_n <- as.integer(grepl("sarcoma", tumour_phe$histology_coded_desc))

# Thyroid
tumour_phe$thyr_n <- as.integer(grepl("thyroid", tumour_phe$icd10text))
tumour_phe$thyr_n2 <- as.integer(grepl("thyroid", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(thyr_n == 1, "Thyroid", secprim_n),
         secprim_n2 = ifelse(thyr_n2 == 1, "Thyroid", secprim_n2))

# Uterus
tumour_phe$uter_n <- as.integer(grepl("uterus", tumour_phe$icd10text) |
                                  grepl("uteri", tumour_phe$icd10text))
tumour_phe$uter_n2 <- as.integer(grepl("uterus", tumour_phe$site_coded_desc) |
                                   grepl("uteri", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(uter_n == 1, "Uterus", secprim_n),
         secprim_n2 = ifelse(uter_n == 1, "Uterus", secprim_n2))

# Others coding of cancer not in study
# Brain
tumour_phe$brai_n <- as.integer(grepl("brain", tumour_phe$icd10text))
tumour_phe$brai_n2 <- as.integer(grepl("brain", tumour_phe$site_coded_desc))

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(brai_n == 1, "Brain", secprim_n),
         secprim_n2 = ifelse(brai_n2 == 1, "Brain", secprim_n2))

# Female genital
tumour_phe$geni_n <- as.integer(grepl("female genital", tumour_phe$icd10text) |
                                  grepl("vulva", tumour_phe$icd10text) |
                                  grepl("vagina", tumour_phe$icd10text))
tumour_phe$geni_n2 <- as.integer(grepl("female genital", tumour_phe$site_coded_desc) |
                                  grepl("vulva", tumour_phe$site_coded_desc) |
                                  grepl("vagina", tumour_phe$site_coded_desc))

freq(tumour_phe$secprim_n)
freq(tumour_phe$secprim_n2)

tumour_phe <- tumour_phe %>%
  mutate(secprim_n = ifelse(geni_n == 1, "Female genital", secprim_n),
         secprim_n2 = ifelse(geni_n2 == 1, "Female genital", secprim_n2),
         secprim_n = ifelse(skin_n == 1, "Skin", secprim_n),
         secprim_n2 = ifelse(skin_n == 1, "Skin", secprim_n2),
         # others
         secprim_n = ifelse(is.na(secprim_n) & !is.na(diagnosisdatebest) & br_flag !=1, "Other", secprim_n),
         secprim_n2 = ifelse(is.na(secprim_n2) & !is.na(diagnosisdatebest) & br_flag !=1, "Other", secprim_n2),
         secprim_n = ifelse(secprim_n == "Ovary" | secprim_n == "retroperitoneum or peritoneum", 
                            "Ovary/Peritoneal/Fallopian Tube", secprim_n)
  )

# some adjustments on comparison when looking at study POETIC data may need
# additional updates long term to decide on coding sites to be used

tumour_phe$secprim_nfin <- tumour_phe$secprim_n

tumour_phe <- tumour_phe %>%
  mutate(secprim_nfin = ifelse(
    grepl("small intestine", icd10text), "Colon", secprim_nfin),
    secprim_nfin = ifelse(sarchist_n == 1, "Sarcoma", secprim_nfin),
    secprim_nfin = ifelse(secprim_n == "Brain", "Other", secprim_nfin))

# set those that have behviour_coded_desc == "MALIGNANT, METASTATIC/SECONDARY SITE" &
# site_coded_desc == "unknown primary site" to breast cancer and not other.
tumour_phe <- tumour_phe %>%
  mutate(br_flag = ifelse(behaviour_coded_desc == "malignant malignant, metastatic/ secondary site" &
                            site_coded_desc == "unknown primary site" & histology_coded_desc == "adenocarcinoma, metastatic, nos", 1, br_flag))

freq(tumour_phe$secprim_n)
freq(tumour_phe$secprim_n2)
freq(tumour_phe$br_flag)
# ------------------------------------------------------------------------------
# Tidying stage data

table(tumour_phe$stage_best)

# based on data dictionary we will use stage_best for time being
tumour_phe <- tumour_phe %>%
  mutate(tumourmet_n = ifelse((stage_best == "4" | stage_best == "4B"), 1, NA),
         tumourmet_n = ifelse((stage_best == "0" | stage_best == "1" | stage_best == "1A" |
                                stage_best == "1B" | stage_best == "1C" |
                                stage_best == "1E" | stage_best == "2" |
                                stage_best == "2A" | stage_best == "2B" |
                                stage_best == "2C" | stage_best == "3" |
                                stage_best == "3A" | stage_best == "3B" |
                                stage_best == "3C"), 0, tumourmet_n))

freq(tumour_phe$tumourmet_n)

# ------------------------------------------------------------------------------
# STATA code then combines with POETIC event code to 
# Identify whether diagnosis is pre or post randomisation 
# and add date of surgery and randomisation from IMPORT HIGH data
# calculate difference in time from NCRAS diagnosis to surgery study date
# Identify non breast cancer diagnoses post randomisation and classify

tumour_phe <- left_join(tumour_phe, BC_dat, by = "bgs_id", relationship = "many-to-many")

tumour_phe <- tumour_phe %>%
  filter(!is.na(tumourid))

n_distinct(tumour_phe$bgs_id)

# calculate difference in time from diagnosis NCRAS to rand date (BC date)
tumour_phe <- tumour_phe %>%
  mutate(tdiagrand_n = (diagnosisdatebest - bcdat)/7)

# create field with time to diagnosis - study diagnosis to NCRAS diagnosis
tumour_phe <- tumour_phe %>%
  mutate(tdiag_n = (diagnosisdatebest-bcdat)/7)

freq(tumour_phe$tdiag_n)

# ------------------------------------------------------------------------------
# identify non breast cancer diagnosis post randomisation

# time in relation to randomisation
tumour_phe <- tumour_phe %>%
  mutate(trand_n = ifelse(tdiagrand_n <= 0, "Pre", NA),
         trand_n = ifelse(tdiagrand_n > 0 & !is.na(tdiagrand_n), "Post", trand_n))


# ------------------------------------------------------------------------------
# type of diagnosis

# may be different for different studies
tumour_phe <- tumour_phe %>%
  mutate(typediag_tum = ifelse(br_flag == 1 & tdiag_n <= 12, "InitBr", NA),
         typediag_tum = ifelse(br_flag == 1 & tdiag_n > 12, "SecBr", typediag_tum),
         typediag_tum = ifelse(br_flag == 0, "Secprim", typediag_tum))

freq(tumour_phe$typediag_tum)

# ------------------------------------------------------------------------------
# create data set coded data for breast and non breast

# remove all study data fields

tumour_phe <- tumour_phe %>%
  select(-bcdat, -ICDt)

freq(tumour_phe$typediag_tum)
freq(tumour_phe$br_flag)


############################
#### END OF CODE BLOCK #####
############################