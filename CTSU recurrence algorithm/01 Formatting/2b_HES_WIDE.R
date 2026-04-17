# ********************* CTSU Recurrence Algorithm  ************************


# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 17/06/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# ------------------------------------------------------------------------------
# # Code converted from the stata script 6b POETIC HES wide not coded
# Purpose: Keeping HES ICR and OP codes, dates and speciality data.
# Rename HES ICD and op codes.

# ------------------------------------------------------------------------------

hes_OP$appdate_hesop <- hes_OP$apptdate

hes_OP <- hes_OP %>%
  arrange(bgs_id, appdate_hesop) %>%
  group_by(bgs_id) %>% mutate(OPid_n = row_number()) %>% ungroup()

hes_OP$mainspef_hesop <- hes_OP$mainspef
hes_OP$tretspef_hesop <- hes_OP$tretspef

# Create medical oncology flags and radiology/radiotherapy flags
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

hes_OP <- hes_OP %>%
  select(bgs_id, OPid_n, appdate_hesop, 
         diag_01, diag_02, diag_03, diag_04, diag_05, diag_06, diag_07, 
         diag_08, diag_09, diag_10, diag_11, diag_12, 
         starts_with("opertn_"), mainspef_nc, tretspef_nc, 
         HEStrtoncol_flag, HESmainoncol_flag) 

hes_OP <- hes_OP %>%
  mutate(data = "HES_OP",
         temp = "_OP",
         HESid_n = paste0(OPid_n, temp)) %>%
  select(-temp) %>%
  rename(
    HESdat_n = appdate_hesop
  )

vars <- c('diag_01', 'diag_02', 'diag_03', 'diag_04', 'diag_05', 'diag_06',
       'diag_07', 'diag_08', 'diag_09', 'diag_10', 'diag_11', 'diag_12')

for (i in 1:12) {
  old_name <- vars[i]
  new_name <- paste0("diag4_", i)
  hes_OP <- hes_OP %>% rename(!!new_name := !!sym(old_name))
}

vars2 <- c('opertn_01', 'opertn_02', 'opertn_03', 'opertn_04', 'opertn_05', 'opertn_06',
           'opertn_07', 'opertn_08', 'opertn_09', 'opertn_10', 'opertn_11', 'opertn_12',
           'opertn_13', 'opertn_14', 'opertn_15', 'opertn_16', 'opertn_17', 'opertn_18',
           'opertn_19', 'opertn_20', 'opertn_21', 'opertn_22', 'opertn_23', 'opertn_24')

for (i in 1:24) {
  old_name1 <- vars2[i]
  new_name1 <- paste0("oper4_", i)
  hes_OP <- hes_OP %>% rename(!!new_name1 := !!sym(old_name1))
  
}


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

hes_APC <- hes_APC %>%
  select(bgs_id, APCid_n, admidate_hesapc, starts_with("diag"), 
         starts_with("opertn"), starts_with("opdate"),
         mainspef_nc, tretspef_nc, HEStrtoncol_flag, HESmainoncol_flag) 

hes_APC <- hes_APC %>%
  mutate(data = "HES_APC",
         temp = "_APC",
         HESid_n = paste0(APCid_n, temp)) %>%
  select(-temp) %>%
  rename(
    HESdat_n = admidate_hesapc
  )

vars3 <- c('diag_01', 'diag_02', 'diag_03', 'diag_04', 'diag_05', 'diag_06',
          'diag_07', 'diag_08', 'diag_09', 'diag_10', 'diag_11', 'diag_12',
          'diag_13', 'diag_14', 'diag_15', 'diag_16', 'diag_17', 'diag_18',
          'diag_19')

for (i in 1:19) {
  old_name2 <- vars3[i]
  new_name2 <- paste0("diag4_", i)
  hes_APC <- hes_APC %>% rename(!!new_name2 := !!sym(old_name2))
}

for (i in 1:24) {
  old_name1 <- vars2[i]
  new_name1 <- paste0("oper4_", i)
  hes_APC <- hes_APC %>% rename(!!new_name1 := !!sym(old_name1))
  
}

vars4 <- c('opdate_01', 'opdate_02', 'opdate_03', 'opdate_04', 'opdate_05',
           'opdate_06', 'opdate_07', 'opdate_08', 'opdate_09', 'opdate_10',
           'opdate_11', 'opdate_12', 'opdate_13', 'opdate_14', 'opdate_15',
           'opdate_16', 'opdate_17', 'opdate_18', 'opdate_19', 'opdate_20',
           'opdate_21', 'opdate_22', 'opdate_23', 'opdate_24')

for (i in 1:24) {
  old_name2 <- vars4[i]
  new_name2 <- paste0("operdat_", i)
  hes_APC <- hes_APC %>% rename(!!new_name2 := !!sym(old_name2))
  
}

HESwidenotcoded <- bind_rows(hes_APC, hes_OP)


############################
#### END OF CODE BLOCK #####
###########################