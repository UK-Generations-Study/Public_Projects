# ********************* CTSU Recurrence Algorithm  ************************


# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 17/06/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# ------------------------------------------------------------------------------

# Code converted from the stata script 9. Identify excision data

# ------------------------------------------------------------------------------

# 1 Identify excision data from TREATMENT data set (assume this means trt summary)

treatment <- treatment %>%
  mutate(opcs4_code = ifelse(opcs4_code == "ABLATION", NA, opcs4_code)) %>%
  mutate(TRTopcs_nc = substr(opcs4_code, 1,1),
         TRTopcs_nn2 = substr(opcs4_code, 2,3),
         TRTopcs_nn3 = substr(opcs4_code, 2,4))

treatment$TRTopcs_nn2 <- as.numeric(treatment$TRTopcs_nn2)
treatment$TRTopcs_nn3 <- as.numeric(treatment$TRTopcs_nn3)

# identification
treatment <- treatment %>%
  mutate(TRTBrex_n1 = 0,
         TRTBrex_n2 = 0,
         TRTBrex_n3 = 0,
         TRTBrex_n4 = 0,
         TRTlymphl_n1 = 0,
         TRTlymphl_n2 = 0,
         TRTlymphl_n3 = 0,
         temp = "")

# creating flag for key procedures relating to excision and dissection for 
# breast and relevant breast lymph nodes

treatment <- treatment %>%
  mutate(temp = str_sub(opcs4_code, 1, -2)) %>%
  mutate(TRTBrex_n1 = ifelse(temp == "B27", 1, NA_integer_),
         TRTBrex_n2 = ifelse(temp == "B28", 10, NA_integer_),
         TRTBrex_n3 = ifelse(temp == "B41", 100, NA_integer_),
         TRTBrex_n4 = ifelse(temp == "T01", 1000, NA_integer_),
         TRTlymphl_n1 = ifelse(opcs4_code == "T852", 2, NA_integer_),
         TRTlymphl_n2 = ifelse(opcs4_code == "T862", 20, NA_integer_),
         TRTlymphl_n3 = ifelse(opcs4_code == "T873", 200, NA_integer_)) 

freq(treatment$TRTBrex_n1)
freq(treatment$TRTBrex_n2)
freq(treatment$TRTBrex_n3)
freq(treatment$TRTBrex_n4)

freq(treatment$TRTlymphl_n1)
freq(treatment$TRTlymphl_n2)
freq(treatment$TRTlymphl_n3)

treatment <- treatment %>%
  mutate(exdat_n = eventdate)

class(treatment$exdat_n)

# add in tumour type so we know if it is initial primary, second primary, second breast primary

tumourcantype <- tumour_phe %>%
  select(bgs_id, tumourid, tumdiag_inv, benign_n, insitu_n, typediag_tum) 

treatment <- left_join(treatment, tumourcantype, by = c("bgs_id", "tumourid"), relationship = "many-to-many")

BC_dat <- BCinv_cohort_conf %>%
  select(bgs_id, DIAGNOSISDATEBEST) %>%
  rename(
    bcdat = DIAGNOSISDATEBEST
  )

treatment <- left_join(treatment, BC_dat, by = "bgs_id", relationship = "many-to-many")


treatment <- treatment %>%
  mutate(trTRTex_52wk = ifelse(eventdate < bcdat+365, 1, NA),
         trTRTex_52wk = ifelse(eventdate >= bcdat+365, 2, trTRTex_52wk))

freq(treatment$trTRTex_52wk)

#  mutate(TRTcount_excbr = ifelse(temp %in% c("B27", "B28", "B41", "T01"), 1, 0),
 #        TRTcount_excl = ifelse(opcs4_code %in% c("T852", "T862", "T873"), 1, 0)) 
#  mutate(TRTcount_excbr = cumsum(TRTcount_excbr),
 #        TRTcount_excl = cumsum(TRTcount_excl))


treatment <- treatment %>%
  filter(TRTBrex_n1 !=0 |
           TRTBrex_n2 !=0 |
           TRTBrex_n3 !=0 |
           TRTBrex_n4 !=0 |
           TRTlymphl_n1 !=0 |
           TRTlymphl_n2 !=0 |
           TRTlymphl_n3 !=0) %>%
  filter(!(is.na(TRTBrex_n1) &
             is.na(TRTBrex_n2) &
             is.na(TRTBrex_n3) &
             is.na(TRTBrex_n4) &
             is.na(TRTlymphl_n1) &
             is.na(TRTlymphl_n2) &
             is.na(TRTlymphl_n3)))

# identifying max means the events are kept across episode as when event not 
# present the value is 0
varlist <- c('TRTBrex_n1', 'TRTBrex_n2', 'TRTBrex_n3', 'TRTBrex_n4', 'TRTlymphl_n1',
             'TRTlymphl_n2', 'TRTlymphl_n3')

rm(v)
replace_max <- function(treatment, varlist) {
  for (v in varlist) {
    treatment <- treatment %>%
      group_by(bgs_id, exdat_n) %>%
      mutate(!!sym(v) := ifelse(all(is.na(!!sym(v))), NA, max(!!sym(v), na.rm = TRUE))) %>%
      ungroup()
  }
  return(treatment)
}

treatment <- replace_max(treatment, varlist)

treatment <- treatment %>%
  distinct(bgs_id, exdat_n, .keep_all = TRUE)

# Creating a field with all the events of interest in one field
treatment <- treatment %>%
  mutate(TRTBrex_count = rowSums(select(., TRTBrex_n1, TRTBrex_n2, TRTBrex_n3, TRTBrex_n4), na.rm = TRUE)) %>%
  mutate(TRTlymph_count = rowSums(select(., TRTlymphl_n1, TRTlymphl_n2, TRTlymphl_n3), na.rm = TRUE))

freq(treatment$TRTlymph_count)
freq(treatment$TRTBrex_count)

trtexcision <- treatment %>%
  select(bgs_id, tumourid, exdat_n, TRTBrex_count, TRTBrex_n1, TRTBrex_n2,
         TRTBrex_n3, TRTBrex_n4, TRTlymph_count, TRTlymphl_n1, TRTlymphl_n2,
         TRTlymphl_n3, tumdiag_inv, typediag_tum)

n_distinct(trtexcision$bgs_id)
# ------------------------------------------------------------------------------
# Identify excision data from HES data set

# STATA code merges HESwidenotcoded and eventsstudy this is POETIC study data


# ------------------------------------------------------------------------------
# Identify breast and lymph node excisions

HESexcision <- HESwidenotcoded %>%
  mutate(HESBrex_n1 = 0,
         HESBrex_n2 = 0,
         HESBrex_n3 = 0,
         HESBrex_n4 = 0,
         HESlymphl_n1 = 0,
         HESlymphl_n2 = 0,
         HESlymphl_n3 = 0,
         temp = "")

process_vars_oper <- function(HESexcision) {
  
  for (v in grep("^oper4_", colnames(HESexcision), value = TRUE)) { 
    
    temp <- substr(HESexcision[[v]], 1, nchar(HESexcision[[v]]) - 1)
    
    # breast excisions
    HESexcision$HESBrex_n1[temp == "B27"] <- 1
    HESexcision$HESBrex_n2[temp == "B28"] <- 10
    HESexcision$HESBrex_n3[temp == "B41"] <- 100
    HESexcision$HESBrex_n4[temp == "T01"] <- 1000
    
    # lymph node block dissection breast
    HESexcision$HESlymphl_n1[HESexcision[[v]] == "T852"] <- 2
    # lymph node sample breast
    HESexcision$HESlymphl_n2[HESexcision[[v]] == "T862"] <- 20
    # lymph node excision breast
    HESexcision$HESlymphl_n3[HESexcision[[v]] == "T873"] <- 200
    
  }
  return(HESexcision)
  }

HESexcision <- process_vars_oper(HESexcision)
print(HESexcision)

HESexcision <- HESexcision %>% select(-temp)


# Distant lymph node excisions -------------------------------------------------
HESexcision <- HESexcision %>%
  mutate(HESlymphd_n1 = 0,
         HESlymphd_n2 = 0,
         HESlymphd_n3 = 0)

process_vars_oper_lymph <- function(HESexcision) {
  
  for (v in grep("^oper4_", colnames(HESexcision), value = TRUE)) { 
    
    temp <- substr(HESexcision[[v]], 1, nchar(HESexcision[[v]]) - 1)
    
    # lymph node block dissection non breast
    HESexcision$HESlymphd_n1[temp == "T85" & !(HESexcision[[v]] %in% c("T852"))] <- 1
    # lymph node sample non breast
    HESexcision$HESlymphd_n2[temp == "T86" & !(HESexcision[[v]] %in% c("T862"))] <- 10
    # lymph node excision non breast 
    HESexcision$HESlymphd_n3[temp == "T87" & !(HESexcision[[v]] %in% c("T873"))] <- 100
    
  }
  return(HESexcision)
  }

HESexcision <- process_vars_oper_lymph(HESexcision)
print(HESexcision)

# Cancer Diagnosis Details -----------------------------------------------------
HESexcision <- HESexcision %>%
  mutate(HESCan_br = 0,
         HESCan_ld = 0,
         HESCan_ll = 0,
         HESCan_lu = 0,
         HESCan_m = 0,
         HESCan_s = 0,
         HESCan_binbr = 0,
         HESCan_bino = 0,
         HESCan_binuk = 0)

# Should update this so d=100 l=10

# Function to perform string operations
process_vars_diag <- function(HESexcision) {
  
  for (v in grep("^diag4_", colnames(HESexcision), value = TRUE)) {
    
    temp <- substr(HESexcision[[v]], 1, nchar(HESexcision[[v]]) - 1)
    temp2 <- substr(HESexcision[[v]], 2, nchar(HESexcision[[v]]) - 1)
    temp3 <- substr(HESexcision[[v]], 1, nchar(HESexcision[[v]]) - 3)
    
    temp2 <- as.numeric(temp2)
    
    HESexcision$HESCan_br[temp == "C50"] <- 1
    HESexcision$HESCan_ld[temp == "C77" & !(HESexcision[[v]] %in% c("C773", "C779"))] <- 100
    HESexcision$HESCan_ll[HESexcision[[v]] == "C773"] <- 10
    HESexcision$HESCan_lu[HESexcision[[v]] == "C779"] <- 1000
    HESexcision$HESCan_m[temp == "C78" | temp == "C79"] <- 10000
    HESexcision$HESCan_s[temp3 == "C" & (temp2 > 00 & temp2 < 77) & temp2!=50] <- 100000
    HESexcision$HESCan_binbr[temp == "D05" | temp == "D24"] <- 1
    HESexcision$HESCan_bino[temp3 == "D" & (temp2 >= 0 & temp2 <= 36) & temp !="D05" & temp !="D24"] <- 1
    HESexcision$HESCan_binuk[temp3 == "D" & (temp2 >= 37 & temp2 <= 48) & temp !="D05" & temp !="D24"] <- 1
    
  }
  
  return(HESexcision)
}

HESexcision <- process_vars_diag(HESexcision)
print(HESexcision)

HESexcision <- HESexcision %>%
  mutate(exdat_n = HESdat_n)
class(HESexcision$exdat_n)

HESexcision <- HESexcision %>%
  filter(HESBrex_n1 !=0 | HESBrex_n2 !=0 | HESBrex_n3 !=0 | HESBrex_n4 !=0 |
           HESlymphl_n1 !=0 | HESlymphl_n2 !=0 | HESlymphl_n3 !=0 |
           HESlymphd_n1 !=0 | HESlymphd_n2 !=0 | HESlymphd_n3 !=0)

# create maximum for each field as duplicates by date will be related to 
# reporting over year probably 

var_names <- c('HESBrex_n1', 'HESBrex_n2', 'HESBrex_n3', 'HESBrex_n4', 'HESlymphl_n1',
                 'HESlymphl_n2', 'HESlymphl_n3', 'HESlymphd_n1', 'HESlymphd_n2', 
                 'HESlymphd_n3', 'HESCan_br', 'HESCan_ld', 'HESCan_ll', 'HESCan_lu',
                 'HESCan_m', 'HESCan_s', 'HESCan_binbr', 'HESCan_bino', 'HESCan_binuk')

replace_max2 <- function(HESexcision, var_names) {
  for (v in var_names) {
    HESexcision <- HESexcision %>%
      group_by(bgs_id, exdat_n) %>%
      mutate(!!sym(v) := ifelse(all(is.na(!!sym(v))), NA, max(!!sym(v), na.rm = TRUE))) %>%
      ungroup()
  }
  return(HESexcision)
}

HESexcision <- replace_max(HESexcision, var_names)

HESexcision <- HESexcision %>%
  distinct(bgs_id, exdat_n, .keep_all = TRUE)

# Row Total Variables

HESexcision <- HESexcision %>%
  # Breast excision (B27, B28, B41, T01) - HES data
  mutate(HESBrex_count = rowSums(select(., HESBrex_n1, HESBrex_n2, HESBrex_n3, HESBrex_n4), na.rm = TRUE)) %>%
  # Lymph node excision local (axillary) (T852, T862, T873) - HES data
  mutate(HESlymphl_count = rowSums(select(., HESlymphl_n1, HESlymphl_n2, HESlymphl_n3), na.rm = TRUE)) %>%
  # Lymph node excision distant (T852, T862, T873) - HES data
  mutate(HESlymphd_count = rowSums(select(., HESlymphd_n1, HESlymphd_n2, HESlymphd_n3), na.rm = TRUE)) %>%
  # HES_br (Breast cancer diag - HES), HESCan_ld (Distant lymph nodes diag - HES)
  # HESCan_ll (Axillary lymph nodes diag - HES), HESCan_lu (Unknown lymph nodes diag - HES),
  # HESCan_m (Metastatic diag - HES) and HESCan_s (Second non-breast cancer diag -HES)
  mutate(HESCan_count = rowSums(select(., HESCan_br, HESCan_ll, HESCan_ld, HESCan_lu, HESCan_m, HESCan_s), na.rm = TRUE)) %>%
  # HESCan_binbr (Benign/in situ breast cancer diag - HES), HESCan_bino (Benign non breast cancer diag - HES)
  # HESCan_binuk (Benign unknown cancer diag - HES)
  mutate(HESCanbi_count = rowSums(select(., HESCan_binbr, HESCan_bino, HESCan_binuk), na.rm = TRUE))

freq(HESexcision$HESBrex_count)
freq(HESexcision$HESlymphl_count)
freq(HESexcision$HESlymphd_count)

HESexcision <- left_join(HESexcision, BC_dat, by = "bgs_id", relationship = "many-to-many")

# Keep those with breast and lymph node excisions
HESexcision <- HESexcision %>%
  filter(HESBrex_count !=0 | HESlymphl_count !=0 | HESlymphd_count !=0) %>%
  mutate(tHESrand = ifelse(exdat_n <= bcdat, 1, NA),
         tHESrand = ifelse((exdat_n > bcdat) & !is.na(exdat_n), 2, tHESrand)) 

freq(HESexcision$tHESrand)

HESexcision <- HESexcision %>%
  filter(!(exdat_n <= bcdat))

HESexcision <- HESexcision %>%
  select(bgs_id, exdat_n, HESBrex_n1, HESBrex_n2, HESBrex_n3, HESBrex_n4,
         HESBrex_count, HESlymphl_n1, HESlymphl_n2, HESlymphl_n3, HESlymphl_count,
         HESlymphd_n1, HESlymphd_n2, HESlymphd_n3, HESlymphd_count, 
         HESCan_br, HESCan_ld, HESCan_ll, HESCan_lu, HESCan_m, HESCan_s,
         HESCan_binbr, HESCan_bino, HESCan_binuk, HESCan_count, HESCanbi_count)

# merge HESexcision and trtexcision to get excisionfull by id and exdat_n

excisionfull <- left_join(HESexcision, trtexcision, 
                          by = c("bgs_id", "exdat_n"), relationship = "many-to-many")

# ------------------------------------------------------------------------------
# Identify breast and lymph node excisions

excisionfull <- left_join(excisionfull, BC_dat, by = "bgs_id", relationship = "many-to-many")

excisionfull <- excisionfull %>%
  filter(!(exdat_n < bcdat)) %>%
  arrange(bgs_id, exdat_n) %>%
  mutate(trevent_n = (exdat_n-bcdat)/7) %>%
  mutate(trevent_52wk = ifelse(trevent_n <= 52, 1, NA),
         trevent_52wk = ifelse((trevent_n > 52) & !is.na(trevent_n), 2, trevent_52wk))

freq(excisionfull$HESBrex_count)
freq(excisionfull$HESlymphd_count)

# over all HES excisions combined in one field
excisionfull <- excisionfull %>%
  mutate(HESEXCrec1_n = 0,
         # HES breast excisions and HES breast cancer diagnosis
         HESEXCrec1_n = ifelse(!is.na(HESBrex_count) & (HESBrex_count > 0 & HESBrex_count < 10000) &
                                 HESCan_br !=0 & !is.na(HESCan_br), 1, HESEXCrec1_n),
         # HES Breast lymph node excisions and HES breast cancer axillary lymph node diagnosis
         HESEXCrec1_n = ifelse(!is.na(HESlymphl_count) & (HESlymphl_count > 0 & HESlymphl_count < 1000) &
                                 HESCan_ll !=0 & !is.na(HESCan_ll), (HESEXCrec1_n + 10), HESEXCrec1_n),
         #HES breast lymph node excisions and HES breast cancer unknown lymph node diagnosis
         HESEXCrec1_n = ifelse(!is.na(HESlymphl_count) & (HESlymphl_count > 0 & HESlymphl_count < 1000) &
                                 HESCan_lu !=0 & !is.na(HESCan_lu), (HESEXCrec1_n + 100), HESEXCrec1_n),
         # HES breast lymph node excisions and HES breast cancer lymph node diagnosis
         HESEXCrec1_n = ifelse(!is.na(HESlymphd_count) & (HESlymphd_count > 0 & HESlymphd_count < 1000) &
                                 HESCan_ld !=0 & !is.na(HESCan_ld), (HESEXCrec1_n + 1000), HESEXCrec1_n),
         # HES distant lymph node excisions and HES unknown lymph node diagnosis
         HESEXCrec1_n = ifelse(!is.na(HESlymphd_count) & (HESlymphd_count > 0 & HESlymphd_count < 1000) &
                                 HESCan_lu !=0 & !is.na(HESCan_lu), (HESEXCrec1_n + 300), HESEXCrec1_n),
         HESEXCrec1_n = ifelse(!is.na(HESlymphd_count) & (HESlymphd_count > 0 & HESlymphd_count < 1000) & 
                                 HESCan_lu !=0 & !is.na(HESCan_lu), (HESEXCrec1_n + 10000), HESEXCrec1_n),
         HESEXCrec1_n = ifelse(HESCan_m !=0 & !is.na(HESCan_m), (HESEXCrec1_n + 100000), HESEXCrec1_n),
         HESEXCrec1_n = ifelse(HESCan_s !=0 & !is.na(HESCan_s), (HESEXCrec1_n + 1000000), HESEXCrec1_n),
         HESEXCrec1_n = ifelse(HESCan_binbr == 1 | HESCan_bino == 1 | HESCan_binuk == 1, (HESEXCrec1_n + 10000000), HESEXCrec1_n))

freq(excisionfull$HESEXCrec1_n)

# 1 - breast cancer, 10 - axillary lymph nodes, 100 - unknown lymph nodes, 
# 1000 - distant lymph nodes, 300 - unknown lymph nodes and distant,
# 1000000 - second primary, 10000000 - benign/in situ

excisionfull <- excisionfull %>%
  mutate(TRTEXCrec1_n = 0,
         TRTEXCrec1_n = ifelse(!is.na(TRTBrex_count) & (TRTBrex_count > 0 & TRTBrex_count < 1000000), (TRTEXCrec1_n + 3), TRTEXCrec1_n),
         TRTEXCrec1_n = ifelse(!is.na(TRTlymph_count) & (TRTlymph_count > 0 & TRTlymph_count < 1000000), (TRTEXCrec1_n + 30), TRTEXCrec1_n)) %>%
  mutate(EXCrec1_n = 0) %>%
  mutate(EXCrec1_n = ifelse(!is.na(HESBrex_count) & (HESBrex_count > 0 & HESBrex_count < 1000000) & 
                              HESCan_br == 1 & !is.na(HESCan_br), 1, EXCrec1_n),
         EXCrec1_n = ifelse((TRTBrex_count > 0 & TRTBrex_count < 1000000), EXCrec1_n + 3, EXCrec1_n),
         EXCrec1_n = ifelse((HESlymphl_count > 0 & HESlymphl_count < 1000000) &
                              HESCan_ll == 10 & !is.na(HESCan_ll), (EXCrec1_n + 10), EXCrec1_n),
         EXCrec1_n = ifelse((TRTlymph_count > 0 & TRTlymph_count < 1000000), (EXCrec1_n + 30), EXCrec1_n),
         EXCrec1_n = ifelse((HESlymphl_count > 0 & HESlymphl_count < 1000000) &
                              HESCan_lu == 1000 & !is.na(HESCan_lu), (EXCrec1_n + 100), EXCrec1_n),
         EXCrec1_n = ifelse((HESlymphd_count > 0 & HESlymphd_count < 1000000) & 
                              HESCan_ld == 100 & !is.na(HESCan_ld), (EXCrec1_n + 1000), EXCrec1_n),
         EXCrec1_n = ifelse((HESlymphd_count > 0 & HESlymphd_count < 1000000) &
                              HESCan_lu == 1000 & !is.na(HESCan_lu), (EXCrec1_n + 300), EXCrec1_n),
         EXCrec1_n = ifelse((HESlymphd_count > 0 & HESlymphd_count < 1000000) & 
                              HESCan_lu == 1000 & !is.na(HESCan_lu), (EXCrec1_n + 10000), EXCrec1_n),
         EXCrec1_n = ifelse(HESCan_m == 10000 & !is.na(HESCan_m), (EXCrec1_n + 100000), EXCrec1_n),
         EXCrec1_n = ifelse(HESCan_s == 100000 & !is.na(HESCan_s), (EXCrec1_n + 1000000), EXCrec1_n),
         EXCrec1_n = ifelse(HESCan_binbr == 1 | HESCan_bino == 1 | HESCan_binuk == 1, (EXCrec1_n + 10000000), EXCrec1_n))

freq(excisionfull$TRTEXCrec1_n)
freq(excisionfull$EXCrec1_n)

# create types of recurrence field
excisionfull <- excisionfull %>%
  mutate(EXCrec1_type = ifelse(HESCan_s !=0 & !is.na(HESCan_s), "Sec", NA),
         EXCrec1_type = ifelse(is.na(EXCrec1_type) & HESCan_m !=0 & !is.na(HESCan_m), "Met", EXCrec1_type),
         EXCrec1_type = 
           ifelse(is.na(EXCrec1_type) & HESCan_ld !=0 & !is.na(HESCan_ld) & !is.na(HESlymphd_count) &
                    (HESlymphd_count > 0 & HESlymphd_count < 1000000), "Metbrl", EXCrec1_type),
         EXCrec1_type = 
           ifelse(is.na(EXCrec1_type) & HESCan_lu !=0 & !is.na(HESCan_lu) & !is.na(HESlymphd_count) &
                    (HESlymphd_count > 0 & HESlymphd_count < 1000000), "Metbrlu", EXCrec1_type),
         EXCrec1_type = 
           case_when(
             is.na(EXCrec1_type) & HESCan_ll !=0 & !is.na(HESCan_ll) & !is.na(HESlymphl_count) &
                    (HESlymphl_count > 0 & HESlymphl_count < 1000000) & trevent_52wk == 2 ~ "locbrl",
             TRUE ~ EXCrec1_type),
         EXCrec1_type = 
           case_when(is.na(EXCrec1_type) & (TRTlymph_count > 0 & TRTlymph_count < 1000000) &
                    trevent_52wk == 2 ~ "locbrl",
                    TRUE ~ EXCrec1_type),
         EXCrec1_type = 
           case_when(is.na(EXCrec1_type) & HESCan_br !=0 & !is.na(HESCan_br) & 
                    (HESBrex_count > 0 & HESBrex_count < 1000000) & trevent_52wk == 2 ~ "locbr",
                    TRUE ~ EXCrec1_type),
         EXCrec1_type = 
           case_when(is.na(EXCrec1_type) & (TRTBrex_count > 0 & TRTBrex_count < 1000000) &
                    trevent_52wk == 2 ~ "locbr",
                    TRUE ~ EXCrec1_type),
         EXCrec1_type = case_when(is.na(EXCrec1_type) & HESCan_ll !=0 & !is.na(HESCan_ll) & 
                                 (HESlymphl_count > 0 & HESlymphl_count < 1000000) &
                                 trevent_52wk == 1 ~ "locbrpre52",
                                TRUE ~ EXCrec1_type),
         EXCrec1_type = case_when(is.na(EXCrec1_type) & 
                                 (TRTlymph_count > 0 & TRTlymph_count < 1000000) &
                                 trevent_52wk == 1 ~ "locbrpre52",
                                 TRUE ~ EXCrec1_type),
         EXCrec1_type = 
           case_when(is.na(EXCrec1_type) & HESCan_br !=0 & !is.na(HESCan_br) & 
                    (HESBrex_count > 0 & HESBrex_count < 1000000) & trevent_52wk==1 ~ "locbrpre52",
                  TRUE ~ EXCrec1_type),
         EXCrec1_type = 
           case_when(is.na(EXCrec1_type) & (TRTBrex_count > 0 & TRTBrex_count < 1000000) &
                    trevent_52wk == 1 ~ "locbrpre52",
                    TRUE ~ EXCrec1_type))

freq(excisionfull$EXCrec1_type)

excisionfull <- excisionfull %>%
  arrange(bgs_id, exdat_n)

# if tumour data indicates 2nd breast primary then update code
excisionfull <- excisionfull %>%
  mutate(EXCrec1_type = 
           ifelse(typediag_tum == "SecBr", "SecBr", EXCrec1_type))

freq(excisionfull$EXCrec1_type)

# calculate if any cancer diagnosis
excisionfull <- excisionfull %>%
  mutate(HEScan_any = rowSums(select(., HESCan_br, HESCan_ld,HESCan_ll, HESCan_lu, HESCan_m, HESCan_s), na.rm = TRUE))

excisionfull <- excisionfull %>%
  filter(!(HEScan_any == 0 & is.na(tumourid)))

excisionfull <- excisionfull %>%
  mutate(trexc_n = (exdat_n - bcdat)/7)

boxdf <- excisionfull %>% filter(HEScan_any !=0 & !is.na(HEScan_any))
ggplot(boxdf, aes(x = factor(trevent_52wk), y = trexc_n)) +
  geom_boxplot()

excisionfull2 <- excisionfull %>%
  filter(trevent_52wk == 2)

# ------------------------------------------------------------------------------
# check those with duplicates, not many within 1st 52 weeks most 2nd events close together

excisionfull2 <- excisionfull2 %>%
  group_by(bgs_id) %>%
  mutate(noext_n = row_number()) %>% ungroup()

excisionfull2 <- excisionfull2 %>%
  select(bgs_id, noext_n, tumourid, exdat_n, bcdat, typediag_tum, tumdiag_inv,
         TRTBrex_count, TRTlymph_count, HESBrex_count, HESlymphl_count, 
         HESlymphd_count, trevent_52wk, HESEXCrec1_n, TRTEXCrec1_n, EXCrec1_n,
         EXCrec1_type, HEScan_any)

excisionfull2 <- excisionfull2 %>%
  pivot_wider(
    id_cols = c(bgs_id, bcdat),
    names_from = noext_n,
    values_from = c(
      tumourid, exdat_n, typediag_tum, tumdiag_inv, TRTBrex_count, TRTlymph_count,
      HESBrex_count, HESlymphl_count, HESlymphd_count, trevent_52wk, HESEXCrec1_n,
      TRTEXCrec1_n, EXCrec1_n, EXCrec1_type, HEScan_any
    )
  )

excisionfull2 <- excisionfull2 %>%
  distinct(bgs_id, .keep_all = TRUE)

n_distinct(excisionfull2$bgs_id)

# number of excision details from HES and Treatment data sets
excisionfull2 <- excisionfull2 %>%
  mutate(HESEXCrec1_yn = 
           ifelse((HESEXCrec1_n_1 > 0 & HESEXCrec1_n_1 < 10000000), 1, NA),
         HESEXCrec1_yn = ifelse(is.na(HESEXCrec1_yn), 0, HESEXCrec1_yn),
         HESEXCrec1_yn = ifelse(HEScan_any_1 == 0, 0, HESEXCrec1_yn))

freq(excisionfull2$HESEXCrec1_yn)

excisionfull2 <- excisionfull2 %>%
  mutate(TRTEXCrec1_yn = 
           ifelse((TRTEXCrec1_n_1 > 0 & TRTEXCrec1_n_1 < 10000000), 1, NA),
         TRTEXCrec1_yn = 
           ifelse(is.na(TRTEXCrec1_yn), 0, TRTEXCrec1_yn))

freq(excisionfull2$TRTEXCrec1_yn)

excisionfull2 <- excisionfull2 %>%
  mutate(EXCrec1_yn = 
           ifelse(EXCrec1_n_1 > 0 & EXCrec1_n_1 < 10000000, 1, NA),
         EXCrec1_yn = 
           ifelse(is.na(EXCrec1_yn), 0, EXCrec1_yn),
         EXCrec1_yn = 
           ifelse(HEScan_any_1 == 0, 0, EXCrec1_yn))

freq(excisionfull2$EXCrec1_yn)

excisionshort <- excisionfull2 %>%
  mutate(EXCrec1_type1new = EXCrec1_type_1,
         EXCrec1_type1new = ifelse(EXCrec1_yn == 0, EXCrec1_type1new == "", EXCrec1_type1new))


freq(excisionshort$EXCrec1_type1new)

############################
#### END OF CODE BLOCK #####
############################