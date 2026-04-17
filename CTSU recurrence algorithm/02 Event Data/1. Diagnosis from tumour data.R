# ********************* CTSU Recurrence Algorithm  ************************


# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 25/06/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# ------------------------------------------------------------------------------

# Code converted from the stata script 1. Diagnosis from tumour data

# Purpose:
# 1 Identify all those with a breast cancer diagnosis
# 2 Match to study data
# 3 Look to see how many duplicates with same date and check multifocal or bi-lateral
# 4 Look to see if ER positive
# 5 Look to see if breast cancer diagnosis is within 3 months prior to surgery
# 6 Remove duplicates when different date (unless dates close together)

# ------------------------------------------------------------------------------
# Just those with breast cancer looking at timing 
# ------------------------------------------------------------------------------

# add date of surgery and (randomisation study date) so cancer summary diag date
init <- left_join(tumour_phe, BC_dat, by = "bgs_id", relationship = "many-to-many")

init <- init %>%
  filter(br_flag == 1)

n_distinct(init$bgs_id)

# difference between study rand date and surgery date

# ------------------------------------------------------------------------------
# Initial breast cancer defined as those with diagnosis date within 12 weeks randomisation
# could be 12 weeks of initial surgery
# ------------------------------------------------------------------------------

# keep just those with breast diagnosis in site_icd10_3 and site_coded_desc
init <- tumour_phe %>%
  filter(br_flag == 1) %>%
  filter(tdiag_n <=12)

n_distinct(init$bgs_id)  

# any breast duplicates
init <- init %>%
  mutate(anybr_mult = duplicated(bgs_id) | duplicated(bgs_id, fromLast = TRUE))

init$anybr_mult <- as.integer(init$anybr_mult)
freq(init$anybr_mult)

init <- init %>%
  arrange(bgs_id, diagnosisdatebest)

head(init)

Brinitdiag_pt <- init %>%
  filter(tumdiag_inv == 1)

n_distinct(Brinitdiag_pt$bgs_id)

Brinitdiag_pt <- Brinitdiag_pt %>%
  select(bgs_id, tumourid, diagnosisdatebest, site_icd10_o2, laterality)

Brinitdiag_pt <- Brinitdiag_pt %>%
  mutate(brinitdiag_org = "OK") %>%
  mutate(brinitdiag_mult = duplicated(bgs_id) | duplicated(bgs_id, fromLast = TRUE))

Brinitdiag_pt$brinitdiag_mult <- as.integer(Brinitdiag_pt$brinitdiag_mult)
freq(Brinitdiag_pt$brinitdiag_mult)

Brinitdiag_pt <- Brinitdiag_pt %>%
  group_by(bgs_id) %>%
  mutate(brinitdiag_order = row_number(),
         brinitdiag_no = n()) %>%
  ungroup()

Brinitdiag_pt <- Brinitdiag_pt %>%
  rename(
    brinitdiag_dat = diagnosisdatebest,
    brinitdiag_4icd = site_icd10_o2
  )

head(Brinitdiag_pt)

# ------------------------------------------------------------------------------

Brinitdiag_pt <- Brinitdiag_pt %>%
  pivot_wider(
    id_cols = c(bgs_id),
    names_from = brinitdiag_order,
    values_from = c(
      brinitdiag_dat, brinitdiag_4icd, tumourid, brinitdiag_org, laterality
    )
  )

head(Brinitdiag_pt)

# ------------------------------------------------------------------------------

# update if other non-invasive types

Brinitdiaginv_rec <- init %>%
  filter(tumdiag_inv == 0) %>%
  select(bgs_id, tumourid, diagnosisdatebest, site_icd10_o2) %>%
  rename(
    brinitdiaginv_dat = diagnosisdatebest,
    brinitdiaginv_4icd = site_icd10_o2)

Brinitdiaginv_rec <- Brinitdiaginv_rec %>%
  arrange(bgs_id, brinitdiaginv_dat)

Brinitdiaginv_rec <- Brinitdiaginv_rec %>%
  mutate(brinitdiaginv_mult = duplicated(bgs_id) | duplicated(bgs_id, fromLast = TRUE))

Brinitdiaginv_rec$brinitdiaginv_mult <- as.integer(Brinitdiaginv_rec$brinitdiaginv_mult)
freq(Brinitdiaginv_rec$brinitdiaginv_mult)

Brinitdiaginv_rec <- Brinitdiaginv_rec %>%
  group_by(bgs_id) %>%
  mutate(brinitdiaginv_order = row_number()) %>% ungroup()

head(Brinitdiaginv_rec)

# ------------------------------------------------------------------------------

Brinitdiaginv_pt <- Brinitdiaginv_rec %>%
  pivot_wider(
    id_cols = c(bgs_id),
    names_from = brinitdiaginv_order,
    values_from = c(
      brinitdiaginv_dat, brinitdiaginv_4icd, tumourid, 
    )
  )

head(Brinitdiaginv_pt)

# ------------------------------------------------------------------------------
# 2nd breast diagnosis
# ------------------------------------------------------------------------------

# keep just those with breast diagnosis in site_icd10_3 site_coded_desc
init2nbr <- tumour_phe %>%
  filter(br_flag == 1) %>%
  filter(tdiag_n > 12)

n_distinct(init2nbr$bgs_id)

# any breast duplicates
init2nbr <- init2nbr %>%
  mutate(anybr_mult = duplicated(bgs_id) | duplicated(bgs_id, fromLast = TRUE))

init2nbr$anybr_mult <- as.integer(init2nbr$anybr_mult)
freq(init2nbr$anybr_mult)

init2nbr <- init2nbr %>%
  arrange(bgs_id, diagnosisdatebest)

head(init2nbr)
freq(init2nbr$tumdiag_inv)

# update if other non-invasive types
Br2nddiag_rec <- init2nbr %>%
  filter(tumdiag_inv == 1) %>%
  select(bgs_id, tumourid, diagnosisdatebest, site_icd10_o2, behaviour_coded_desc) %>%
  arrange(bgs_id, diagnosisdatebest)

Br2nddiag_rec <- Br2nddiag_rec %>%
  mutate(br2nddiag_mult = duplicated(bgs_id) | duplicated(bgs_id, fromLast = TRUE))

Br2nddiag_rec$br2nddiag_mult <- as.integer(Br2nddiag_rec$br2nddiag_mult)
freq(Br2nddiag_rec$br2nddiag_mult)

Br2nddiag_rec <- Br2nddiag_rec %>%
  group_by(bgs_id) %>%
  mutate(br2nddiag_order = row_number()) %>%
  ungroup() %>%
  rename(
    br2nddiag_dat = diagnosisdatebest,
    br2nddiag_4icd = site_icd10_o2,
    br2nddiag_behav = behaviour_coded_desc
  )

head(Br2nddiag_rec)

Br2nddiag_pt <- Br2nddiag_rec %>%
  pivot_wider(
    id_cols = c(bgs_id),
    names_from = br2nddiag_order,
    values_from = c(br2nddiag_dat, br2nddiag_4icd, br2nddiag_behav, tumourid)
  )

head(Br2nddiag_pt)

# update if other non-invasive types

Br2nddiagninv_rec <- init2nbr %>%
  filter(tumdiag_inv == 0) %>%
  select(bgs_id, tumourid, diagnosisdatebest, site_icd10_o2, behaviour_coded_desc) %>%
  arrange(bgs_id, diagnosisdatebest)

Br2nddiagninv_rec <- Br2nddiagninv_rec %>%
  mutate(br2nddiagninv_mult = duplicated(bgs_id) | duplicated(bgs_id, fromLast = TRUE))

Br2nddiagninv_rec$br2nddiagninv_mult <- as.integer(Br2nddiagninv_rec$br2nddiagninv_mult)
freq(Br2nddiagninv_rec$br2nddiagninv_mult)

Br2nddiagninv_rec <- Br2nddiagninv_rec %>%
  group_by(bgs_id) %>%
  mutate(br2nddiagninv_order = row_number()) %>%
  ungroup() %>%
  rename(
    br2nddiagninv_dat = diagnosisdatebest,
    br2nddiagninv_4icd = site_icd10_o2,
    br2nddiagninv_behav = behaviour_coded_desc
  )

head(Br2nddiagninv_rec)

Br2nddiagninv_pt <- Br2nddiagninv_rec %>%
  pivot_wider(
    id_cols = c(bgs_id),
    names_from = br2nddiagninv_order,
    values_from = c(br2nddiagninv_dat, br2nddiagninv_4icd, br2nddiagninv_behav, tumourid)
  )

head(Br2nddiagninv_pt)

# ------------------------------------------------------------------------------
# Non breast diagnosis based on tumour data and invasive
# ------------------------------------------------------------------------------

freq(tumour_phe$br_flag)

# keep just those with breast diag in site_icd10_3 site_coded_desc
nonbr <- tumour_phe %>%
  filter(br_flag == 0)

# note also interested if diagnosis is before randomisation as patient could have
# a recurrence of this cancer

n_distinct(nonbr$bgs_id)

# relationship to randomisation
nonbr <- left_join(nonbr, BC_dat, by = "bgs_id", relationship = "many-to-many")

nonbr <- nonbr %>%
  arrange(bgs_id, diagnosisdatebest)

n_distinct(nonbr$bgs_id)

# post randomisation only
NonBrdiag_rec <- nonbr %>%
  filter((diagnosisdatebest >= bcdat) & !is.na(diagnosisdatebest)) %>%
  filter(!is.na(diagnosisdatebest))

n_distinct(NonBrdiag_rec$bgs_id)

# invasive only
freq(NonBrdiag_rec$tumdiag_inv)
NonBrdiag_rec <- NonBrdiag_rec %>%
  filter(tumdiag_inv == 1)

# any breast duplicates
NonBrdiag_rec <- NonBrdiag_rec %>%
  select(bgs_id, tumourid, diagnosisdatebest, site_icd10_o2, secprim_nfin,
         secprim_n, secprim_n2, behaviour_coded_desc) %>%
  rename(
    secprim1_n = secprim_n,
    secprim2_n = secprim_n2
  ) %>%
  arrange(bgs_id, diagnosisdatebest) 

NonBrdiag_rec <- NonBrdiag_rec %>%
  mutate(nonbrdiag_mult = duplicated(bgs_id) | duplicated(bgs_id, fromLast = TRUE))

NonBrdiag_rec$nonbrdiag_mult <- as.integer(NonBrdiag_rec$nonbrdiag_mult)
freq(NonBrdiag_rec$nonbrdiag_mult)

NonBrdiag_rec <- NonBrdiag_rec %>%
  group_by(bgs_id) %>%
  mutate(nonbrdiag_order = row_number()) %>%
  ungroup() %>%
  rename(
    nonbrdiag_dat = diagnosisdatebest,
    nonbrdiag_4icd = site_icd10_o2,
    nonbrdiag_behav = behaviour_coded_desc
  )

n_distinct(NonBrdiag_rec$bgs_id)

NonBrdiag_pt <- NonBrdiag_rec %>%
  pivot_wider(
    id_cols = c(bgs_id),
    names_from = nonbrdiag_order,
    values_from = c(nonbrdiag_dat, nonbrdiag_4icd, nonbrdiag_behav, 
                    secprim_nfin, secprim1_n, secprim2_n,
                    tumourid, nonbrdiag_mult)
  )

head(NonBrdiag_pt)

# ------------------------------------------------------------------------------
# Non BC second primary: update with death details related to 2nd primary where
# tumour data has no second primary
# ------------------------------------------------------------------------------

# death data

NonBrdiag_pt <- NonBrdiag_pt %>%
  rename(
    tumourid = tumourid_1
  )

tumourcatype <- tumour_phe %>%
  select(bgs_id, tumourid, site_icd10_o2, site_coded_desc, histology_coded_desc)

df_list <- list(NonBrdiag_pt, tumourcatype)

NonBrdiag_pt <- df_list %>% reduce(full_join, by = c("bgs_id", "tumourid"))


names(deathwide_full) <- tolower(names(deathwide_full))

deathca <- deathwide_full %>%
  select(bgs_id, deathdatebest, dth2prim_n, dthcan_n, dthcanbr_n, undcan_n,
         undcanbr_n, undcanunk_n, dthmet_n, deathcausetext_1a, deathcausetext_1b,
         deathcausetext_1c, deathcausetext_2, deathcausetext_underlying)

df_list2 <- list(NonBrdiag_pt, deathca)

NonBrdiag_pt <- df_list2 %>% reduce(full_join, by = c("bgs_id"))

# update with death information
NonBrdiag_pt <- NonBrdiag_pt %>%
  mutate(nonbrdiag_flagn = ifelse(!is.na(nonbrdiag_dat_1), 1, NA),
  nonbrdiag_typen = ifelse(nonbrdiag_flagn == 1, "Tumour", NA),
  nonbrdiag_flagn = ifelse(dth2prim_n == "OC", 1, nonbrdiag_flagn))  %>%
  # count as prob breast cancer if unknown cancer in death details
  mutate(nonbrdiag_typen = ifelse(nonbrdiag_flagn == 1 & is.na(nonbrdiag_typen), "Death", nonbrdiag_typen),
         nonbrdiag_dat_1 = as.Date(ifelse(dth2prim_n == "OC" & is.na(nonbrdiag_dat_1), deathdatebest, nonbrdiag_dat_1)))

class(NonBrdiag_pt$nonbrdiag_dat_1)
freq(NonBrdiag_pt$nonbrdiag_flagn)
freq(NonBrdiag_pt$nonbrdiag_typen)

NonBrdiag_pt <- NonBrdiag_pt %>%
  rename(
    tumourid1 = tumourid
  ) %>%
  select(bgs_id, tumourid1, nonbrdiag_dat_1, nonbrdiag_4icd_1, secprim1_n_1, secprim2_n_1,
         secprim_nfin_1, nonbrdiag_mult_1, nonbrdiag_flagn, nonbrdiag_typen,
         nonbrdiag_behav_1)

NonBrdiag_pt <- NonBrdiag_pt %>% filter(!is.na(nonbrdiag_dat_1))

n_distinct(NonBrdiag_pt$bgs_id)

# ------------------------------------------------------------------------------
# Non breast diagnosis based on tumour data pre randomisation and invasive
# ------------------------------------------------------------------------------

# post randomisation only
NonBrdiagprerand_rec <- nonbr %>%
  filter((diagnosisdatebest < bcdat) & !is.na(diagnosisdatebest))

# invasive only
NonBrdiagprerand_rec <- NonBrdiagprerand_rec %>%
  filter(tumdiag_inv == 1)

NonBrdiagprerand_rec <- NonBrdiagprerand_rec %>%
  select(bgs_id, tumourid, diagnosisdatebest, site_icd10_o2, secprim_nfin,
         secprim_n, secprim_n2, behaviour_coded_desc) %>%
  arrange(bgs_id, diagnosisdatebest)

NonBrdiagprerand_rec <- NonBrdiagprerand_rec %>%
  mutate(nonbrdiag_mult = duplicated(bgs_id) | duplicated(bgs_id, fromLast = TRUE))

NonBrdiagprerand_rec$nonbrdiag_mult <- as.integer(NonBrdiagprerand_rec$nonbrdiag_mult)
freq(NonBrdiagprerand_rec$nonbrdiag_mult)

NonBrdiagprerand_rec <- NonBrdiagprerand_rec %>%
  group_by(bgs_id) %>%
  mutate(nonbrdiag_order = row_number()) %>%
  ungroup() %>%
  rename(
    tumourid_prerand = tumourid,
    nonbrdiag_dat_prerand = diagnosisdatebest,
    nonbrdiagn_icd4_prerand = site_icd10_o2,
    secprim_nfin_prerand = secprim_nfin,
    secprim_n_prerand = secprim_n,
    secprim_n2_prerand = secprim_n2,
    nonbrdiag_behav_prerand = behaviour_coded_desc
  )

head(NonBrdiagprerand_rec)

NonBrdiagprerand_pt <- NonBrdiagprerand_rec %>%
  pivot_wider(
    id_cols = c(bgs_id),
    names_from = nonbrdiag_order,
    values_from = c(tumourid_prerand, nonbrdiag_dat_prerand, nonbrdiagn_icd4_prerand,
                    secprim_nfin_prerand, secprim_n_prerand, secprim_n2_prerand,
                    nonbrdiag_behav_prerand)
  )

head(NonBrdiagprerand_pt)

# update if other non-invasive types

NonBrdiagninv_rec <- nonbr %>%
  filter(tumdiag_inv == 0) %>%
  select(bgs_id, tumourid, diagnosisdatebest, site_icd10_o2, secprim_nfin, 
         secprim_n, secprim_n2, behaviour_coded_desc) %>%
  arrange(bgs_id, diagnosisdatebest)

NonBrdiagninv_rec <- NonBrdiagninv_rec %>%
  mutate(nonbrdiag_mult = duplicated(bgs_id) | duplicated(bgs_id, fromLast = TRUE))

NonBrdiagninv_rec$nonbrdiag_mult <- as.integer(NonBrdiagninv_rec$nonbrdiag_mult)
freq(NonBrdiagninv_rec$nonbrdiag_mult)

NonBrdiagninv_rec <- NonBrdiagninv_rec %>%
  group_by(bgs_id) %>%
  mutate(nonbrdiagninv_order = row_number()) %>%
  ungroup() %>%
  rename(
    tumourid_ninv = tumourid,
    nonbrdiag_dat_ninv = diagnosisdatebest,
    nonbrdiag_icd4_ninv = site_icd10_o2,
    secprim_nfin_ninv = secprim_nfin,
    secprim_n_ninv = secprim_n,
    secprim_n2_ninv = secprim_n2,
    nonbrdiag_behav_ninv = behaviour_coded_desc
  )

head(NonBrdiagninv_rec)

NonBrdiagninv_pt <- NonBrdiagninv_rec %>%
  pivot_wider(
    id_cols = c(bgs_id),
    names_from = nonbrdiagninv_order,
    values_from = c(tumourid_ninv, nonbrdiag_dat_ninv, nonbrdiag_icd4_ninv,
                    secprim_nfin_ninv, secprim_n_ninv, secprim_n2_ninv,
                    nonbrdiag_behav_ninv)
  )


############################
#### END OF CODE BLOCK #####
############################