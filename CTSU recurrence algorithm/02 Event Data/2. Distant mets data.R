# ********************* CTSU Recurrence Algorithm  ************************


# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 25/06/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# ------------------------------------------------------------------------------

# Code converted from the stata script 2. distant mets data

# Purpose: Identify metastatic disease

# ------------------------------------------------------------------------------
# Death from breast cancer indicates disease recurrence
# ------------------------------------------------------------------------------

dth <- deathwide_full %>%
  select(bgs_id, deathdatebest, dth_flagn, dthcan_n, dthcanbr_n, dthmet_n,
         undcan_n, undcanbr_n, dth2prim_n, deathcausetext_1a, deathcausetext_1b,
         deathcausetext_1c, deathcausetext_2, deathcausetext_underlying,
         dthdat_nonbc_n, dthdat_bc_n)

freq(dth$dthmet_n)


# ------------------------------------------------------------------------------
# look at RTDS data
# ------------------------------------------------------------------------------

RTDS <- left_join(RTDSreg, BC_dat, by = "bgs_id")

RTDS <- RTDS %>%
  filter(RTDSprstdat_n >= bcdat)

RTDS <- RTDS %>%
  filter(rttreatmentregion == "M") %>%
  arrange(bgs_id, RTDSprstdat_n)

RTDS <- RTDS %>%
  mutate(tempdup = duplicated(bgs_id) | duplicated(bgs_id, fromLast = TRUE))

RTDS$tempdup <- as.integer(RTDS$tempdup)
freq(RTDS$tempdup)

RTDS <- RTDS %>%
  select(bgs_id, RTDSprstdat_n, RTDSprednddat, RTDSsite_flag, RTDSsite_flags, rttreatmentregion) %>%
  group_by(bgs_id) %>%
  mutate(RTDSmet_order = row_number(), # create field for order of multiple events
         RTDSmet_no = n()) %>% # total number of treatments
  ungroup()

# reshape so can see were multiple eventsstudy
RTDS <- RTDS %>%
  pivot_wider(
    id_cols = c(bgs_id),
    names_from = RTDSmet_order,
    values_from = c(RTDSprstdat_n, RTDSprednddat, RTDSsite_flag, RTDSsite_flags, 
                    rttreatmentregion, RTDSmet_no)
  )

# ------------------------------------------------------------------------------
# disease recurrence from HES secondary malignancy (diagnosis program)
# ------------------------------------------------------------------------------

HESsec <- HESwidediag_can %>%
  filter(!is.na(HEScandat_sec_1))

HESsec <- HESsec %>%
  select(bgs_id, HEScandat_sec_1, HESdiagsite_sec_1, HESdiagsite_secc_1,
         HEScandat_sec_2, HESdiagsite_sec_2, HESdiagsite_secc_2,
         HESdrecno_n) %>%
  mutate(HEScan_secflag = 1)

# ------------------------------------------------------------------------------
# disease recurrence from HES lymphnodes classified as distant (diagnosis program)
# ------------------------------------------------------------------------------

HESld <- HESwidediag_can %>%
  filter(!is.na(HEScandat_lympd_1))

HESld <- HESld %>%
  select(bgs_id, HEScandat_lympd_1, HESdiagsite_lympd_1, HESdiagsite_lympdc_1) %>%
  mutate(HEScan_ldflag = 1)

# ------------------------------------------------------------------------------
# COSD data
# ------------------------------------------------------------------------------

# COSD <- COSD_pt

# ------------------------------------------------------------------------------
# 2nd primary
# ------------------------------------------------------------------------------

# renamed 2ndprim in STATA code

Second2ndprim <- NonBrdiag_pt %>%
  select(bgs_id, nonbrdiag_dat_1, nonbrdiag_typen, secprim_nfin_1)

# ------------------------------------------------------------------------------
# SACT
# ------------------------------------------------------------------------------

SACT <- SACTreg_metwide

# ------------------------------------------------------------------------------
# Combine metastases data
# some metastases may be related to a second primary so need to work out way to
# replace subsequent metastases events if initial is suspected to be second primary
# ------------------------------------------------------------------------------

# create new combined dataset called metfull 
n_distinct(BC_dat$bgs_id)
n_distinct(HESsec$bgs_id)
n_distinct(HESld$bgs_id)
n_distinct(RTDS$bgs_id)
n_distinct(Second2ndprim$bgs_id)
n_distinct(dth$bgs_id)
n_distinct(SACT$bgs_id)


df_list_new <- list(BC_dat, HESsec, HESld, RTDS, Second2ndprim, dth, SACT)

metfull <- df_list_new %>% reduce(full_join, by = 'bgs_id')
n_distinct(metfull$bgs_id) 

#metfull <- left_join(BC_dat, HESsec, by = "bgs_id")
#n_distinct(metfull$bgs_id)

#metfull <- left_join(metfull, HESld, by = "bgs_id")

#metfull <- left_join(metfull, RTDS, by = "bgs_id")

# metfull <- left_join(metfull, COSD, by = "bgs_id")

#metfull <- left_join(metfull, Second2ndprim, by = "bgs_id")

#metfull <- left_join(metfull, dth, by = "bgs_id")

#metfull <- left_join(metfull, SACT, by = "bgs_id")

head(metfull)
summary(metfull)

sum(is.na(metfull$bcdat))

metfull <- metfull %>%
  filter(!is.na(bcdat))
n_distinct(metfull)

metfull <- metfull %>%
  select(bgs_id, bcdat, HEScandat_sec_1, HESdiagsite_sec_1, HEScandat_lympd_1,
         HESdiagsite_lympd_1, deathdatebest, dthmet_n, nonbrdiag_dat_1, nonbrdiag_typen, everything())

# ------------------------------------------------------------------------------
# calculate met date
# ------------------------------------------------------------------------------

# check any dates do not occur before randomisation/diagnosis date
metfull <- metfull %>%
  mutate(check1 = ifelse(HEScandat_sec_1 < bcdat, 1, 0),
         check2 = ifelse(HEScandat_lympd_1 < bcdat, 1, 0))

freq(metfull$check1)
freq(metfull$check2)

metfull <- metfull %>%
  select(-check1, - check2)

# check which is first date indicating possible mets
metfull <- metfull %>%
  mutate(NCRASmet_dat = HEScandat_sec_1,
         NCRASmet_source = ifelse(!is.na(NCRASmet_dat), "Sec", NA),
         NCRASmet_site = HESdiagsite_secc_1)

class(metfull$NCRASmet_dat) # check it is date
freq(metfull$NCRASmet_source)

# add death if no mets
metfull <- metfull %>%
  mutate(NCRASmet_dat = as.Date(
    ifelse(dthmet_n == "BC mets" & is.na(NCRASmet_source), deathdatebest, NCRASmet_dat))) %>%
  mutate(NCRASmet_source = 
           ifelse(dthmet_n == "BC mets" & is.na(NCRASmet_source), "Dth", NCRASmet_source))

freq(metfull$NCRASmet_dat)
freq(metfull$NCRASmet_source)

# create variable for time between 2nd primary diagnosis and metastases diagnosis
metfull <- metfull %>%
  mutate(tmet2prim_diff1 = (nonbrdiag_dat_1 - NCRASmet_dat)/7)

metfull <- metfull %>%
  mutate(NCRASmet_time1 = ifelse(!is.na(NCRASmet_dat), 1, NA),
         NCRASmet_time1 = ifelse(tmet2prim_diff1 >= -8 & tmet2prim_diff1 <= 8, 2, NCRASmet_time1),
         NCRASmet_time1 = ifelse(tmet2prim_diff1 < -8, 3, NCRASmet_time1),
         NCRASmet_time1 = ifelse(is.na(NCRASmet_time1), 0, NCRASmet_time1))

freq(metfull$NCRASmet_time1)

# add 2nd primary

# update add subsequent 2nd primary if they exist

# if within 8 weeks of met diagnosis assume met diagnosis is for 2nd primary
# if second primary is before assume met is due to other cancer
# if patient died with breast cancer metastatic disease then rest to br met disease unless
# within 8 weeks of second primary event

# update so not counted as met if due to 2nd primary
# but if patient died with met breast cancer assume this br mets
# could look at improving this

freq(metfull$NCRASmet_time1)

metfull <- metfull %>%
  mutate(
    NCRASmet_dat = as.Date(case_when(
      NCRASmet_time1 !=1 & dthmet_n !="BC mets" & !is.na(NCRASmet_dat) ~ NA_Date_,
      TRUE ~ NCRASmet_dat))) %>%
  mutate(NCRASmet_source = case_when(
    is.na(NCRASmet_dat) & !is.na(NCRASmet_source) ~ "",
    TRUE ~ NCRASmet_source
  ))

# ------------------------------------------------------------------------------

freq(metfull$NCRASmet_source)
freq(metfull$NCRASmet_dat)

metfull$year_met <- as.numeric(format(metfull$NCRASmet_dat, "%Y"))
freq(metfull$year_met)

metfull <- metfull %>%
  mutate(NCRASmet_dat1 = NCRASmet_dat,
         NCRASmet_source1 = NCRASmet_source,
         NCRASmet_source1add = NCRASmet_source)

class(metfull$NCRASmet_dat1)
freq(metfull$NCRASmet_source1)

metfull <- metfull %>%
  mutate(NCRASmet_source1 = ifelse(!is.na(HEScandat_lympd_1) & !is.na(NCRASmet_dat) &
                                     HEScandat_lympd_1 < NCRASmet_dat, "Lymphd", NCRASmet_source1),
         NCRASmet_source1add = ifelse(!is.na(HEScandat_lympd_1),
                                      paste0(NCRASmet_source, "/Lymphd"), NCRASmet_source1add),
         NCRASmet_dat1 = as.Date(ifelse(HEScandat_lympd_1 < NCRASmet_dat, HEScandat_lympd_1, NCRASmet_dat1)))

freq(metfull$NCRASmet_source1)
freq(metfull$NCRASmet_source1add)
class(metfull$NCRASmet_dat1)

metfull <- metfull %>%
  mutate(NCRASmet_dat3 = NCRASmet_dat,
         NCRASmet_source3 = NCRASmet_source,
         NCRASmet_source3add = NCRASmet_source)

class(metfull$NCRASmet_dat3)
freq(metfull$NCRASmet_source3)

metfull <- metfull %>%
  mutate(NCRASmet_source3 = ifelse(!is.na(HEScandat_lympd_1) & !is.na(NCRASmet_dat) & 
                                     HEScandat_lympd_1 < NCRASmet_dat, "Lymphd", NCRASmet_source3),
         NCRASmet_source3add = ifelse(!is.na(HEScandat_lympd_1),
                                      paste0(NCRASmet_source, "/Lymphd"), NCRASmet_source3add),
         NCRASmet_dat1 = as.Date(ifelse(HEScandat_lympd_1 < NCRASmet_dat, HEScandat_lympd_1, NCRASmet_dat1)))

freq(metfull$NCRASmet_source3)
freq(metfull$NCRASmet_source3add)

# also does same for NCRASmet_source3 and NCRASmet_source3add for COSD

# for other diagnosis
# only if RTDS is before second primary
metfull <- metfull %>%
  mutate(RTDSprstdat_n_1 = ifelse(((nonbrdiag_dat_1-RTDSprstdat_n_1)/7) < 8, NA, RTDSprstdat_n_1),
         start_date_of_regimen_1 = ifelse(((nonbrdiag_dat_1-start_date_of_regimen_1)/7) < 8, NA, start_date_of_regimen_1))

freq(metfull$RTDSprstdat_n_1)

metfull <- metfull %>%
  mutate(NCRASmet_dat2 = NCRASmet_dat,
         NCRASmet_source2 = NCRASmet_source,
         NCRASmet_source2add = NCRASmet_source)

class(metfull$NCRASmet_dat2)
freq(metfull$NCRASmet_source2)

metfull <- metfull %>%
  mutate(NCRASmet_source2add = ifelse(!is.na(HEScandat_lympd_1),
                                      paste0(NCRASmet_source, "/Lymphd"), NCRASmet_source2add),
         NCRASmet_source2 = ifelse(!is.na(HEScandat_lympd_1) & !is.na(NCRASmet_dat) & 
                                     HEScandat_lympd_1 < NCRASmet_dat, "Lymphd", NCRASmet_source2),
         NCRASmet_dat2 = as.Date(ifelse(HEScandat_lympd_1 < NCRASmet_dat, HEScandat_lympd_1, NCRASmet_dat2)),
         NCRASmet_source2add = ifelse(!is.na(RTDSprstdat_n_1),
                                      paste0(NCRASmet_source, "/RTDS"), NCRASmet_source2add),
         NCRASmet_source2 = ifelse(!is.na(RTDSprstdat_n_1) & !is.na(NCRASmet_dat) & 
                                     RTDSprstdat_n_1 < NCRASmet_dat, "RTDS", NCRASmet_source2),
         NCRASmet_dat2 = as.Date(ifelse(RTDSprstdat_n_1 < NCRASmet_dat, RTDSprstdat_n_1, NCRASmet_dat2)))

freq(metfull$NCRASmet_source2)
freq(metfull$NCRASmet_source2add)

metfull <- metfull %>%
  mutate(NCRASmet_dat4 = NCRASmet_dat,
         NCRASmet_source4 = NCRASmet_source,
         NCRASmet_source4add = NCRASmet_source)

class(metfull$NCRASmet_dat4)
freq(metfull$NCRASmet_source4)

metfull <- metfull %>%
  mutate(NCRASmet_source4add = ifelse(!is.na(HEScandat_lympd_1),
                                      paste0(NCRASmet_source, "/Lymphd"), NCRASmet_source4add),
         NCRASmet_source4 = ifelse(!is.na(HEScandat_lympd_1) & !is.na(NCRASmet_dat) &
                                     HEScandat_lympd_1 < NCRASmet_dat, "Lymphd", NCRASmet_source4),
         NCRASmet_dat4 = as.Date(ifelse(HEScandat_lympd_1 < NCRASmet_dat, HEScandat_lympd_1, NCRASmet_dat4)),
         NCRASmet_source4add = ifelse(!is.na(RTDSprstdat_n_1),
                                      paste0(NCRASmet_source, "/RTDS"), NCRASmet_source4add),
         NCRASmet_source4 = ifelse(!is.na(RTDSprstdat_n_1) & !is.na(NCRASmet_dat) &
                                     RTDSprstdat_n_1 < NCRASmet_dat, "RTDS", NCRASmet_source4),
         NCRASmet_dat4 = as.Date(ifelse(RTDSprstdat_n_1 < NCRASmet_dat, RTDSprstdat_n_1, NCRASmet_dat4)),
         NCRASmet_source4add = ifelse(!is.na(start_date_of_regimen_1),
                                      paste0(NCRASmet_source, "/SACT"), NCRASmet_source4add),
         NCRASmet_source4 = ifelse(!is.na(start_date_of_regimen_1) & !is.na(NCRASmet_dat) &
                                     (start_date_of_regimen_1 < NCRASmet_dat) & 
                                     !is.na(start_date_of_regimen_1), "SACT", NCRASmet_source4),
         NCRASmet_dat4 = as.Date(ifelse(start_date_of_regimen_1 < NCRASmet_dat, start_date_of_regimen_1, NCRASmet_dat4)))

freq(metfull$NCRASmet_source4)
freq(metfull$NCRASmet_source4add)

freq(metfull$NCRASmet_dat)
# ------------------------------------------------------------------------------
# create final met date based on presuming data complete up to maximum diagnosis
# date as some HES date and deaths post this date but assume not complete

metfull <- metfull %>%
  mutate(NCRASmet_dat_fin = NCRASmet_dat, # metastatic cancer date
         NCRASmet_source_fin = NCRASmet_source,
         NCRASmet_source_fin = ifelse(is.na(NCRASmet_dat_fin), NA, NCRASmet_source_fin))

class(metfull$NCRASmet_dat_fin)
freq(metfull$NCRASmet_source_fin)

metshort <- metfull %>%
  arrange(bgs_id, NCRASmet_dat_fin, NCRASmet_source_fin, NCRASmet_dat, NCRASmet_source,
          NCRASmet_site, NCRASmet_source2add, NCRASmet_source2, NCRASmet_dat2, NCRASmet_time1,
          tmet2prim_diff1) %>%
  select(bgs_id, NCRASmet_dat_fin, NCRASmet_source_fin, NCRASmet_dat, NCRASmet_source,
         NCRASmet_site, NCRASmet_source2add, NCRASmet_source2, NCRASmet_dat2, NCRASmet_time1,
         tmet2prim_diff1)

metshort <- metshort %>%
  filter(!is.na(NCRASmet_dat) | !is.na(NCRASmet_dat2))


############################
#### END OF CODE BLOCK #####
############################