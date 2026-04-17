# ********************* CTSU Recurrence Algorithm  ************************

# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 26/06/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# ------------------------------------------------------------------------------

# Code converted from the stata script 4. events TTR PFS v2

# Purpose: Identify metastatic disease

# ------------------------------------------------------------------------------
# Death from breast cancer indicates disease recurrence
# ------------------------------------------------------------------------------

deathwide_merge <- deathwide_full %>%
  select(bgs_id, deathdatebest, dthcan_n, dthcanbr_n, dthmet_n, undcan_n, 
         undcanbr_n, dth2prim_n, dth_flagn, undcanunk_n, dthcanbr_combn, dthdat_nonbc_n, dthdat_bc_n)

NonBrdiag_pt_merge <- NonBrdiag_pt %>%
  select(bgs_id, nonbrdiag_dat_1, nonbrdiag_typen, secprim_nfin_1)

metshort_merge <- metshort %>%
  select(bgs_id, NCRASmet_dat_fin, NCRASmet_source_fin, NCRASmet_source2add)

df_list_new2 <- list(BC_dat, Brinitdiag_pt, Br2nddiag_pt, NonBrdiag_pt_merge, metshort_merge, lrec, deathwide_merge)

NCRAS_TTR_PFS_full_WOCOSD <- df_list_new2 %>% reduce(full_join, by = 'bgs_id')
n_distinct(metfull$bgs_id) 


#NCRAS_TTR_PFS_full_WOCOSD <- left_join(NCRAS_TTR_PFS_full_WOCOSD, COSD_pt, by = "bgs_id")



# ------------------------------------------------------------------------------

# first recurrence events, recurrence-free survival

class(NCRAS_TTR_PFS_full_WOCOSD$NCRASmet_dat_fin)
class(NCRAS_TTR_PFS_full_WOCOSD$locrec_dat_n)
class(NCRAS_TTR_PFS_full_WOCOSD$dthdat_bc_n)
#class(NCRAS_TTR_PFS_full_WOCOSD$COSDrec_dat1)

NCRAS_TTR_PFS_full_WOCOSD <- NCRAS_TTR_PFS_full_WOCOSD %>%
  mutate(firstrec_dat_n = pmin(NCRASmet_dat_fin, locrec_dat_n, dthdat_bc_n, na.rm = TRUE))

class(NCRAS_TTR_PFS_full_WOCOSD$firstrec_dat_n)


NCRAS_TTR_PFS_full_WOCOSD <- NCRAS_TTR_PFS_full_WOCOSD %>%
  mutate(temptype = "",
         temptype2 = "")


for (v in c("NCRASmet_dat_fin", "locrec_dat_n", "dthdat_bc_n")) {
  
  NCRAS_TTR_PFS_full_WOCOSD <- NCRAS_TTR_PFS_full_WOCOSD %>%
    mutate(
      temptype = ifelse(
        firstrec_dat_n == NCRAS_TTR_PFS_full_WOCOSD[[v]] & !is.na(NCRAS_TTR_PFS_full_WOCOSD[[v]]),
                                               paste0(NCRAS_TTR_PFS_full_WOCOSD[[v]], "/", temptype), temptype),
        temptype2 = ifelse(
          !is.na(NCRAS_TTR_PFS_full_WOCOSD[[v]]),
                 paste0(NCRAS_TTR_PFS_full_WOCOSD[[v]], "/", temptype2), temptype2))

}

# update to make better identifier

NCRAS_TTR_PFS_full_WOCOSD <- NCRAS_TTR_PFS_full_WOCOSD %>%
  mutate(firstrec_type_n = temptype,
         firstrec_alltype_n = temptype2) %>%
  select(-temptype, -temptype2)

NCRAS_TTR_PFS_full_WOCOSD <- NCRAS_TTR_PFS_full_WOCOSD %>%
  mutate(firstrec_type_n = ifelse(is.na(firstrec_dat_n), "No event", firstrec_type_n),
         firstrec_alltype_n = ifelse(is.na(firstrec_dat_n), "No event", firstrec_alltype_n))

# firstrec_dat_n = First breast cancer recurrence date - NCRAS
# firstrec_type_n = First breast cancer recurrence source - NCRAS
# firstrec_type_n = First breast cancer recurrence source all possible events - NCRAS

# First progression free survival DFS-DCIS - disease-free survival ductal carcinoma in situ
class(NCRAS_TTR_PFS_full_WOCOSD$br2nddiag_dat_1)
class(NCRAS_TTR_PFS_full_WOCOSD$nonbrdiag_dat_1)

NCRAS_TTR_PFS_full_WOCOSD <- NCRAS_TTR_PFS_full_WOCOSD %>%
  mutate(firstevent_dat_n = pmin(NCRASmet_dat_fin, locrec_dat_n, dthdat_bc_n,
                               br2nddiag_dat_1, nonbrdiag_dat_1, deathdatebest, na.rm = TRUE)) # and COSDrec_dat1

class(NCRAS_TTR_PFS_full_WOCOSD$firstevent_dat_n)

NCRAS_TTR_PFS_full_WOCOSD <- NCRAS_TTR_PFS_full_WOCOSD %>%
  mutate(temptype = "",
         temptype2 = "")

rm(v)

for (v in c("NCRASmet_dat_fin", "locrec_dat_n", "dthdat_bc_n",
            "br2nddiag_dat_1", "nonbrdiag_dat_1", "deathdatebest")) {
  
  NCRAS_TTR_PFS_full_WOCOSD <- NCRAS_TTR_PFS_full_WOCOSD %>%
    mutate(
      temptype = ifelse(
        firstevent_dat_n == NCRAS_TTR_PFS_full_WOCOSD[[v]] & !is.na(NCRAS_TTR_PFS_full_WOCOSD[[v]]),
        paste0(NCRAS_TTR_PFS_full_WOCOSD[[v]], "/", temptype), temptype),
      temptype2 = ifelse(
        firstevent_dat_n == NCRAS_TTR_PFS_full_WOCOSD[[v]] & !is.na(NCRAS_TTR_PFS_full_WOCOSD[[v]]),
        paste0(NCRAS_TTR_PFS_full_WOCOSD[[v]], "/", temptype2), temptype2))
  
}

# update to make better identifier
NCRAS_TTR_PFS_full_WOCOSD <- NCRAS_TTR_PFS_full_WOCOSD %>%
  mutate(firsteventtype_n = temptype,
         firstevent_alltype_n = temptype2) %>%
  select(-temptype, -temptype2)

NCRAS_TTR_PFS_full_WOCOSD <- NCRAS_TTR_PFS_full_WOCOSD %>%
  mutate(firsteventtype_n = ifelse(is.na(firstevent_dat_n), "No event", firsteventtype_n),
         firstevent_alltype_n = ifelse(is.na(firstevent_dat_n), "No event", firstevent_alltype_n))

# firstevent_dat_n = first disease event date - NCRAS
# firstevent_type_n = first disease event source - NCRAS
# firstevent_allytype_n = first disease event source all possible events - NCRAS

NCRAS_TTR_PFS_full_WOCOSD <- NCRAS_TTR_PFS_full_WOCOSD %>%
  select(bgs_id, firsteventtype_n, firstevent_alltype_n, firstevent_dat_n, 
         firstrec_type_n, firstrec_alltype_n, firstrec_dat_n, locrec_dat_n, 
         NCRASmet_dat_fin, dthdat_bc_n, NCRASmet_source_fin)

freq(NCRAS_TTR_PFS_full_WOCOSD$firsteventtype_n)
freq(NCRAS_TTR_PFS_full_WOCOSD$firstrec_type_n)

n_distinct(NCRAS_TTR_PFS_full_WOCOSD)

# ------------------------------------------------------------------------------

recurrence_data <- NCRAS_TTR_PFS_full_WOCOSD %>%
  filter(firstrec_type_n !="No event")
n_distinct(recurrence_data$bgs_id)

############################
#### END OF CODE BLOCK #####
############################