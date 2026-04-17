# ********************* CTSU Recurrence Algorithm  ************************


# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 12/06/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# ------------------------------------------------------------------------------

# Code converted from the stata script 8. POETIC COSID tidy
# Purpose: Tidy COSID data

# ------------------------------------------------------------------------------

# merge recurrence and eventstudy

COSD_pt <- COSD_pt %>%
  filter((cosd_rec_date_f > DIAGNOSISDATEBEST) & !is.na(cosd_rec_date) |
           (legacy_date_of_recurrence > DIAGNOSISDATEBEST) & 
           !is.na(legacy_date_of_recurrence))

COSD_pt <- COSD_pt %>%
  distinct(bgs_id, .keep_all = TRUE)

COSD_pt <- COSD_pt %>%
  filter(!(is.na(cosd_rec_date_f) & is.na(legacy_date_of_recurrence) & 
             legacy_recurrence == "" & is.na(legacy_recurrencesource) &
             is.na(legacy_keyworkerseenindicator) &
             is.na(legacy_palliativespecseenind)))

# not clear what following fields are so removing them
COSD_pt <- COSD_pt %>%
  select(-legacy_recurrence, -legacy_recurrencesource, -legacy_keyworkerseenindicator,
         -legacy_palliativespecseenind)

COSD_pt <- COSD_pt %>%
  mutate(COSDrec_dat = cosd_rec_date,
         COSDrec_dat = ifelse(is.na(COSDrec_dat), 
                              as.character(legacy_date_of_recurrence), as.character(COSDrec_dat)))

COSD_pt$COSDrec_dat <- as.Date(COSD_pt$COSDrec_dat)

COSD_pt <- COSD_pt %>%
  mutate(COSDrec_source = ifelse(!is.na(cosd_rec_date), "COSD", NA),
         COSDrec_source = ifelse(is.na(cosd_rec_date) & is.na(legacy_date_of_recurrence),
                                 "COSDlegacy", COSDrec_source),
         COSDrec_source = ifelse(!is.na(cosd_rec_date) & !is.na(legacy_date_of_recurrence),
                                 "Both", COSDrec_source)) %>%
  distinct(bgs_id, .keep_all = TRUE) %>%
  select(bgs_id, tumour_pseudo_id, COSDrec_dat, COSDrec_source) %>%
  arrange(bgs_id, COSDrec_dat) %>%
  group_by(bgs_id) %>%
  mutate(COSD_order = row_number(),
         COSD_no = n()) %>% # total number of treatments
  ungroup()

COSD_pt <- COSD_pt %>%
  pivot_wider(names_from = COSD_order,
              values_from = c(COSDrec_dat, COSDrec_source, COSD_no))

# merge with tumourcantyperec data








