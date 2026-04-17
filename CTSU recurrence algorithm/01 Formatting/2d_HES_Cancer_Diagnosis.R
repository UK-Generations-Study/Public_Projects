# ********************* CTSU Recurrence Algorithm  ************************


# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 05/06/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# ------------------------------------------------------------------------------

# Code converted from the stata script 6g. HES Cancer diagnosis
# Purpose: Create flags to indicate cancer diagnosis and treatment related to Ca

# ------------------------------------------------------------------------------
# Large section of stata script using eventstudy data - can't apply for our cohort

# 1 Non breast primary cancer diagnosis post randomisation
# patient could have recurrence of 2nd primary so keep separating them out

# Use HESlongflag data
n_distinct(HESlongflag$bgs_id)

nonbrprim <- HESlongflag %>%
  filter(HESanycan_flag == 1 & HESanycan_site !="C50") %>%
  filter(HESdiag_nc1 == "C")

freq(nonbrprim$HESanycan_site)
n_distinct(nonbrprim$bgs_id)
head(nonbrprim)

# add randomisation date as only interested in 2nd primary and recurrence post this date
# used breast cancer diagnosis date instead of randomisation date for cohort study application
eventdta <- BCinv_cohort_conf %>% select(bgs_id, DIAGNOSISDATEBEST, her2_Status)

nonbrprim <- left_join(nonbrprim, eventdta, by = "bgs_id")

nonbrprim <- nonbrprim %>% select(1:5, DIAGNOSISDATEBEST, everything())

# Calculate time from randomisation
nonbrprim <- nonbrprim %>%
  mutate(tranddiag_n = (HESdat_n-DIAGNOSISDATEBEST)/7)

#ggplot(nonbrprim, aes(y = tranddiag_n)) +
#  geom_boxplot() + labs(title = "Boc plot of tranddiag_n", y = "tranddiag_n") +
#  theme_minimal()

# Keep post diag date (randomisation in ctsu code) only
nonbrprim <- nonbrprim %>%
  filter(!(HESdat_n < DIAGNOSISDATEBEST))

n_distinct(nonbrprim$bgs_id)

#ggplot(nonbrprim, aes(y = tranddiag_n)) +
#  geom_boxplot() + labs(title = "Boc plot of tranddiag_n", y = "tranddiag_n") +
#  theme_minimal()

# remove non-invasive cancer types
nonbrprim$skin_n <- as.integer(grepl("skin", nonbrprim$HESdiag_d4))
nonbrprim$melanoma_n <- as.integer(grepl("melanoma", nonbrprim$HESdiag_d4))
nonbrprim$insitu_n <- as.integer(grepl("in situ", nonbrprim$HESdiag_d4))
nonbrprim$benign_n <- as.integer(grepl("benign", nonbrprim$HESdiag_d4))

freq(nonbrprim$skin_n)
freq(nonbrprim$melanoma_n)
freq(nonbrprim$insitu_n)
freq(nonbrprim$benign_n)

# remove non invasive cancers
nonbrprim <- nonbrprim %>%
  filter(!(skin_n == 1 & melanoma_n == 0)) %>%
  filter(!(insitu_n == 1 | benign_n == 1))

# post randomisation
nonbrprim <- nonbrprim %>%
  filter(tranddiag_n >= 0)

n_distinct(nonbrprim$bgs_id)

nonbrprim <- nonbrprim %>%
  arrange(bgs_id, HESdiag_nc4, HESdat_n) %>%
  distinct(bgs_id, HESdiag_nc4, .keep_all = TRUE) %>%
  distinct(bgs_id, HESdiag_nc3, .keep_all = TRUE)

nonbrprim <- nonbrprim %>%
  rename(
    HEScandat_prim = HESdat_n,
    HESdiagsite_prim = HESdiag_nc3,
    HESdiagsite_primc = HESdiag_d3
  )

# Create number different diagnosis sites
nonbrprim <- nonbrprim %>%
  arrange(bgs_id, HEScandat_prim) %>%
  group_by(bgs_id) %>%
  mutate(HESdiagorder_prim = row_number()) %>%
  ungroup()

nonbrprim <- nonbrprim %>%
  select(bgs_id, HEScandat_prim, HESdiagsite_prim, HESdiagorder_prim, HESdiagsite_primc)

nonbrprim <- nonbrprim %>%
  pivot_wider(
    id_cols = bgs_id,
    names_from = HESdiagorder_prim,
    values_from = c(HEScandat_prim, HESdiagsite_prim, HESdiagsite_primc)
  )
head(nonbrprim)

# ------------------------------------------------------------------------------
# 2 Non Breast Primary Cancer diagnosis pre rand - can't do either

nonbrprimpre <- HESlongflag %>%
  filter(HESanycan_flag == 1 & HESanycan_site !="C50") %>%
  filter(HESdiag_nc1 == "C")

nonbrprimpre <- left_join(nonbrprimpre, eventdta, by = "bgs_id")

# patient could have recurrence of 2nd primary so keep pre randomization too
nonbrprimpre <- nonbrprimpre %>%
  filter(!(HESdat_n > DIAGNOSISDATEBEST))

# calculate time from randomisation
nonbrprimpre <- nonbrprimpre %>%
  mutate(tranddiag_n = (HESdat_n-DIAGNOSISDATEBEST)/7)

#ggplot(nonbrprimpre, aes(y = tranddiag_n)) +
#  geom_boxplot() + labs(title = "Boc plot of tranddiag_n", y = "tranddiag_n") +
#  theme_minimal()

# remove non-invasive cancer types
nonbrprimpre$skin_n <- as.integer(grepl("skin", nonbrprimpre$HESdiag_d4))
nonbrprimpre$melanoma_n <- as.integer(grepl("melanoma", nonbrprimpre$HESdiag_d4))
nonbrprimpre$insitu_n <- as.integer(grepl("in situ", nonbrprimpre$HESdiag_d4))
nonbrprimpre$benign_n <- as.integer(grepl("benign", nonbrprimpre$HESdiag_d4))

freq(nonbrprimpre$skin_n)
freq(nonbrprimpre$melanoma_n)
freq(nonbrprimpre$insitu_n)
freq(nonbrprimpre$benign_n)

# remove non invasive cancers
nonbrprimpre <- nonbrprimpre %>%
  filter(!(skin_n == 1 & melanoma_n == 0)) %>%
  filter(!(insitu_n == 1 | benign_n == 1))

nonbrprimpre <- nonbrprimpre %>%
  arrange(bgs_id, HESdiag_nc4, HESdat_n) %>%
  distinct(bgs_id, HESdiag_nc4, .keep_all = TRUE) %>%
  distinct(bgs_id, HESdiag_nc3, .keep_all = TRUE)

nonbrprimpre <- nonbrprimpre %>%
  rename(
    HEScandat_prim = HESdat_n,
    HESdiagsite_prim = HESdiag_nc3,
    HESdiagsite_primc = HESdiag_d3
  )

# Create number different diagnosis sites
nonbrprimpre <- nonbrprimpre %>%
  arrange(bgs_id, HEScandat_prim) %>%
  group_by(bgs_id) %>%
  mutate(HESdiagorder_prim = row_number()) %>% ungroup()

nonbrprimpre <- nonbrprimpre %>%
  select(bgs_id, HEScandat_prim, HESdiagsite_prim, HESdiagorder_prim, HESdiagsite_primc)

nonbrprimpre <- nonbrprimpre %>%
  pivot_wider(
    id_cols = bgs_id,
    names_from = HESdiagorder_prim,
    values_from = c(HEScandat_prim, HESdiagsite_prim, HESdiagsite_primc)
  )
head(nonbrprimpre)
# ------------------------------------------------------------------------------
# 3 Primary diagnosis initial breast

initbrprim <- HESlongflag %>%
  filter(HESdiag_nc3 == "C50")

n_distinct(initbrprim$bgs_id)

initbrprim <- left_join(initbrprim, eventdta, by = "bgs_id")

# Calculate time from randomisation
initbrprim <- initbrprim %>%
  mutate(tranddiag_n = (HESdat_n-DIAGNOSISDATEBEST)/7)

#ggplot(initbrprim, aes(y = tranddiag_n)) +
#  geom_boxplot() + labs(title = "Boc plot of tranddiag_n", y = "tranddiag_n") +
#  theme_minimal()

# Assumes reference to C50 within set time period is the initial diagnosis
freq(initbrprim$her2_Status)

initbrprim <- initbrprim %>%
  filter(!(HESdat_n > (DIAGNOSISDATEBEST + 7*wks) & !is.na(HESdat_n) & her2_Status == "Positive")) %>%
  filter(!(HESdat_n > (DIAGNOSISDATEBEST + 7*wks + 7*her2wks) & !is.na(HESdat_n) & her2_Status == "Negative")) %>%
  filter(!(HESdat_n > (DIAGNOSISDATEBEST + 7*wks + 7*her2wks) & !is.na(HESdat_n) & her2_Status == "Borderline"))

#ggplot(initbrprim, aes(y = tranddiag_n)) +
#  geom_boxplot() + labs(title = "Boc plot of tranddiag_n", y = "tranddiag_n") +
#  theme_minimal()

initbrprim <- initbrprim %>%
  arrange(bgs_id, HESdiag_nc4, HESdat_n) %>%
  distinct(bgs_id, HESdiag_nc4, .keep_all = TRUE) %>%
  distinct(bgs_id, HESdiag_nc3, .keep_all = TRUE)

initbrprim <- initbrprim %>%
  rename(
    HEScandat_brprim = HESdiag_nc4,
    HESdiagsite_brprim = HESdiag_nc3,
    HESdiagsite_brprimc = HESdiag_d3
  )

initbrprim <- initbrprim %>%
  group_by(bgs_id) %>%
  mutate(HESdiagorder_brprim = row_number()) %>%
  ungroup()

initbrprim <- initbrprim %>%
  select(bgs_id, HEScandat_brprim, HESdiagsite_brprim, HESdiagorder_brprim, HESdiagsite_brprimc)

n_distinct(initbrprim$bgs_id)
head(initbrprim)

# ------------------------------------------------------------------------------
# 4 Unknown lymph node information post initial diagnosis

lympunk <- HESlongflag %>%
  filter(HESdiag_nc4 == "C779")

n_distinct(lympunk$bgs_id)

# add randomisation date as only interested in 2nd primary and recurrence post this date

lympunk <- left_join(lympunk, eventdta, by = "bgs_id")

lympunk <- lympunk %>% 
  filter(HESdat_n > DIAGNOSISDATEBEST)

# assume if her2 status is unknown then it is negative as most common
freq(lympunk$her2_Status)

lympunk$her2_Status <- ifelse(is.na(lympunk$her2_Status), "Negative", lympunk$her2_Status)

lympunk <- lympunk %>%
  mutate(temptime = ifelse(HESdat_n <= (DIAGNOSISDATEBEST + 7 * wks) & 
                             !is.na(HESdat_n) & (her2_Status == "Positive"), 1, NA),
         temptime = ifelse(HESdat_n <= (DIAGNOSISDATEBEST + 7 * wks + 7 * her2wks) &
                             !is.na(HESdat_n) & (her2_Status == "Negative" | her2_Status == "Borderline"), 1, temptime),
         temptime = ifelse(HESdat_n > (DIAGNOSISDATEBEST + 7 * wks) &
                             !is.na(HESdat_n) & (her2_Status == "Positive"), 2, temptime),
         temptime = ifelse(HESdat_n > (DIAGNOSISDATEBEST + 7 * wks + 7 * her2wks) &
                             !is.na(HESdat_n) & (her2_Status == "Negative" | her2_Status == "Borderline"), 2, temptime))

freq(lympunk$temptime)

lympunk <- lympunk %>% filter(temptime == 2)  

lympunk <- lympunk %>%
  rename(HEScandat_lympunk = HESdat_n) %>%
  arrange(bgs_id, HEScandat_lympunk) %>%
  distinct(bgs_id, .keep_all = TRUE) %>%
  select(bgs_id, HEScandat_lympunk)

head(lympunk)
# ------------------------------------------------------------------------------
# 5 Unknown lymph node info initial diagnosis

lympunki <- HESlongflag %>%
  filter(HESdiag_nc4 == "C779")

n_distinct(lympunki$bgs_id)

# add randomisation date as only interested in 2nd primary and recurrence post this date
lympunki <- left_join(lympunki, eventdta, by = "bgs_id")

lympunki <- lympunki %>% 
  filter(HESdat_n > DIAGNOSISDATEBEST)

# assume if her2_Status unknown then it is negative as most common
freq(lympunki$her2_Status)
lympunki$her2_Status <- ifelse(is.na(lympunki$her2_Status), "Negative", lympunki$her2_Status)

lympunki <- lympunki %>%
  mutate(temptime = ifelse(HESdat_n <= (DIAGNOSISDATEBEST + 7 * wks) & 
                             !is.na(HESdat_n) & (her2_Status == "Positive"), 1, NA),
         temptime = ifelse(HESdat_n <= (DIAGNOSISDATEBEST + 7 * wks + 7 * her2wks) &
                             !is.na(HESdat_n) & (her2_Status == "Negative" | her2_Status == "Borderline"), 1, temptime),
         temptime = ifelse(HESdat_n > (DIAGNOSISDATEBEST + 7 * wks) &
                             !is.na(HESdat_n) & (her2_Status == "Positive"), 2, temptime),
         temptime = ifelse(HESdat_n > (DIAGNOSISDATEBEST + 7 * wks + 7 * her2wks) &
                             !is.na(HESdat_n) & (her2_Status == "Negative" | her2_Status == "Borderline"), 2, temptime))
freq(lympunki$temptime)

# unknown lymph node in temptime 1 considered to be local and initial diagnosis
lympunki <- lympunki %>% filter(temptime == 1)  

lympunki <- lympunki %>%
  rename(HEScandat_lympunki = HESdat_n) %>% 
  arrange(desc(HEScandat_lympunki), bgs_id) %>%
  distinct(bgs_id, .keep_all = TRUE) %>% 
  select(bgs_id, HEScandat_lympunki)

head(lympunki)

# ------------------------------------------------------------------------------
# 6 Lymph nodes local diagnosis - first post initial diagnosis

lympl <- HESlongflag %>%
  filter(HESdiag_nc4 == "C773")

n_distinct(lympl$bgs_id)

# add randomisation date as only interested in 2nd primary and recurrence post this date
lympl <- left_join(lympl, eventdta, by = "bgs_id")

lympl <- lympl %>% filter(HESdat_n > DIAGNOSISDATEBEST)

# assume if her2 status unknown then it is negative as most common
freq(lympl$her2_Status)
lympl$her2_Status <- ifelse(is.na(lympl$her2_Status), "Negative", lympl$her2_Status)

lympl <- lympl %>%
  mutate(temptime = ifelse(HESdat_n <= (DIAGNOSISDATEBEST + 7 * wks) & 
                             !is.na(HESdat_n) & (her2_Status == "Positive"), 1, NA),
         temptime = ifelse(HESdat_n <= (DIAGNOSISDATEBEST + 7 * wks + 7 * her2wks) &
                             !is.na(HESdat_n) & (her2_Status == "Negative" | her2_Status == "Borderline"), 1, temptime),
         temptime = ifelse(HESdat_n > (DIAGNOSISDATEBEST + 7 * wks) &
                             !is.na(HESdat_n) & (her2_Status == "Positive"), 2, temptime),
         temptime = ifelse(HESdat_n > (DIAGNOSISDATEBEST + 7 * wks + 7 * her2wks) &
                             !is.na(HESdat_n) & (her2_Status == "Negative" | her2_Status == "Borderline"), 2, temptime))
freq(lympl$temptime)

lympl <- lympl %>% filter(temptime == 2) 

lympl <- lympl %>%
  rename(HEScandat_lympl = HESdat_n) %>%
  arrange(bgs_id, HEScandat_lympl) %>%
  distinct(bgs_id, .keep_all = TRUE) %>%
  select(bgs_id, HEScandat_lympl)

head(lympl)

# ------------------------------------------------------------------------------
# 7 Lymph nodes distant diagnosis

lympd <- HESlongflag %>%
  filter((HESdiag_nn2 >= 77 & HESdiag_nn2 < 78) & HESdiag_nc1 == "C") %>%
  filter(HESdiag_nc4 !="C773" & HESdiag_nc4 !="C779")

n_distinct(lympd$bgs_id)

# add randomisation date as only interested in 2nd primary and recurrence post this date
lympd <- left_join(lympd, eventdta, by = "bgs_id")

lympd <- lympd %>% 
  filter(HESdat_n > DIAGNOSISDATEBEST)

n_distinct(lympd$bgs_id)

lympd <- lympd %>%
  distinct(bgs_id, HESdiag_nc4, .keep_all = TRUE) %>%
  arrange(bgs_id, HESdat_n) %>%
  rename(
    HEScandat_lympd = HESdat_n,
    HESdiagsite_lympd = HESdiag_nc4,
    HESdiagsite_lympdc = HESdiag_d4
  )

lympd <- lympd %>%
  group_by(bgs_id) %>%
  mutate(HESdiagorder_lympd = row_number()) %>%
  ungroup()

lympd <- lympd %>%
  select(bgs_id, HEScandat_lympd, HESdiagsite_lympd, HESdiagorder_lympd, HESdiagsite_lympdc)

lympd <- lympd %>%
  pivot_wider(
    id_cols = bgs_id,
    names_from = HESdiagorder_lympd, 
    values_from = c(HEScandat_lympd, HESdiagsite_lympd, HESdiagsite_lympdc))

head(lympd)

# ------------------------------------------------------------------------------
# 8 Lymph nodes diagnosis at initial diagnosis

# update if want to have separate fields than in 2nd primary and recurrence post this date
lympil <- HESlongflag %>%
  filter(HESdiag_nc4 == "C773")

# add randomisation date as only interested in 2nd primary and recurrence post this date
lympil <- left_join(lympil, eventdta, by = "bgs_id")

lympil <- lympil %>% 
  filter(HESdat_n > DIAGNOSISDATEBEST)

# assume if her2 status unknown then it is negative as most common
freq(lympil$her2_Status)
lympil$her2_Status <- ifelse(is.na(lympil$her2_Status), "Negative", lympil$her2_Status)

lympil <- lympil %>%
  mutate(temptime = ifelse(HESdat_n <= (DIAGNOSISDATEBEST + 7 * wks) & 
                             !is.na(HESdat_n) & (her2_Status == "Positive"), 1, NA),
         temptime = ifelse(HESdat_n <= (DIAGNOSISDATEBEST + 7 * wks + 7 * her2wks) &
                             !is.na(HESdat_n) & (her2_Status == "Negative" | her2_Status == "Borderline"), 1, temptime),
         temptime = ifelse(HESdat_n > (DIAGNOSISDATEBEST + 7 * wks) &
                             !is.na(HESdat_n) & (her2_Status == "Positive"), 2, temptime),
         temptime = ifelse(HESdat_n > (DIAGNOSISDATEBEST + 7 * wks + 7 * her2wks) &
                             !is.na(HESdat_n) & (her2_Status == "Negative" | her2_Status == "Borderline"), 2, temptime))

freq(lympil$temptime)

# initial diagnosis of lymph node involvement based on pre set time period
lympil <- lympil %>% filter(temptime == 1)

lympil <- lympil %>%
  rename(HEScandat_lympil = HESdat_n) %>%
  arrange(bgs_id, HEScandat_lympil) %>%
  distinct(bgs_id, .keep_all = TRUE) %>%
  select(bgs_id, HEScandat_lympil)

head(lympil)
# ------------------------------------------------------------------------------
# 9 Metastatic disease

sec <- HESlongflag %>%
  filter((HESdiag_nn2 >= 78 & HESdiag_nn2 <80) & HESdiag_nc1 == "C")

n_distinct(sec$bgs_id)

# secondary diagnosis when lymph nodes involved only
sec$temp <- as.integer(grepl("lymph", sec$HESoper_d4))

freq(sec$temp)

sec <- sec %>% filter(temp != 1)

n_distinct(sec$bgs_id)

# add randomisation date as only interested in 2nd primary and recurrence post this date
sec <- left_join(sec, eventdta, by = "bgs_id")

sec <- sec %>% filter(!(HESdat_n < DIAGNOSISDATEBEST))

n_distinct(sec$bgs_id)

# possible seen in different treatment specialties in one day but not considering this here
sec <- sec %>%
  group_by(bgs_id) %>%
  mutate(HESdrecno_n = n()) %>% ungroup()

sec <- sec %>%
  group_by(bgs_id) %>%
  mutate(tempmindat = min(HESdat_n),
         tempmaxdat = max(HESdat_n)) %>%
  ungroup()

class(sec$tempmindat)
class(sec$tempmaxdat)

sec <- sec %>%
  mutate(HESdrecrate_n =(tempmaxdat-tempmindat)/n()) %>%
  arrange(bgs_id, HESdiag_nc3, HESdat_n) %>%
  distinct(bgs_id, HESdiag_nc4, .keep_all = TRUE) %>%
  rename(
    HEScandat_sec = HESdat_n,
    HESdiagsite_sec = HESdiag_nc4,
    HESdiagsite_secc = HESdiag_d4
  ) %>%
  # create number different diagnosis sites
  arrange(bgs_id, HEScandat_sec, HESdiagsite_secc) %>%
  group_by(bgs_id) %>%
  mutate(HES_diagorder_sec = row_number()) %>% ungroup()

freq(sec$HES_diagorder_sec)
freq(sec$HESdiagsite_secc)
# shorten the text code

sec <- sec %>%
  mutate(temp_text = ifelse(grepl("secondary malignant neoplasm", HESdiagsite_secc), 1, NA),
         temp_text = ifelse(grepl("secondary malignant neoplasm of", HESdiagsite_secc), 2, temp_text),
         temp_len = nchar(HESdiagsite_secc))

freq(sec$temp_text)

sec <- sec %>%
  mutate(HESdiagsite_secc = ifelse(temp_text == 1, substr(HESdiagsite_secc, 31, temp_len), HESdiagsite_secc),
         HESdiagsite_secc = ifelse(temp_text == 2, substr(HESdiagsite_secc, 33, temp_len), HESdiagsite_secc))

freq(sec$HESdiagsite_secc)

sec <- sec %>%
  select(bgs_id, HEScandat_sec, HESdiagsite_sec, HES_diagorder_sec, HESdiagsite_secc,
         HESdrecno_n, HESdrecrate_n) # mainspef_nc, tretspef_nc


sec <- sec %>%
  pivot_wider(
    id_cols = c(bgs_id, HESdrecno_n, HESdrecrate_n),
    names_from = HES_diagorder_sec, 
    values_from = c(HEScandat_sec, HESdiagsite_sec, HESdiagsite_secc)) # mainspef_nc, tretspef_nc


head(sec)
n_distinct(sec$bgs_id)

# ------------------------------------------------------------------------------
# 10 Merge Data

head(nonbrprim)
head(initbrprim)
head(lympil)
head(lympl)
head(lympd)
head(lympunk)
head(lympunki)
head(sec)

df_list <- list(nonbrprim, initbrprim, lympil, lympl, lympd, lympunk, lympunki, sec)

HESwidediag_can <- df_list %>% reduce(full_join, by = 'bgs_id')
n_distinct(HESwidediag_can$bgs_id)


##########################
#### END OF CODE BLOCK ###
##########################