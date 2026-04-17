# ********************* CTSU Recurrence Algorithm  ************************

# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 15/05/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# Prepare datasets used in CTSU recurrence algorithm

# ------------------------------------------------------------------------------
# ICD-10
# create two columns in icd10 one for icd10 code and second for icd10 description

icd10 <- separate(icd10, text, into = c("first", "second"), sep = " ", remove = FALSE)

icd10 <- icd10 %>%
  select(-second) %>%
  mutate(dashflag = ifelse(grepl("-", first), 1, NA)) 

icd10 <- icd10 %>% filter(is.na(dashflag)) %>% select(-dashflag)

# ------------------------------------------------------------------------------
# Death data
# BGS_ID = unqiue participant identifier in the Breast Cancer Now Generations Study
n_distinct(end_status$BGS_ID)

end_status <- end_status %>%
  mutate_at(c('VITALSTATUSDATE'), dmy) %>%
  mutate_at(c('BGS_ID'), as.character)
end_status

death_dataL <- end_status %>%
  pivot_longer(cols = starts_with("DEATHCAUSECODE_"),
               names_to = "DEATHCAUSECODE",
               values_to = "ICD10")

death_dataL <- death_dataL %>%
  mutate(first = substr(ICD10, 1, 3))

death_dataL <- left_join(death_dataL, icd10, by = "first")
death_dataL <- death_dataL %>%
  select(-first) %>%
  mutate(text = tolower(text)) %>%
  rename(
    icd10text = text
  )

death_dataW <- death_dataL %>%
  pivot_wider(names_from = "DEATHCAUSECODE", values_from = c("ICD10", "icd10text"))
death_dataW

death_data <- death_dataW %>%
  rename(
    DEATHCAUSECODE_1A = ICD10_DEATHCAUSECODE_1A,
    DEATHCAUSECODE_1B = ICD10_DEATHCAUSECODE_1B,
    DEATHCAUSECODE_1C = ICD10_DEATHCAUSECODE_1C,
    DEATHCAUSECODE_2 = ICD10_DEATHCAUSECODE_2,
    DEATHCAUSECODE_UNDERLYING = ICD10_DEATHCAUSECODE_UNDERLYING,
    DEATHCAUSETEXT_1A = icd10text_DEATHCAUSECODE_1A,
    DEATHCAUSETEXT_1B = icd10text_DEATHCAUSECODE_1B,
    DEATHCAUSETEXT_1C = icd10text_DEATHCAUSECODE_1C,
    DEATHCAUSETEXT_2 = icd10text_DEATHCAUSECODE_2,
    DEATHCAUSETEXT_UNDERLYING = icd10text_DEATHCAUSECODE_UNDERLYING
  )

rm(death_dataW, death_dataL)

# ------------------------------------------------------------------------------
# OPCS 4
# create two columns in opcs4 one for opcs4 code and second for opcs4 description
opcs4 <- separate(opcs4, text, into = c("first", "second"), sep = " ", remove = FALSE)

opcs4$opcs4code <- opcs4$first
opcs4$opcs4text <- opcs4$text

opcs4 <- opcs4 %>%
  select(-first, -second, -text) %>%
  mutate(dashflag = ifelse(grepl("-", opcs4code), 1, NA)) 

opcs4 <- opcs4 %>% filter(is.na(dashflag)) %>% select(-dashflag)

# ------------------------------------------------------------------------------
# HES APC
hes_APC <- hes_APC %>%
  mutate_at(vars(starts_with("OPDATE"), starts_with("ADMIDATE")), dmy) %>%
  mutate_at(c('BGS_ID'), as.character)
hes_APC

names(hes_APC) <- tolower(names(hes_APC))
hes_APC

# ------------------------------------------------------------------------------
# HES OP
hes_OP <- hes_OP %>%
  mutate_at(c('apptdate'), dmy) %>%
  mutate_at(c('bgs_id'), as.character)
hes_OP

# ------------------------------------------------------------------------------
# SACT
# convert to lower case
names(sact) <- tolower(names(sact))

sact <- sact %>%
  mutate_at(c('start_date_of_regimen', 'start_date_of_cycle', 
              'administration_date', 'date_of_final_treatment'), dmy) %>%
  mutate_at(c('bgs_id', 'tumourid', 'administration_route',
              'perf_stat_start_of_cycle', 'perf_stat_start_of_regimen'), as.character)


# ------------------------------------------------------------------------------
# RTDS
# convert to lower case
names(rtds) <- tolower(names(rtds))

# remove those with no data
rtds <- rtds %>%
  filter(!is.na(radiotherapyepisodeid))

rtds <- rtds %>%
  mutate_at(c('proceduredate', 'treatmentstartdate', 'apptdate'), dmy) %>%
  mutate_at(c('attendid', 'radiotherapyepisodeid', 'bgs_id', 'radiotherapyintent'), as.character)

# ------------------------------------------------------------------------------
# Tumour PHE
names(tumour_phe) <- tolower(names(tumour_phe))

tumour_phe <- tumour_phe %>%
  mutate_at(c('diagnosisdatebest'), dmy) %>%
  mutate_at(c('bgs_id', 'tumourid'), as.character)
tumour_phe

# ------------------------------------------------------------------------------
# Treatment

names(treatment) <- tolower(names(treatment))

treatment <- treatment %>%
  mutate_at(c('eventdate'), dmy) %>%
  mutate_at(c('bgs_id', 'tumourid'), as.character)
treatment

############################
#### END OF CODE BLOCK #####
############################