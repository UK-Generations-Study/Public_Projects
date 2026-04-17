# ********************* CTSU Recurrence Algorithm  ************************

# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 15/05/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# preparation of the death/vitalstatus data - determine breast cancer deaths

# ------------------------------------------------------------------------------

# Creating metastatic cancer diagnosis field

n_distinct(death_data$BGS_ID)

# dth_flagn = NCRAS death flag
death_data <- death_data %>%
  mutate(dth_flagn = ifelse(str_starts(VITALSTATUS, "D") & !is.na(VITALSTATUSDATE), 1, 0),
         deathdatebest = as.Date(ifelse(dth_flagn == 1, VITALSTATUSDATE, NA)))

freq(death_data$dth_flagn)
class(death_data$deathdatebest)


# breast cancer death based on codes -------------------------------------------

death_data$temp <- 0

varlist <- c("DEATHCAUSECODE_1A", "DEATHCAUSE_1B", "DEATHCAUSECODE_1C",
             "DEATHCAUSECODE_2", "DEATHCAUSECODE_UNDERLYING")

for (v in varlist) {
  death_data$temp <- death_data$temp + 1
  
  dthcan_n <- paste0("dthcan_n", death_data$temp[1])
  dthcanbr_n <- paste0("dthcanbr_n", death_data$temp[1])
  temp_C80 <- paste0("temp_C80", death_data$temp[1])
  
  death_data[[dthcan_n]] <- ""
  death_data[[dthcanbr_n]] <- ""
  death_data[[temp_C80]] <- ""
  
  death_data[[dthcan_n]][grepl("C", death_data[[v]], ignore.case = TRUE) | 
                           grepl("c", death_data[[v]], ignore.case = TRUE)] <- "Cancer"
  death_data[[dthcanbr_n]][grepl("C50", death_data[[v]], ignore.case = TRUE) | 
                             grepl("c50", death_data[[v]], ignore.case = TRUE)] <- "Breast"
}

freq(death_data$dthcan_n1)
freq(death_data$dthcanbr_n1)

# breast cancer deaths based on text/descriptions ------------------------------

death_data <- death_data %>% select(-temp)

death_data$temp <- 0

varlist2 <- c("DEATHCAUSETEXT_1A", "DEATHCAUSETEXT_1B", "DEATHCAUSETEXT_1C",
              "DEATHCAUSETEXT_2", "DEATHCAUSETEXT_UNDERLYING")

for (v in varlist2) {
  death_data$temp <- death_data$temp + 1
  
  dthcan_nc <- paste0("dthcan_nc", death_data$temp[1])
  
  death_data[[dthcan_nc]][grepl("cancer", death_data[[v]], ignore.case = TRUE) |
                            grepl("carcinoma", death_data[[v]], ignore.case = TRUE) |
                            grepl("malignancy", death_data[[v]], ignore.case = TRUE) |
                            grepl("malignant", death_data[[v]], ignore.case = TRUE)] <- "Cancer"
  
  dthcanbr_nc <- paste0("dthcanbr_nc", death_data$temp[1])
  
  death_data[[dthcanbr_nc]][grepl("breast", death_data[[v]], ignore.case = TRUE) & 
                              death_data[[dthcan_nc]] == "Cancer"] <- "Breast"
}

freq(death_data$dthcanbr_nc1)
freq(death_data$dthcan_nc1)

# overall - per patient summary ------------------------------------------------

death_data <- death_data %>%
  mutate(dthcan_n = ifelse(
    (dthcan_n1 == "Cancer" | dthcan_n2 == "Cancer" | dthcan_n3 == "Cancer" |
      dthcan_n4 == "Cancer" | dthcan_n5 == "Cancer"), "Cancer", NA),
    dthcanbr_nn = ifelse(
      (dthcanbr_n1 == "Breast" | dthcanbr_n2 == "Breast" | dthcanbr_n3 == "Breast" |
        dthcanbr_n4 == "Breast" | dthcanbr_n5 == "Breast"), "Breast", NA),
    dthcanbr_nc = ifelse(
      (dthcanbr_nc1 == "Breast" | dthcanbr_nc2 == "Breast" | dthcanbr_nc3 == "Breast" |
        dthcanbr_nc4 == "Breast" | dthcanbr_nc5 == "Breast"), "Breast", NA)
  )

freq(death_data$dthcan_n) # <- cancer diagnosis
freq(death_data$dthcanbr_nn)
freq(death_data$dthcanbr_nc)

# sometimes text says breast cancer but codes don't
death_data <- death_data %>%
  mutate(dthcanbr_n = ifelse(dthcanbr_nn == "Breast" | dthcanbr_nc == "Breast", "Breast", NA))

freq(death_data$dthcanbr_n) # <- Breast cancer diagnosis

death_data <- death_data %>%
  mutate(dthcan_not2n = 
           ifelse((dthcan_n1 == "Cancer" | dthcan_n2 == "Cancer" | dthcan_n3 == "Cancer") | 
                    (dthcan_nc1 == "Cancer" | dthcan_nc2 == "Cancer" | dthcan_nc3 == "Cancer"), "Cancer", NA),
         dthcanbr_not2n = 
           ifelse((dthcanbr_n1 == "Breast" | dthcanbr_n2 == "Breast" | dthcanbr_n3 == "Breast") |
                    (dthcanbr_nc1 == "Breast" | dthcanbr_nc2 == "Breast" | dthcanbr_nc3 == "Breast"), "Breast", NA))

freq(death_data$dthcan_not2n) # <- Cancer mentioned 1a to 1c
freq(death_data$dthcanbr_not2n) # <- Breast mentioned 1a to 1c

# Metastatic disease death based on text/description ---------------------------

# metastatic disease codes not suitable generally if main cause breast cancer
# then metastatic disease just C509 code used

death_data <- death_data %>% select(-temp)

death_data$temp <- 0

death_data <- death_data %>%
  mutate_at(vars(DEATHCAUSETEXT_1A, DEATHCAUSETEXT_1B, DEATHCAUSETEXT_1C,
                 DEATHCAUSETEXT_2, DEATHCAUSETEXT_UNDERLYING), as.character)

varlist2 <- c("DEATHCAUSETEXT_1A", "DEATHCAUSETEXT_1B", "DEATHCAUSETEXT_1C",
              "DEATHCAUSETEXT_2", "DEATHCAUSETEXT_UNDERLYING")


for (v in varlist2) {
  
  death_data$temp <- death_data$temp + 1
  
  dthmet_n <- paste0("dthmet_n", death_data$temp[1])
  
  death_data[[dthmet_n]][grepl("metastatic", death_data[[v]], ignore.case = TRUE) | 
                            grepl("metastases", death_data[[v]], ignore.case = TRUE) |
                            grepl("metastasis", death_data[[v]], ignore.case = TRUE) |
                            grepl("secondaries", death_data[[v]], ignore.case = TRUE) |
                            grepl("carcinomatosis", death_data[[v]], ignore.case = TRUE) |
                            grepl("progressive", death_data[[v]], ignore.case = TRUE) |
                            grepl("advanced", death_data[[v]], ignore.case = TRUE) |
                            grepl("disseminated", death_data[[v]], ignore.case = TRUE) |
                            grepl("C79", death_data[[v]], ignore.case = TRUE) |
                            grepl("dissemminated", death_data[[v]], ignore.case = TRUE) |
                            grepl("terminal", death_data[[v]], ignore.case = TRUE)] <- "OC mets"
  
  death_data[[dthmet_n]][(grepl("metastatic", death_data[[v]], ignore.case = TRUE) |
                             grepl("metastases", death_data[[v]], ignore.case = TRUE) |
                             grepl("metastasis", death_data[[v]], ignore.case = TRUE) |
                             grepl("secondaries", death_data[[v]], ignore.case = TRUE) |
                             grepl("carcinomatosis", death_data[[v]], ignore.case = TRUE) |
                             grepl("progressive", death_data[[v]], ignore.case = TRUE) |
                             grepl("advanced", death_data[[v]], ignore.case = TRUE) |
                             grepl("disseminated", death_data[[v]], ignore.case = TRUE) |
                             grepl("dissemminated", death_data[[v]], ignore.case = TRUE) |
                             grepl("terminal", death_data[[v]], ignore.case = TRUE) |
                             grepl("C79", death_data[[v]], ignore.case = TRUE)) &
                            grepl("breast", death_data[[v]], ignore.case = TRUE) |
                            grepl("c50", death_data[[v]], ignore.case = TRUE)] <- "BC mets"
  
}

freq(death_data$dthmet_n1)
freq(death_data$dthmet_n2)

# Metastatic disease based on text ---------------------------------------------

# dthmet_n = metastatic cancer diagnosis
death_data <- death_data %>%
  mutate(dthmet_n = ifelse(dthmet_n1 == "BC mets" | dthmet_n2 == "BC mets" |
                             dthmet_n3 == "BC mets" | dthmet_n4 == "BC mets" |
                             dthmet_n5 == "BC mets", "BC mets", NA))

death_data <- death_data %>%
  mutate(dthmet_n = ifelse((dthmet_n1 == "OC mets" | dthmet_n2 == "OC mets" |
                             dthmet_n3 == "OC mets" | dthmet_n4 == "OC mets" |
                             dthmet_n5 == "OC mets") & is.na(dthmet_n), "OC mets", dthmet_n))

freq(death_data$dthmet_n)

# cancer based underlying cause possible met disease ---------------------------

death_data <- death_data %>%
  # malignant neoplasm mentioned in underlying cause
  mutate(undcan_n = 
           ifelse(grepl("malignant neoplasm", DEATHCAUSETEXT_UNDERLYING), 1, NA),
         # breast mentioned in underlying cause
         undcanbr_n =
           ifelse(grepl("breast", DEATHCAUSETEXT_UNDERLYING) & undcan_n == 1, 1, NA),
         # unknown cancer
         undcanunk_n = 
           ifelse((grepl("unknown", DEATHCAUSETEXT_UNDERLYING) | 
                    grepl("without specification", DEATHCAUSETEXT_UNDERLYING) |
                    grepl("unspecification", DEATHCAUSETEXT_UNDERLYING) | 
                    grepl("unspecified", DEATHCAUSETEXT_UNDERLYING)) &
                    grepl("site", DEATHCAUSETEXT_UNDERLYING) &
                    undcan_n == 1, 1, NA))

freq(death_data$undcan_n)
freq(death_data$undcanbr_n)
freq(death_data$undcanunk_n)

# other cancer non-breast ------------------------------------------------------

death_data$dth2prim_n <- NA
death_data$dth2prim_n[death_data$dthcan_n == "Cancer" & is.na(death_data$dthcanbr_n)] <- "OC"
death_data$dth2prim_n[death_data$undcanunk_n == 1] <- "UNK"

freq(death_data$dth2prim_n) # <- non breast cancer diagnosis

# create combined death details - mets or underlying cause ---------------------

death_data <- death_data %>%
  # death BC mets and 1a and 1c is breast cancer
  mutate(dthcanbr_combn = 
           ifelse(dthmet_n == "BC mets" | dthcanbr_not2n == "Breast", 1, NA))

freq(death_data$dthcanbr_combn)

# breast cancer death and not breast cancer death ------------------------------

death_data <- death_data %>%
  # nreast cancer death - NCRAS
  mutate(dthdat_BC_n = as.Date(ifelse(dthcanbr_combn == 1, deathdatebest, NA)),
         dthdat_nonBC_n = as.Date(ifelse(is.na(dthdat_BC_n), deathdatebest, NA)))

class(death_data$dthdat_BC_n) # <- breast cancer death NCRAS

# checking data - update based on review ---------------------------------------

# if breast cancer and mets mentioned and breast mentioned then update to 
# metastatic breast cancer once checked manually

death_data <- death_data %>%
  mutate(dthmet_n = ifelse(dthmet_n == "OC mets" & dthcanbr_combn == 1, "BC mets", dthmet_n),
         # updating those that have indication of other primary but no met disease
         dth2prim_n = ifelse(is.na(dth2prim_n) & undcan_n == 1 & undcanbr_n !=1, "OC", dth2prim_n))

freq(death_data$dthmet_n)
freq(death_data$dth2prim_n)

# ------------------------------------------------------------------------------

# format data
deathwide_full <- death_data %>%
  select(BGS_ID, deathdatebest, DEATHCAUSECODE_1A, 
         DEATHCAUSECODE_1B, DEATHCAUSECODE_1C,
         DEATHCAUSECODE_2, DEATHCAUSECODE_UNDERLYING, DEATHCAUSETEXT_1A,
         DEATHCAUSETEXT_1B, DEATHCAUSETEXT_1C, DEATHCAUSETEXT_2,
         DEATHCAUSETEXT_UNDERLYING, dthcan_n, dthcanbr_n, dthmet_n, undcan_n,
         undcanbr_n, dth2prim_n, dth_flagn, undcanunk_n, dthcanbr_combn,
         dthdat_nonBC_n, dthdat_BC_n) %>%
  filter(!is.na(deathdatebest))
deathwide_full

n_distinct(deathwide_full$BGS_ID)


############################
#### END OF CODE BLOCK #####
############################