# ********************* CTSU Recurrence Algorithm  ************************


# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 26/06/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# ------------------------------------------------------------------------------

# Code converted from the stata script 3. Local rec data

# Purpose: Identify metastatic disease

# ------------------------------------------------------------------------------

metfull_local <- metfull %>%
  select(bgs_id, NCRASmet_dat_fin, NCRASmet_source_fin)

df_list_new3 <- list(excisionshort, metfull_local)

lrec <- df_list_new3 %>% reduce(full_join, by = 'bgs_id')
n_distinct(lrec$bgs_id) 

# creating field to include local and met information
lrec <- lrec %>%
  mutate(locrec_dat_n = 
           as.Date(case_when(
             EXCrec1_type1new == "locbr" | EXCrec1_type1new == "locbrl" ~ exdat_n_1,
             TRUE ~ NA)))

class(lrec$locrec_dat_n)
freq(lrec$EXCrec1_type1new)

lrec <- lrec %>%
  mutate(metlocNCRAS_diff = (NCRASmet_dat_fin - locrec_dat_n)/7,
         locrec_type_n = case_when(!is.na(locrec_dat_n) ~ "Loc",
                                   TRUE ~ NA),
         locrec_type_n = 
           case_when(!is.na(locrec_dat_n) & !is.na(NCRASmet_dat_fin) & metlocNCRAS_diff < 8 & metlocNCRAS_diff > -8 ~ "Loc+met",
                     TRUE ~ locrec_type_n))

# locrec_type_n = indicates if patient had metastatic disease date within 8 weeks +/- of local recurrence date

lrec <- lrec %>%
  select(bgs_id, locrec_dat_n, metlocNCRAS_diff, locrec_type_n) %>%
  filter(!is.na(locrec_dat_n))

############################
#### END OF CODE BLOCK #####
############################