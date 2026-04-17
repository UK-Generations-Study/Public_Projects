# ********************* CTSU Recurrence Algorithm  ************************


# Purpose:  ICR-code the CTSU STATA recurrence algorithm used to identify
# breast cancer recurrence events. 
# R implementation and translation of original Stata code.
# Original Stata code developed by: ICR-CTSU 
# ICR-CTSU contact: Lucy Kilburn (Lucy.Kilburn@icr.ac.uk)
# Date: 12/06/2024
# Imogen Sawyer (imogen.sawyer@icr.ac.uk)

# ------------------------------------------------------------------------------

# Code converted from the stata script 7. POETIC RT & SACT
# Purpose: Sorting out RT and SACT data

# ------------------------------------------------------------------------------
# 1 RT data
# ------------------------------------------------------------------------------

# (a) Formatting
# stata code uses decision to treat date, in cohort data we only have treatment start date
RTDSreg <- rtds %>%
  arrange(bgs_id, treatmentstartdate, proceduredate)


# tidying site field
RTDSreg$rttreatmentregion <- str_to_upper(RTDSreg$rttreatmentregion)
RTDSreg$rttreatmentregion <- gsub(" ", "", RTDSreg$rttreatmentregion)

RTDSreg$rttreatmentanatomicalsite <- str_to_upper(RTDSreg$rttreatmentanatomicalsite)

# codes with 0 at end don't exist so updating to ones that do (estimated)
RTDSreg <- RTDSreg %>%
  mutate(templen = nchar(rttreatmentanatomicalsite))

freq(RTDSreg$templen)

RTDSreg <- RTDSreg %>%
  mutate(templastc = ifelse(templen == 4, substr(rttreatmentanatomicalsite, 4, 4), NA),
         templastc = ifelse(templen == 3, substr(rttreatmentanatomicalsite, 3, 3), templastc))

freq(RTDSreg$templastc)

RTDSreg <- RTDSreg %>%
  mutate(rttreatmentanatomicalsite = ifelse(templastc == "0", 
                                            substr(rttreatmentanatomicalsite, 1, 3),
                                            rttreatmentanatomicalsite))
# some terms include text - remove text
RTDSreg <- RTDSreg %>%
  mutate(rttreatmentanatomicalsite = ifelse(templen == 10 | templen == 17, 
                                            substr(rttreatmentanatomicalsite, 1, 4),
                                            rttreatmentanatomicalsite))
# ------------------------------------------------------------------------------
# (b) linking operation codes

opcs4_RTDS <- opcs4
opcs4_RTDS$rttreatmentanatomicalsite <- opcs4_RTDS$HESoper_nc3
opcs4_RTDS$rttreatmentanatomicasite_fd <- opcs4_RTDS$HESoper_d4

opcs4_RTDS <- opcs4_RTDS %>% select(-HESoper_nc3, -HESoper_d4)

# merge opcs4 data table with the RTDSreg by rttreatmentantomicalsite
RTDSreg <- left_join(RTDSreg, opcs4_RTDS, by = "rttreatmentanatomicalsite")

# ------------------------------------------------------------------------------

# (c) create one row per regimen, identify regimen

RTDSreg <- RTDSreg %>%
  arrange(bgs_id, treatmentstartdate, proceduredate)

RTDSreg <- RTDSreg %>%
  group_by(bgs_id, treatmentstartdate) %>%
  mutate(RTDSprstdat_n = min(proceduredate),
         RTDSprednddat = max(proceduredate)) %>%
  ungroup() %>%
  mutate(difftemp = treatmentstartdate - RTDSprstdat_n)

class(RTDSreg$RTDSprstdat_n)
class(RTDSreg$RTDSprednddat)

RTDSreg <- RTDSreg %>%
  group_by(bgs_id, treatmentstartdate) %>%
  mutate(tempst = min(treatmentstartdate)) %>% ungroup() %>%
  mutate(RTDSprstdat_n = as.Date(ifelse(is.na(RTDSprstdat_n), tempst, RTDSprstdat_n)),
         tempdays = (prescribedfractions/5)*7,
         RTDSprednddat = as.Date(ifelse(is.na(RTDSprednddat), RTDSprstdat_n + tempdays, RTDSprednddat))) %>%
  arrange(bgs_id, rttreatmentanatomicalsite, treatmentstartdate, RTDSprstdat_n)

class(RTDSreg$RTDSprstdat_n)

# check whether site of treatment is missing within 'planned trt' or missing for complete planned treatment
RTDSreg <- RTDSreg %>%
  mutate(temp = ifelse(rttreatmentanatomicalsite == "" & 
                         lag(rttreatmentanatomicalsite) !="" & bgs_id == lag(bgs_id), 1, NA))

# identifying unique treatment regimens. Also possible to trt different sites
# during the same treatment regimen and different types (e.g. primary and mets)
# so also use site of treatment as identifier

RTDSreg <- RTDSreg %>%
  arrange(bgs_id, treatmentstartdate, RTDSprstdat_n, proceduredate) %>%
  distinct(bgs_id, rttreatmentanatomicalsite, treatmentstartdate, RTDSprstdat_n, 
           rttreatmentregion, .keep_all = TRUE)

# calculate the time from previous RTD
RTDSreg <- RTDSreg %>%
  arrange(bgs_id, RTDSprstdat_n) %>%
  mutate(tprevRT_n = ifelse(bgs_id == lag(bgs_id),
                            (RTDSprstdat_n-lag(RTDSprednddat))/7, NA))

ggplot(RTDSreg, aes(y = tprevRT_n)) +
  geom_boxplot() + labs(title = "Box plot of tprevRT_n", y = "tprevRT_n") +
  theme_minimal()

RTDSreg <- RTDSreg %>%
  arrange(bgs_id, treatmentstartdate) %>%
  mutate(RTDSsiteopcs_nc1 = substr(rttreatmentanatomicalsite, 1,1),
         RTDSsiteopcs_n3 = substr(rttreatmentanatomicalsite, 1,3),
         RTDSsiteopcs_nn2 = substr(rttreatmentanatomicalsite, 2,4),
         RTDSsiteopcs_nn2 = as.numeric(RTDSsiteopcs_nn2),
         # Breast
         RTDSsite_flag = ifelse(RTDSsiteopcs_n3 == "Z15", 1, NA),
         # Lymph nodes
         RTDSsite_flag = ifelse(RTDSsiteopcs_n3 == "Z61", 2, RTDSsite_flag),
         # Liver
         RTDSsite_flag = ifelse(rttreatmentanatomicalsite == "Z301", 3, RTDSsite_flag),
         # Lung
         RTDSsite_flag = ifelse(rttreatmentanatomicalsite == "Z246" | 
                                  RTDSsiteopcs_n3 == "Z52" |
                                  rttreatmentanatomicalsite == "Z245", 4, RTDSsite_flag),
         # Bladder
         RTDSsite_flag = ifelse(rttreatmentanatomicalsite == "Z421", 5, RTDSsite_flag),
         # Brain
         RTDSsite_flag = ifelse(RTDSsiteopcs_n3 == "Z01" | RTDSsiteopcs_n3 == "Z02" |
                                  RTDSsiteopcs_n3 == "Z03" | RTDSsiteopcs_n3 == "Z04" |
                                  RTDSsiteopcs_n3 == "Z05", 6, RTDSsite_flag),
         # Bone
         RTDSsite_flag = ifelse(RTDSsiteopcs_nc1 == "Z" & 
                                  (RTDSsiteopcs_nn2>=63 & RTDSsiteopcs_nn2 <=80), 7, RTDSsite_flag),
         # Skin
         RTDSsite_flag = ifelse(RTDSsiteopcs_n3 == "Z47" | RTDSsiteopcs_n3 == "Z48" | 
                                  RTDSsiteopcs_n3 == "Z49" | RTDSsiteopcs_n3 == "Z50", 8, RTDSsite_flag),
         # Other
         RTDSsite_flag = ifelse(is.na(RTDSsite_flag) & rttreatmentanatomicalsite !="", 100, RTDSsite_flag),
         # No Code
         RTDSsite_flag = ifelse(rttreatmentanatomicalsite == "", 101, RTDSsite_flag))

RTDSreg <- RTDSreg %>%
  mutate(RTDSsite_flags = case_when(
    RTDSsite_flag == 1 ~ "Breast",
    RTDSsite_flag == 2 ~ "Lymph nodes",
    RTDSsite_flag == 3 ~ "Liver",
    RTDSsite_flag == 4 ~ "Lung",
    RTDSsite_flag == 5 ~ "Bladder",
    RTDSsite_flag == 6 ~ "Brain",
    RTDSsite_flag == 7 ~ "Bone",
    RTDSsite_flag == 8 ~ "Skin",
    RTDSsite_flag == 100 ~ "Other",
    RTDSsite_flag == 101 ~ "No code",
    TRUE ~ NA
  ))

freq(RTDSreg$RTDSsite_flag)

RTDSreg <- RTDSreg %>%
  select(bgs_id, rttreatmentregion, RTDSprstdat_n, tprevRT_n, RTDSprednddat,
         rttreatmentanatomicalsite, RTDSsiteopcs_n3, RTDSsiteopcs_nn2, RTDSsite_flag, RTDSsite_flags)

# ------------------------------------------------------------------------------
# 2 SACT data
# ------------------------------------------------------------------------------

SACTreg <- sact %>%
  arrange(start_date_of_regimen, bgs_id) %>%
  arrange(bgs_id, start_date_of_regimen, start_date_of_cycle, cycle_number)

# calculate start and end of regimen and number of cycles
class(SACTreg$date_of_final_treatment)

SACTreg <- SACTreg %>%
  group_by(bgs_id, analysis_group, start_date_of_regimen) %>%
  mutate(regenddat_n = max(date_of_final_treatment)) %>% ungroup()

sum(is.na(SACTreg$regenddat_n))

# final date is missing a lot so substitute with last start of cycle
SACTreg <- SACTreg %>%
  group_by(bgs_id, analysis_group, start_date_of_regimen) %>%
  mutate(templastcycle = max(start_date_of_cycle)) %>%
  mutate(templastadmin = max(administration_date)) %>%
  ungroup()

class(SACTreg$templastadmin)
class(SACTreg$templastcycle)

# create field to indicate when estimate used
SACTreg <- SACTreg %>%
  mutate(regenddat_ne = ifelse(!is.na(regenddat_n), "regenddat_n", NA),
         regenddat_n = as.Date(ifelse(is.na(regenddat_n), templastcycle, regenddat_n)),
         regenddat_ne = ifelse(!is.na(regenddat_n) & is.na(regenddat_ne), 
                               "max start cycle", regenddat_ne),
         regenddat_n = as.Date(ifelse(is.na(regenddat_n), templastadmin, regenddat_n)),
         regenddat_ne = ifelse(!is.na(regenddat_n) & is.na(regenddat_ne), 
                               "max admin", regenddat_ne),
         regenddat_n = as.Date(ifelse(is.na(regenddat_n), start_date_of_regimen, regenddat_n)),
         regenddat_ne = ifelse(!is.na(regenddat_n) & is.na(regenddat_ne),
                               "start reg", regenddat_ne))

class(SACTreg$regenddat_n)
sum(is.na(SACTreg$regenddat_n))


SACTreg <- SACTreg %>%
  group_by(bgs_id, analysis_group, start_date_of_regimen) %>%
  mutate(totcycleno_n = max(cycle_number)) %>% ungroup()

# export data to check type of treatment
SACTreg <- SACTreg %>%
  distinct(bgs_id, analysis_group, start_date_of_regimen, .keep_all = TRUE)

freq(SACTreg$regenddat_ne)
freq(SACTreg$totcycleno_n)

# some treatments the same bench group and close together so remove if within 8 weeks
SACTreg <- SACTreg %>%
  mutate(tempclose = 
           ifelse((start_date_of_regimen - lag(start_date_of_regimen)) < 182 &
                    !is.na(start_date_of_regimen) & bgs_id == lag(bgs_id) &
                    benchmark_group == lag(benchmark_group), 1, NA),
         temp2 = 
           ifelse(bgs_id == lag(bgs_id) & benchmark_group == lag(benchmark_group),
                  (start_date_of_regimen - lag(start_date_of_regimen)), NA))

ggplot(SACTreg, aes(y = temp2)) +
  geom_boxplot() + labs(title = "Box plot of temp2", y = "temp2") +
  theme_minimal()

# remove data where within 26 weeks of previous treatment and benchgroup the same
freq(SACTreg$tempclose)

SACTreg <- SACTreg %>%
  filter(is.na(tempclose))

SACTreg <- SACTreg %>%
  arrange(bgs_id, start_date_of_regimen) %>%
  group_by(bgs_id) %>%
  mutate(regnoorder_n = row_number()) %>% ungroup()

# create a unique ID for each 'regimen'
SACTreg <- SACTreg %>%
  mutate(SACTuniqueidreg_n = row_number())

freq(SACTreg$regnoorder_n)

# calculate time from previous RTD
SACTreg <- SACTreg %>%
  arrange(bgs_id, start_date_of_regimen, start_date_of_cycle) %>%
  mutate(tprevCHEMO_n = ifelse(bgs_id == lag(bgs_id),
                               (start_date_of_regimen - lag(regenddat_n))/7, NA))

ggplot(SACTreg, aes(y = tprevCHEMO_n)) +
  geom_boxplot() + 
  labs(title = "Time from end previous chemo to start next (wks)", 
       y = "tprevCHEMO_n") +
  theme_minimal()

SACTreg <- SACTreg %>%
  select(bgs_id, analysis_group, benchmark_group, start_date_of_regimen, 
         regenddat_n, regenddat_ne, totcycleno_n, regnoorder_n, SACTuniqueidreg_n,
         tprevCHEMO_n)

freq(SACTreg$benchmark_group)

# ------------------------------------------------------------------------------

# Code SACT data

# use collapse when doing initial coding
# collapse (sum) CHEMOtrt_chemo, CHEMOtrt_VEGF, CHEMOtrt_aher2, CHEMOtrt_RT,
# CHEMOtrt_HM, CHEMOtrt_bone, by (benchmark_group)

SACTreg <- SACTreg %>%
  mutate(CHEMOtrt_chemo = 
           ifelse(grepl("CAPECITABINE", benchmark_group) |
                    grepl("CARBO", benchmark_group) |
                    grepl("CARBOPLATIN", benchmark_group) |
                    grepl("CHOP", benchmark_group) |
                    grepl("CISPLATIN", benchmark_group) |
                    grepl("CMF", benchmark_group) |
                    grepl("CTD", benchmark_group) |
                    grepl("CVD", benchmark_group) |
                    grepl("CYCLOPHOSPHAMIDE", benchmark_group) |
                    grepl("DOCETAXEL", benchmark_group) |
                    grepl("DOXORUBICIN", benchmark_group) |
                    grepl("EC", benchmark_group) |
                    grepl("EOX", benchmark_group) |
                    grepl("EPIRUBICIN", benchmark_group) |
                    grepl("ERLOTINIB", benchmark_group) |
                    grepl("EVEROLIMUS", benchmark_group) |
                    grepl("FCR", benchmark_group) |
                    grepl("FEC", benchmark_group) |
                    grepl("FLUOROURACIL", benchmark_group) |
                    grepl("GEMCARBO", benchmark_group) |
                    grepl("GEMCITABINE", benchmark_group) |
                    grepl("HYDROXYCARBAMIDE", benchmark_group) |
                    grepl("IMATINIB", benchmark_group) |
                    grepl("LIPOSOMAL DOXORUBICIN", benchmark_group) |
                    grepl("MITOMYCIN", benchmark_group) |
                    grepl("MITOMYCIN INTRAVESICULAR", benchmark_group) |
                    grepl("MITOXANTRONE", benchmark_group) |
                    grepl("MVP", benchmark_group) |
                    grepl("NAB-PACLITAXEL", benchmark_group) |
                    grepl("OXALIPLATIN + MDG", benchmark_group) |
                    grepl("PACLITAXEL", benchmark_group) |
                    grepl("TAC", benchmark_group) |
                    grepl("TCH", benchmark_group) |
                    grepl("TOPOTECAN", benchmark_group) |
                    grepl("VINORELBINE", benchmark_group) |
                    grepl("DENOSUMAB", benchmark_group) |
                    grepl("ERIBULIN", benchmark_group), 1, NA),
         )

freq(SACTreg$CHEMOtrt_chemo)

# added as just in POETIC but not IH
SACTreg <- SACTreg %>%
  mutate(CHEMOtrt_chemo = 
           ifelse(grepl("AC", benchmark_group) |
                    grepl("ALEMTUZUMAB + FLUDARABINE + MELPHALAN", benchmark_group) |
                    grepl("APHINITY TRIAL", benchmark_group) |
                    grepl("AXITINIB", benchmark_group) |
                    grepl("CHLORAMBUCIL + RITUXIMAB", benchmark_group) |
                    grepl("GEFITINIB", benchmark_group) |
                    grepl("LENALIDOMIDE", benchmark_group) |
                    grepl("NILOTINIB", benchmark_group) |
                    grepl("OXALIPLATIN + MDG", benchmark_group) |
                    grepl("PAZOPANIB", benchmark_group) |
                    grepl("PEGGY TRIAL", benchmark_group) |
                    grepl("PERSEPHONE TRIAL", benchmark_group) |
                    grepl("POETIC TRIAL", benchmark_group) |
                    grepl("POMALIDOMIDE", benchmark_group) |
                    grepl("RITUXIMAB", benchmark_group) |
                    grepl("SAFEHER TRIAL", benchmark_group) |
                    grepl("SUNITINIB", benchmark_group) |
                    grepl("THALIDOMIDE", benchmark_group) |
                    grepl("TRIAL", benchmark_group), 1, CHEMOtrt_chemo))

freq(SACTreg$CHEMOtrt_chemo)

SACTreg <- SACTreg %>%
  mutate(CHEMOtrt_VEGF = 
           ifelse(grepl("BEVACIZUMAB", benchmark_group), 1, NA),
         CHEMOtrt_aher2 = 
           ifelse(grepl("LAPATINIB", benchmark_group) |
                    grepl("TRASTUZUMAB", benchmark_group), 1, NA),
         CHEMOtrt_RT = 
           ifelse(grepl("RT", benchmark_group), 1, NA),
         CHEMOtrt_HM = 
           ifelse(grepl("HORMONES", benchmark_group), 1, NA),
         CHEMOtrt_bone = 
           ifelse(grepl("ZOLEDRONIC ACID", benchmark_group), 1, NA))

# adjuvant treatment
SACTreg <- SACTreg %>%
  mutate(CHEMOtrt_intent = 
           ifelse(grepl("DOCETAXEL + PERTUZUMAB + TRASTUZUMAB", benchmark_group) |
                    grepl("DOCETAXEL + TRASTUZUMAB", benchmark_group) |
                    grepl("EC", benchmark_group) |
                    grepl("EC + DOCETAXEL", benchmark_group) |
                    grepl("EC + DOCETAXEL + TRASTUZUMAB", benchmark_group) |
                    grepl("EPIRUBICIN + CMF", benchmark_group) |
                    grepl("FEC", benchmark_group) |
                    grepl("FEC + DOCETAXEL", benchmark_group) |
                    grepl("FEC + DOCETAXEL", benchmark_group) |
                    grepl("FEC + DOCETAXEL + TRASTUZUMAB", benchmark_group) |
                    grepl("PACLITAXEL + TRASTUZUMAB", benchmark_group) |
                    grepl("PERTUZUMAB + TRASTUZUMAB", benchmark_group) |
                    grepl("TAC", benchmark_group) |
                    grepl("TCH", benchmark_group) |
                    grepl("PERTUZUMAB", benchmark_group), 1, NA))

# identify metastic treatment
SACTreg <- SACTreg %>%
  mutate(CHEMOtrt_intent = 
           ifelse(grepl("ABVD", benchmark_group) |
                    grepl("BEVACIZUMAB", benchmark_group) |
                    grepl("LIPOSOMAL DOXORUBICIN", benchmark_group) |
                    grepl("BEVACIZUMAB + DOCETAXEL + FEC", benchmark_group) |
                    grepl("BEVACIZUMAB + PACLITAXEL", benchmark_group) |
                    grepl("CAPECITABINE", benchmark_group) |
                    grepl("CAPECITABINE + EPIRUBICIN + OXALIPLATIN", benchmark_group) |
                    grepl("CAPECITABINE + LAPATINIB", benchmark_group) |
                    grepl("CAPECITABINE + OXALIPLATIN", benchmark_group) |
                    grepl("CARBOPLATIN", benchmark_group) |
                    grepl("CARBOPLATIN + ETOPOSIDE", benchmark_group) |
                    grepl("CARBOPLATIN + VINORELBINE", benchmark_group) |
                    grepl("CAV", benchmark_group) |
                    grepl("CISPLATIN + ETOPOSIDE", benchmark_group) |
                    grepl("CISPLATIN + GEMCITABINE", benchmark_group) |
                    grepl("CISPLATIN + PEMETREXED", benchmark_group) |
                    grepl("CTD", benchmark_group) |
                    grepl("CYCLOPHOSPHAMIDE + VINORELBINE", benchmark_group) |
                    grepl("CYCLOPHOSPHAMIDE HIGH DOSE", benchmark_group) |
                    grepl("DENOSUMAB", benchmark_group) |
                    grepl("DOCETAXEL", benchmark_group) |
                    grepl("DOXORUBICIN", benchmark_group) |
                    grepl("EPIRUBICIN", benchmark_group) |
                    grepl("ERIBULIN", benchmark_group) |
                    grepl("EVEROLIMUS", benchmark_group) |
                    grepl("GEMCARBO", benchmark_group) |
                    grepl("GEMCITABINE", benchmark_group) |
                    grepl("GEMCITABINE + NAB-PACLITAXEL", benchmark_group) |
                    grepl("IMATINIB", benchmark_group) |
                    grepl("IRINOTECAN + MDG", benchmark_group) |
                    grepl("LIPOSOMAL DOXORUBICIN", benchmark_group) |
                    grepl("MITOXANTRONE + PACLITAXEL", benchmark_group) |
                    grepl("NAB-PACLITAXEL", benchmark_group) |
                    grepl("OXALIPLATIN + MDG", benchmark_group) |
                    grepl("PACLITAXEL", benchmark_group) |
                    grepl("PACLITAXEL + BEVACIZUMAB", benchmark_group) |
                    grepl("PAMIDRONATE", benchmark_group) |
                    grepl("PEMETREXED", benchmark_group) |
                    grepl("TOPOTECAN", benchmark_group) |
                    grepl("TRASTUZUMAB", benchmark_group) |
                    grepl("TRASTUZUMAB EMTANSINE", benchmark_group) |
                    grepl("VINORELBINE", benchmark_group) |
                    grepl("ZOLEDRONIC ACID", benchmark_group) |
                    grepl("CISPLATIN + FLUOROURACIL + RT", benchmark_group) |
                    grepl("CISPLATIN + RT", benchmark_group) |
                    grepl("CMF", benchmark_group) |
                    grepl("CVD", benchmark_group) |
                    grepl("CYCLOPHOSPHAMIDE", benchmark_group) |
                    grepl("CYCLOPHOSPHAMIDE + LENALIDOMIDE", benchmark_group), 2, CHEMOtrt_intent))

# more details on hormone data
# maybe mets
SACTreg <- SACTreg %>%
  mutate(CHEMOtrt_intent = 
           ifelse(grepl("HORMONES", benchmark_group), 3, CHEMOtrt_intent))

freq(SACTreg$CHEMOtrt_intent)

SACTreg <- SACTreg %>%
  arrange(bgs_id, start_date_of_regimen)

# make wide format
SACTreg <- SACTreg %>%
  select(bgs_id, analysis_group, benchmark_group, start_date_of_regimen,
         totcycleno_n, regnoorder_n, CHEMOtrt_chemo, CHEMOtrt_VEGF,
         CHEMOtrt_aher2, CHEMOtrt_RT, CHEMOtrt_HM, CHEMOtrt_bone, CHEMOtrt_intent)

SACTreg_metwide <- SACTreg %>%
  pivot_wider(id_cols = bgs_id,
              names_from = regnoorder_n,
              values_from = c(analysis_group, benchmark_group, start_date_of_regimen,
                              totcycleno_n, CHEMOtrt_chemo, CHEMOtrt_VEGF,
                              CHEMOtrt_aher2, CHEMOtrt_RT, CHEMOtrt_HM,
                              CHEMOtrt_bone, CHEMOtrt_intent))

SACTreg_met <- SACTreg %>%
  filter(CHEMOtrt_intent == 2)

n_distinct(SACTreg$bgs_id)

SACTreg_met <- SACTreg %>%
  select(bgs_id, analysis_group, benchmark_group, start_date_of_regimen, 
         totcycleno_n, regnoorder_n) %>%
  rename(
    analysis_group_nmet = analysis_group,
    benchmark_group_nmet = benchmark_group,
    start_date_of_regimen_nmet = start_date_of_regimen,
    totcycleno_nmet = totcycleno_n
  ) %>%
  group_by(bgs_id) %>%
  mutate(regnoorder_nmet = row_number()) %>% ungroup() %>%
  select(-regnoorder_n)

SACTreg_metwide <- SACTreg %>%
  group_by(bgs_id) %>%
  mutate(regnoorder_nmet = row_number()) %>% ungroup() %>%
  filter(CHEMOtrt_intent == 2) %>%
  distinct(bgs_id, .keep_all = TRUE) %>%
  select(bgs_id, analysis_group, benchmark_group, start_date_of_regimen, 
         totcycleno_n, regnoorder_nmet) %>%
  pivot_wider(id_cols = bgs_id,
              names_from = regnoorder_nmet,
              values_from = c(analysis_group, benchmark_group, start_date_of_regimen,
                              totcycleno_n))

##########################
#### END OF CODE BLOCK ###
##########################