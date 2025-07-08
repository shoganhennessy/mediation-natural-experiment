#!/usr/bin/R
## Senan Hogan-Hennessy, 14 June 2025.
## Script to to extract relevant data from Oregon Health Insurance rep package.
print(Sys.time())
set.seed(47)

## Packages:
# functions for data manipulation and visualisation
library(tidyverse)

# Define folder paths (1) input data (2) clean data.
data.folder <- file.path("..", "..", "data", "oregon-lottery-icspr")
figures.folder <- file.path("..", "..", "text", "sections", "figures")


################################################################################
## Load provided data.

# DS0001 Descriptive Variables
# -> Contains the 2008 lottery treatment variable, and other descriptives
descriptive.data <- data.folder %>%
    file.path("DS0001", "34314-0001-Data.dta") %>%
    haven::read_dta()

# DS0005 Twelve Month Mail Survey
# -> Contains info for >20k survey respondents, some useful outcomes.
survey_12m.data <- data.folder %>%
    file.path("DS0005", "34314-0005-Data.dta") %>%
    haven::read_dta()

# DS0006 In-person Survey
# -> COntains objective outcome measures.
inperson_12m.data <- data.folder %>%
    file.path("DS0006", "34314-0006-Data.dta") %>%
    haven::read_dta()


################################################################################
## Clean different data files.

# Restrict the basic descriptive data (for everyone in the entire data).
descriptive.data <- descriptive.data %>%
    transmute(person_id = person_id,
        household_id = household_id,
        # 0, 1 for the lottery
        lottery_iv = treatment) %>%
    # hh_size = how many people within each household.
    group_by(household_id) %>%
    mutate(hh_size = n()) %>%
    ungroup()

# Restrict the follow-up data to relevant variables.
survey_12m.data <- survey_12m.data %>%
    # Restrict to only those in the 12 month sample
    filter(sample_12m == 1 & returned_12m == 1 & !is.na(weight_12m)) %>%
    # Get the relevant variables.
    transmute(person_id = person_id,
        survey_weight = weight_12m,
        ## Variables for healthcare + insurance.
        any_insurance = ins_any_12m,
        usual_health_location = usual_care_12m,
        any_doc_visits = doc_any_12m,
        any_hospital_visits = er_any_12m,
        health_needs_met = needmet_med_cor_12m,
        ## Variables for (survey responded) outcomes.
        happiness_level_survey = happiness_12m, #: Current overall happiness
        health_level_survey = health_gen_12m, # Overall health
        health_change_survey = health_chg_12m, # How has your health changed: past 6 months
        days_bad_health_survey = baddays_phys_12m, # Number of days (out of past 30) when physical health not goo
        days_bad_mentalhealth_survey = baddays_ment_12m, # Number of days (out of past 30) when mental health not good
        health_limits_work = health_work_12m, # Physical, mental or emotional problem currently limits ability to work
        ## Controls for (presumably prior) health conditions.
        dia_diagnosis = dia_dx_12m, # Ever been told by a health professional that you have: Diabetes/Sugar diabetes
        ast_diagnosis = ast_dx_12m, # Ever been told by a health professional that you have: Asthma
        hbp_diagnosis = hbp_dx_12m, # Ever been told by a health professional that you have: High blood pressure
        emp_diagnosis = emp_dx_12m, # Ever been told by a health professional that you have: COP
        ami_diagnosis = ami_dx_12m, # Ever been told by a health professional that you have: Heart Disease/Angina/hear
        chf_diagnosis = chf_dx_12m, # Ever been told by a health professional that you have: Congestive Heart Failure
        dep_diagnosis = dep_dx_12m, # Ever been told by a health professional that you have: Depression or Anxiety
        chl_diagnosis = chl_dx_12m, # Ever been told by a health professional that you have: High Cholesterol
        kid_diagnosis = kid_dx_12m) # Ever been told by a health professional that you have: Kidney Problems

# Join the relevant data files
oregon.data <- descriptive.data %>%
    left_join(survey_12m.data, by = "person_id")

# Show the merging rates.
oregon.data %>% pull(health_needs_met) %>% table(exclude = NULL)


################################################################################
## Save the resulting (merged) data file.

# Restrict to only those who followed up in the survey.
oregon.data <- oregon.data %>%
    filter(!is.na(survey_weight))

# Save the file.
oregon.data %>% write_csv(file.path(data.folder, "cleaned-oregon-data.csv"))


#! Extra data, not yet used:
quit("no")
#! The in-person follow up data.
DS0006 In-person Survey
-> COntains objective outcome measures.

ast_dx_pre_lottery_inp: Diagnosed with asthma before the lottery (March 10th, 2008)
dia_dx_pre_lottery_inp: Diagnosed with diabetes before the lottery (March 10th, 2008)
hbp_dx_pre_lottery_inp: Diagnosed with hypertension before the lottery (March 10th, 2008)
chl_dx_pre_lottery_inp: Diagnosed with high cholesterol before the lottery (March 10th, 2008)
ami_dx_pre_lottery_inp: Diagnosed with heart attack before the lottery (March 10th, 2008)
chf_dx_pre_lottery_inp: Diagnosed with congestive heart failure before the lottery (March 10th, 2008)
emp_dx_pre_lottery_inp: Diagnosed with emphysema/COPD before the lottery (March 10th, 2008)
kid_dx_pre_lottery_inp: Diagnosed with kidney failure before the lottery (March 10th, 2008)
cancer_dx_pre_lottery_inp: Diagnosed with cancer before the lottery (March 10th, 2008)
dep_dx_pre_lottery_inp: Diagnosed with depression before the lottery (March 10th, 2008)
usual_clinic_inp: Indicated that usual place of care is a doctor office

phqtot_inp: PHQ total severity score
pcs8_score_inp: SF-8 physical component score
mcs8_score_inp: SF-8 mental component score
cvd_risk_point_inp: Cardiovascular Risk, Framingham risk score
doc_any_incl_probe_inp: Any doctor visits in the past 12 months
ed_num_incl_probe_inp: Number of ER visits in the past 12 months
ins_any_inp: Has any medical insurance at time of interview
ins_ohp_inp: Has OHP insurance
ins_private_inp: Has private insurance (either employer provided or privately purchased)
race_white_inp: Race/Ethnicity is White
race_black_inp: Race/Ethnicity is Black