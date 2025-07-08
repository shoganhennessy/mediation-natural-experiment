#!/usr/bin/R
## Senan Hogan-Hennessy, 2 July 2025.
## Script to to extract relevant data from causal claims data
print(Sys.time())
set.seed(47)

## Packages:
# functions for data manipulation and visualisation
library(tidyverse)

# Define folder paths (1) input data (2) clean data.
data.folder <- file.path("..", "..", "data", "causal-claims")
figures.folder <- file.path("..", "..", "text", "sections", "figures")
tables.folder <- file.path("..", "..", "text", "sections", "tables")


################################################################################
## Load provided data.

# Load the arrow library, for parquet files.
library(arrow)
# Read the causal claims data (from Parquet file).
causal.data <- data.folder %>%
    file.path("causal_claims_beta.parquet") %>%
    read_parquet()


################################################################################
## Clean data file, ahead of analysis.

# Since the causal revolution
causal.data <- causal.data %>%
    filter(year >= 2000)

# Restrict to empirical, causal papers.
causal.data <- causal.data %>%
    filter(classification_of_paper == "is_empirical_only" |
        classification_of_paper == "mostly_empirical")

# Restrict to NBER WPs
causal.data <- causal.data %>%
    filter(paper_repo == "nber")

# Restrict to papers that use an established causal inference method.
causal.data <- causal.data %>% filter(
    str_detect(tolower(sources_of_exogenous_variation), "random") |
    str_detect(tolower(sources_of_exogenous_variation), "exogenous") |
        paper_method_RDD == 1  |
        paper_method_DID == 1  |
        paper_method_RCT == 1  |
        paper_method_IV  == 1  |
        paper_method_TWFE == 1 |
        `paper_method_Event Study` == 1) %>%
    filter(
        paper_method_Structural != 1 &
        paper_method_Theoretical != 1)

# Fix column for whether it was published
causal.data <- causal.data %>%
    mutate(publication_outlet = ifelse(
        str_detect(tolower(publication_outlet), "nber") |
        str_detect(tolower(publication_outlet), "bureau"),
            NA, publication_outlet)) %>%
    mutate(published = !is.na(publication_outlet))

# Restrict to the relevant columns.
restricted.data <- causal.data %>%
    select(
        authors,
        title,
        year,                              
        classification_of_paper,
        paper_id,
        paper_author_meta,
        paper_year_meta,                   
        paper_abstract_meta,               
        paper_title_meta,
        paper_jel_coalesced,               
        edge_number,
        paper_edge_id,                     
        claim,
        cause,                             
        effect,
        type_of_causal_relationship,       
        is_evidence_provided_in_paper,
        sign_of_impact,                    
        effect_size,
        statistical_significance,          
        causal_inference_method,
        evidence_method_other_description, 
        sources_of_exogenous_variation,
        level_of_tentativeness,            
        is_tentative,                      
        is_certain,
        published)

# Restrict to only causal claims for mediated effects
restricted.data <- restricted.data %>%
    mutate(mediate_effect = as.integer(
        type_of_causal_relationship == "mediation")) %>%
    filter(!is.na(mediate_effect)) %>%
    group_by(paper_id, year) %>%
    mutate(mediate_effect = max(mediate_effect, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(mediate_effect == 1) %>%
    select(- mediate_effect)


################################################################################
## Count how many papers have a mediated effect, by year

# How many papers have mediation effects?  By year.
yearly_mediate.data <- causal.data %>%
    mutate(mediate_effect = as.integer(
        type_of_causal_relationship == "mediation")) %>%
    filter(!is.na(mediate_effect)) %>%
    group_by(paper_id, year) %>%
    summarise(mediate_effect = max(mediate_effect, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(general_count = 1, mediate_count = mediate_effect) %>%
    group_by(year) %>%
    summarise(
        general_count = sum(general_count, na.rm = TRUE),
        mediate_count = sum(mediate_count, na.rm = TRUE))


################################################################################
## Save the resulting data files.

# Save the data on the mediation effects from working papers as a data file.
restricted.data %>%
    write_csv(file.path(data.folder, "causal-wp-mediate.csv"))

# Save the yearly count as a data file.
yearly_mediate.data %>%
    write_csv(file.path(data.folder, "yearly-causal-mediate.csv"))
