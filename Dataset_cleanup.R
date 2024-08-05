# source("Genomes_labeling.R")
# Update some labels/remove some entries after another assessment (difficult labeling)
# Load necessary libraries
library(dplyr)
# Import data
# assess_differ_label <- read.delim("assess_differ_label.csv", sep=";", header=TRUE)
# Remove entries labeled as "Remove"
rem <- assess_differ_label[assess_differ_label$Final_label_v2 == "Remove", "genome_id"]
Genomes_labeled_v2 <- Genomes_labeled[!Genomes_labeled$genome_id %in% rem, ]

# Modify labels based on "Final_label_v2"
Genomes_labeled_v2 <- merge(Genomes_labeled_v2, assess_differ_label[,c("genome_id","Final_label_v2")], by = "genome_id", all.x = TRUE)
Label_v2 <- ifelse(!is.na(Genomes_labeled_v2$Final_label_v2), Genomes_labeled_v2$Final_label_v2, Genomes_labeled_v2$Label)
Genomes_labeled_v2$Label_v2 <- Label_v2

# Exclude non-informative variables
Genomes_labeled_edit <- Genomes_labeled_v2[, !names(Genomes_labeled_v2) %in% c("pathovar_ncbi", "organism_name", "genome_status", "sequencing_status", "sequences", "brc1_cds", "isolation_site", "latitude", "longitude", "altitude", "depth", "body_sample_subsite", "sporulation", "salinity", "groups", "sample_type", "environment (biome)", "env_feature", "metagenome_source", "Isolation_source", "Patient Sex", "Patient Age", "disease_ncbi", "Result", "Final_label_v2")]

# Merge multiple columns into new columns, eliminating "NA" and "missing" values
merge_multiple_unique <- function(vectors) {
  result <- character(0)
  merge_two_unique <- function(vec1, vec2) {
    mapply(function(x, y) {
      x_norm <- tolower(trimws(ifelse(is.na(x) | x == "missing" , "", x)))
      y_norm <- tolower(trimws(ifelse(is.na(y) | y == "missing", "", y)))
      if (x_norm != y_norm && x_norm != "" && y_norm != "") {
        paste(x, y, sep = "; ")
      } else if (x_norm != "") {
        x
      } else {
        y
      }
    }, vec1, vec2)
  }
  for (vec in vectors) {
    len_diff <- length(vec) - length(result)
    if (len_diff > 0) {
      result <- c(result, rep("", len_diff))
    } else if (len_diff < 0) {
      vec <- c(vec, rep("", -len_diff))
    }
    result <- merge_two_unique(result, vec)
  }
  return(unlist(result, use.names = FALSE))
}

# Merge columns for "isolation source"
isol_source_vectors <- list(Genomes_labeled_edit$isolation_source, Genomes_labeled_edit$`isolation source`, Genomes_labeled_edit$`isolation-source`, Genomes_labeled_edit$isolation_source_ncbi)
isolation_source <- merge_multiple_unique(isol_source_vectors)
Genomes_labeled_edit$isolation_source <- isolation_source
Genomes_labeled_edit <- Genomes_labeled_edit[, !names(Genomes_labeled_edit) %in% c("isolation source", "isolation-source","isolation_source_ncbi")]

# Merge columns for "host health state"
health_state_vectors <- list(Genomes_labeled_edit$`host health state`, Genomes_labeled_edit$health_state, Genomes_labeled_edit$host_health_state, Genomes_labeled_edit$health_disease_stat, Genomes_labeled_edit$host_health_ncbi)
health_state <- merge_multiple_unique(health_state_vectors)
Genomes_labeled_edit$host_health_state <- health_state
Genomes_labeled_edit <- Genomes_labeled_edit[, !names(Genomes_labeled_edit) %in% c("health_state", "host health state","health_disease_stat", "host_health_ncbi", "health_state_2", "merged_health", "host_disease_stat")]

# Merge columns for "disease"
disease_vectors <- list(Genomes_labeled_edit$disease, Genomes_labeled_edit$`host disease`, Genomes_labeled_edit$`Host disease`)
disease <- merge_multiple_unique(disease_vectors)
Genomes_labeled_edit$disease <- disease
Genomes_labeled_edit <- Genomes_labeled_edit[, !names(Genomes_labeled_edit) %in% c("host disease", "Host disease")]

# Merge columns for "biotic relationship"
biotic_vectors <- list(Genomes_labeled_edit$`observed biotic relationship`, Genomes_labeled_edit$biotic_relationship)
biotic_relationship <- merge_multiple_unique(biotic_vectors)
Genomes_labeled_edit$biotic_relationship <- biotic_relationship
Genomes_labeled_edit <- Genomes_labeled_edit[, !names(Genomes_labeled_edit) %in% "observed biotic relationship"]

# Merge columns for "biome"
biome_vectors <- list(Genomes_labeled_edit$biome, Genomes_labeled_edit$env_biome)
biome <- merge_multiple_unique(biome_vectors)
Genomes_labeled_edit$biome <- biome
Genomes_labeled_edit <- Genomes_labeled_edit[, !names(Genomes_labeled_edit) %in% c("env_biome")]

# Merge columns for "host status"
host_status_vectors <- list(Genomes_labeled_edit$host_status, Genomes_labeled_edit$`host status`)
host_status <- merge_multiple_unique(host_status_vectors)
Genomes_labeled_edit$host_status <- host_status
Genomes_labeled_edit <- Genomes_labeled_edit[, !names(Genomes_labeled_edit) %in% c("host status")]

# Merge columns for "subsource note"
subsource_note_vectors <- list(Genomes_labeled_edit$subsrc_note, Genomes_labeled_edit$subsource_note)
subsource_note <- merge_multiple_unique(subsource_note_vectors)
Genomes_labeled_edit$subsource_note <- subsource_note
Genomes_labeled_edit <- Genomes_labeled_edit[, !names(Genomes_labeled_edit) %in% c("subsrc_note")]

## 2 - If "host_disease" is not empty then it means that the bacterium caused the disease - should we relabel 
## inspect ones that were left out
left_out_labeling <- all_metadata_filt[(!all_metadata_filt$genome_id %in% Genomes_labeled_edit$genome_id ) ,] #5824; 5657 (v2); 5648
#x <- left_out_labeling[, c("genome_id","fields_pasted","host health state", "health_state","host_health_state", "health_disease_stat", "host_health_ncbi", "host_health", "host_health_ncbi", "host status", "host disease", "host_disease_stage", "host_status", "Host disease"),]
# inspect with "table"
## ONLY "host_health" / "host_disease" have a considerable amount of entries (3185), and many of them unique
left_out_labeling_with_host_dis <- left_out_labeling[!is.na(left_out_labeling$host_disease),] # 2274; 2950; 2940
# assess which ones are not related with infection
t <- table(left_out_labeling_with_host_dis$host_disease)
#write.table(t, "TERMS_left_out_labeling_with_host_dis.csv", sep="\t", row.names = FALSE, quote = FALSE) 
# Import assessed df
TERMS_left_out_labeling_with_host_dis <- read.delim("TERMS_left_out_labeling_with_host_dis.tsv", skip = 1, header = TRUE)
# 428
# select ones with "yes" as a vector 
yes_rows <- TERMS_left_out_labeling_with_host_dis[TERMS_left_out_labeling_with_host_dis$Infection_associated == "Yes", ] # 143 

# Process left-out labeling
left_out_labeling <- all_metadata_filt[!all_metadata_filt$genome_id %in% Genomes_labeled_edit$genome_id,]
left_out_labeling_with_host_dis <- left_out_labeling[!is.na(left_out_labeling$host_disease),]
terms_HP_host_dis <- yes_rows$Var1
HP_left_out_labeling_with_host_dis <- left_out_labeling_with_host_dis[left_out_labeling_with_host_dis$host_disease %in% terms_HP_host_dis,]

# Correct for mutants and exclude certain terms
effective_mut <- c("470.3044","470.9249","83334.615")
HP_left_out_labeling_with_host_dis <- HP_left_out_labeling_with_host_dis[!HP_left_out_labeling_with_host_dis$genome_id %in% effective_mut,]
effective_mutations <- "1313.15800"
HP_left_out_labeling_with_host_dis <- HP_left_out_labeling_with_host_dis[!HP_left_out_labeling_with_host_dis$genome_id %in% "1313.15800",]

# Add to labeled data
Genomes_labeled_v2_edit_host <- bind_rows(Genomes_labeled_edit, HP_left_out_labeling_with_host_dis %>% select(intersect(names(.), names(Genomes_labeled_edit))))
Genomes_labeled_v2_edit_host <- Genomes_labeled_v2_edit_host %>% mutate(Label = ifelse(is.na(Label), "HP", Label))


# Remove specific terms
rm_carriage <- c("1280.21672", "1352.2737", "158836.977", "2587529.4", "486.36", "571.1455", "573.15360", "573.15361", "573.15364", "573.15366", "573.15369", "573.15370", "573.15372", "573.15376")
Genomes_labeled_v2_edit_host <- Genomes_labeled_v2_edit_host[!Genomes_labeled_v2_edit_host$genome_id %in% rm_carriage,]
rm_coloniz <- c("1134687.258","1134687.342", "1259973.7", "1280.31570", "1280.35619", "1282.2497", "1313.34904", "1352.3594", "1352.3595", "1352.3596", "2058152.231", "225992.16", "2608867.60", "299766.116","562.45436","562.45437","562.50775", "562.79034","562.96044","573.33881", "573.35024", "573.7226", "573.7245", "573.7354", "573.7356", "573.7359", "573.9645", "67824.39" )
Genomes_labeled_v2_edit_host <- Genomes_labeled_v2_edit_host[!Genomes_labeled_v2_edit_host$genome_id %in% rm_coloniz,]
Genomes_labeled_v2_edit_host <- Genomes_labeled_v2_edit_host[!Genomes_labeled_v2_edit_host$genome_id %in% "1134687.133",]

# Rename columns for final export
Genomes_labeled_database <- Genomes_labeled_v2_edit_host
names(Genomes_labeled_database)[names(Genomes_labeled_v2_edit_host) == "Label"] <- "pathogenicity_label"
# Remove irrelevant columns
Genomes_labeled_database <- Genomes_labeled_database[, !names(Genomes_labeled_database) %in% c("host_name", "environment (feature)", "fields_pasted")]

# Export final table
write.table(Genomes_labeled_database, "Genomes_labeled_v3_edit_database.txt", sep="\t", row.names = FALSE, quote = FALSE)

# Descriptive statistics
Genomes_labeled_database %>%
  group_by(family) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

Genomes_labeled_database_top10fam <- Genomes_labeled_database %>%
  group_by(family) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  slice_head(n = 10)
