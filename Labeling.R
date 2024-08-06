# Load necessary libraries
library(dplyr)
library(readr)
library(stringr)

# ------------------------------------------------------------------------------
# Create a combined field for inspection 
# all_metadata_filt as the integrated and preprocessed dataset 
# ------------------------------------------------------------------------------

# Combine multiple fields into one for labeling purposes
all_metadata_filt$fields_pasted <- paste(all_metadata_filt$isolation_source, all_metadata_filt$host_health, all_metadata_filt$disease,
                                         all_metadata_filt$comments, all_metadata_filt$isolation_comments, 
                                         all_metadata_filt$other_clinical, all_metadata_filt$additional_metadata, 
                                         all_metadata_filt$env_broad_scale, all_metadata_filt$env_local_scale, all_metadata_filt$env_medium,
                                         all_metadata_filt$isol_growth_condt, all_metadata_filt$project_name,
                                         all_metadata_filt$host_health_state, all_metadata_filt$health_state,all_metadata_filt$host_disease_outcome,
                                         all_metadata_filt$`host health state`, all_metadata_filt$host_description, all_metadata_filt$host_disease_stage,
                                         all_metadata_filt$subsrc_note, all_metadata_filt$subsource_note, all_metadata_filt$pathotype, all_metadata_filt$pathogenicity, all_metadata_filt$host_disease_stat,
                                         all_metadata_filt$`isolation source`, all_metadata_filt$isolation_source_ncbi, all_metadata_filt$`isolation-source`, all_metadata_filt$Isolation_source,
                                         all_metadata_filt$note, all_metadata_filt$description, all_metadata_filt$biotic_relationship, all_metadata_filt$env_biome,
                                         all_metadata_filt$env_feature, all_metadata_filt$`host status`, all_metadata_filt$`environment (biome)`, all_metadata_filt$`environment (feature)`,
                                         all_metadata_filt$pathovar_ncbi, all_metadata_filt$`host disease`, all_metadata_filt$biome,
                                         all_metadata_filt$`Host disease`, all_metadata_filt$host_status, all_metadata_filt$host_health_ncbi, all_metadata_filt$metagenome_source,
                                         all_metadata_filt$`observed biotic relationship`, all_metadata_filt$`Note on infection mode`, all_metadata_filt$disease_ncbi, all_metadata_filt$health_state,
                                         all_metadata_filt$`risk group`, all_metadata_filt$health_disease_stat)

# ------------------------------------------------------------------------------
# Filter out entries related to wastewater origins
# ------------------------------------------------------------------------------

# Remove entries isolated from wastewater - verified manually
all_metadata_filt <- all_metadata_filt[!grepl("water", all_metadata_filt$env_biome, ignore.case = TRUE), ]


# ------------------------------------------------------------------------------
# Labeling for Pathogenic to humans (HP) and Non-pathogenic to humans (NHP)
# ------------------------------------------------------------------------------

## 1) After assessing iteratively, refine HP keywords
HP_keywords <- c("virulence", "superbug", "waterborne", "foodborne", "outbreak", "infection", "pathogen", "water borne",
              "food borne", "itis\\b", "poisoning", "infectious", "sepsis", "infected", "biofilm", "purulent", "pus", "death",
              "severe", "diseased", "pandemic", "epidemic","Transmission","vector", "toxin", "toxic", "clinical",
              "Biosafety Level 2", "hypervirulent", "diarrhea", "intensive", "disease", "STEC")

HP_pattern <- paste(HP_keywords, collapse = "|")
HP_genomes_a2 <- all_metadata_filt[grepl(HP_pattern, all_metadata_filt$fields_pasted, ignore.case = TRUE), ]

# Remove inconclusive entries and projects (examples, after manual assessment of entries)
HP_disease_term_except <- c("1177574.116","108619.365", "108619.366", "115981.207", "1309.145", "1313.15800", "1313.17716", 
                            "1352.1643", "1490.105", "149390.28", "149539.2309", "149539.2774", "149539.2775", "172045.55",
                            "1812935.462", "192954.128", "199.901", "192954.129", "195.3310", "195.3311", "195.3312", "195.3313",
                            "197.20875", "197.20876", "197.20877", "197.20893", "197.20894", "199.901", "210.10653", "210.10654",
                            "210.10655", "2109685.3", "2109687.3", "2109688.3", "2109690.3", "2109691.3", "2109692.3", "280145.39",
                            "287.29820", "287.30241", "287.30278", "28901.9943", "38323.49", "470.13785", "487.123", "54291.284", "556269.12",
                            "562.104348", "562.104349", "562.104350", "562.104351", "562.104352", "562.104353", "562.39410",
                            "59201.3757", "617123.56", "90371.3192")
limbo_projects <- c("PRJNA478278", "PRJNA316728", "PRJNA393749", "PRJNA358390", "PRJNA689979")
HP_genomes_a2 <- HP_genomes_a2[!HP_genomes_a2$genome_id %in% HP_disease_term_except,] 

## 2) In addition, the same fields cannot include any of the EXCLUSION keywords (Non-Pathological OR INCONCLUSIVE keywords).
## same rationale than 1) to refine the keywords
NHP_keywords<- c("Healthy","probiotic", "Commensal", "microbiome", "microbiota", "nutraceutical", "\\bnormal\\b",
              "asymptomatic", "naturally occurring", "human-associated habitat", "opportunistic")
NHP_keywords<- paste(NHP_keywords, collapse="|")
NHP_entries11 <- HP_genomes_a2[grepl(NHP_keywords, HP_genomes_a2$fields_pasted, ignore.case = TRUE),]
# EXCEPTIONS - ADD TO HP
except_a2_2_HP <- c("28188.16", "1314.966", "1902136.3", "1280.34801", "1280.37108", "1280.37109","487.3961", "1280.19052")
# EXCEPTIONS - ADD TO NHP
except_a2_2_NHP <- c("573.15377", "1679.51", "28037.1359", "28037.1407", "239935.2451", "1682.170", "39488.252","1520.583","1492.235")

## FINAL A2_2 - remove ones with NHPkeywords, except for ones which are HP (manually verified)
HP_genomes_a2_2 <- HP_genomes_a2[!(grepl(NHP_keywords, HP_genomes_a2$fields_pasted, ignore.case = TRUE)),]
# add ones that are still HP
HP_genomes_a2_2 <- rbind(HP_genomes_a2_2, HP_genomes_a2[HP_genomes_a2$genome_id %in% except_a2_2_HP,]) 


## HP - B2 PATHOVAR + PATHOGENICITY COLUMN WITH CONDITIONS
# b2_1) Pathovar as not "" except for "not applicable" (same as previous HP_c) from prev analysis)
# remove 
HP_genomes_b2_1_ncbi <- all_metadata_filt[!(all_metadata_filt$pathovar_ncbi == ""| is.na(all_metadata_filt$pathovar_ncbi)),] 
all_metadata_filt <- all_metadata_filt[!all_metadata_filt$genome_id %in% HP_genomes_b2_1_ncbi$genome_id,]
HP_genomes_b2_1 <- all_metadata_filt[!(all_metadata_filt$pathovar == "" | all_metadata_filt$pathovar == "not applicable" | is.na(all_metadata_filt$pathovar)), ] 
# b2_2) "pathotype" values other than "missing", "not applicable", "NA", "not available" 
excl_keywords_pathotype <- c("missing","not applicable","NA","not available")
HP_genomes_b2_2 <- all_metadata_filt[!is.na(all_metadata_filt$pathotype) & !(all_metadata_filt$pathotype %in% excl_keywords_pathotype), ]
# Union of a and b
HP_genomes_b2 <- rbind(HP_genomes_b2_1, HP_genomes_b2_2) 
nrow(HP_genomes_b2)
# Union of HP_genomes_a2_2 with HP_genomes_b2
HP_genomes <- rbind(HP_genomes_a2_2, HP_genomes_b2) 
# remove duplicated entries
HP_genomes <- HP_genomes[!(duplicated(HP_genomes$genome_id)),] 

# Add bioterrorist+common pathogens if not present (https://argos.igs.umaryland.edu/doc/pdf/wanted-orgnaism-list-Jan2019.pdf)
# remove H37rA (laboratory) "419947.17", 1313.35815 (lab strain) - remove from HP and then do rbind to update
HP_genomes <- HP_genomes[!(HP_genomes$genome_id == "419947.17"), ] 
# ADD : 632.861,  1392.1053; 1392.1031; 1392.1047; 1392.1050; 1392.1064; 1392.1278, 446.641 
HP_genomes_to_add <- c("632.861",  "1392.1053", "1392.1031", "1392.1047", "1392.1050", "1392.1064", "1392.1278", "446.641" )
HP_genomes <- rbind(HP_genomes, all_metadata_filt[which(all_metadata_filt$genome_id %in% HP_genomes_to_add), ])

## EXCLUDE MUTANTS
# search for "mutant", "mutations" 
HP_genomes_mutants <- HP_genomes[grepl("mutant", HP_genomes$fields_pasted, ignore.case = TRUE),c("genome_id","fields_pasted")] # 23
# Verify manually, exclude
HP_genomes_nomutants <- HP_genomes[!HP_genomes$genome_id %in% HP_genomes_mutants$genome_id,] 
# search for "mutations" 
HP_genomes_mutations <- HP_genomes_nomutants[grepl("mutation", HP_genomes_nomutants$fields_pasted, ignore.case = TRUE),c("genome_id","fields_pasted")] 
# search for words starting with "edit"
HP_genomes_edit <- HP_genomes_nomutants[grepl("\\bedit\\w*\\b", HP_genomes_nomutants$fields_pasted, ignore.case = TRUE),c("genome_id","fields_pasted")] #1
HP_genomes_no_edit <- HP_genomes_nomutants[!HP_genomes_nomutants$genome_id %in% HP_genomes_edit$genome_id,] # (-1)
HP_genomes <- HP_genomes_no_edit
# Add to HP
addto_HP <- c("1496.5982", "1496.5983", "1496.5981")
HP_genomes <- rbind(HP_genomes, all_metadata_filt[which(all_metadata_filt$genome_id %in% addto_HP), ])
HP_genomes <- HP_genomes[!(duplicated(HP_genomes$genome_id)),] 


## NHP SELECTION ##
## 1) Exclude Generated HP list:
NHP_genomes_1 <- all_metadata_filt[!(all_metadata_filt$genome_id %in% HP_genomes$genome_id), ] 
# And  ones that were in "limbo":
NHP_genomes_1 <- NHP_genomes_1[!(NHP_genomes_1$genome_id %in% limbo), ] 
NHP_genomes_1 <- NHP_genomes_1[!NHP_genomes_1$bioproject_accession %in% limbo_projects,]
# MANUAL -  Add to NHP one with "commensal" in "pathogenicity"
patho_commensal <- all_metadata_filt[!is.na(all_metadata_filt$pathogenicity) & all_metadata_filt$pathogenicity == "commensal", ]
NHP_genomes_1 <- rbind(NHP_genomes_1, all_metadata_filt[which(all_metadata_filt$genome_id %in% patho_commensal$genome_id), ]) # 

## 2) Genomes are labeled as NHP if they contain NHP keywords
# 1 - Remove ones in limbo (manually assessed)
# remove 199.901, 562.54032, 562.54035, 562.54031 and others
healthy_torem <- c("199.901", "562.54032", "562.54035", "562.54031", 	"573.34237", "573.34302", "573.34252")
NHP_genomes_1 <- NHP_genomes_1[!NHP_genomes_1$genome_id %in% healthy_torem,]
NHP_genomes_1 <- NHP_genomes_1[!NHP_genomes_1$genome_id=="1280.21781",]
NHP_genomes_1 <- NHP_genomes_1[!NHP_genomes_1$genome_id=="1636603.3",] 
nrow(NHP_genomes_1) 
NHP_keywords_NHP <- c("Healthy","probiotic", "Commensal", "microbiome", "microbiota", "symbiotic", "nutraceutical", "\\bnormal\\b", "commercial", "flora")
NHP_keywords_NHP <- paste(NHP_keywords_NHP, collapse="|")
NHP_entries_NHP11 <- NHP_genomes_1[grepl(NHP_keywords_NHP, NHP_genomes_1$fields_pasted, ignore.case = TRUE),] 
dif_genms_NHP11 <- setdiff(NHP_entries_NHP11$genome_id, NHP_entries_NHP10$genome_id) 
manual_assess <- NHP_entries_NHP11[NHP_entries_NHP11$genome_id %in% dif_genms_NHP11,c("genome_id", "fields_pasted", "bioproject_accession")]

## Final
NHP_genomes_2 <- NHP_genomes_1[grepl(NHP_keywords_NHP, NHP_genomes_1$fields_pasted, ignore.case = TRUE),] 
nrow(NHP_genomes_2) 

### 3 ) In addition, the same fields cannot include any of the HPkeywords.
# 1 - Remove ones in limbo (manually assessed) and some add to HP 
# remove
del <- c("1309.125", "573.33903", "699240.12", "485.18123", "2067421.3", "333849.47", "562.50493", "40215.119", "714.686")
NHP_genomes_2 <- NHP_genomes_2[!NHP_genomes_2$genome_id %in% del, ]
addto_HP <- c("1496.5982", "1496.5983", "1496.5981", "485.18123", "562.50493", "40215.119", "714.686" , "1313.35611", "1280.34802", "1280.34803", "1280.34804", "1338.31")
# add
HP_genomes <- rbind(HP_genomes, all_metadata_filt[which(all_metadata_filt$genome_id %in% addto_HP), ])
HP_genomes <- HP_genomes[!(duplicated(HP_genomes$genome_id)),] 

NHP_keywords_nonNHP14 <- c("patient","abscess", "wound", "bacteremia", "pneumonia", "ICU ", "disease", "contaminated", "symptom", "clinic", "opportunistic", "distress") # pathogenic K. pneumoniae
NHP_keywords_nonNHP14 <- paste(NHP_keywords_nonNHP14, collapse="|")
NHP_entries_nonNHP14 <- NHP_genomes_2[grepl(NHP_keywords_nonNHP14, NHP_genomes_2$fields_pasted, ignore.case = TRUE),]
dif_genms_NHP14 <- setdiff(NHP_entries_nonNHP14$genome_id, NHP_entries_nonNHP13$genome_id) 
manual_assess <- NHP_entries_nonNHP14[NHP_entries_nonNHP14$genome_id %in% dif_genms_NHP14,c("genome_id", "fields_pasted", "bioproject_accession")]

# Final
NHP_genomes_2 <- NHP_genomes_2[!grepl(NHP_keywords_nonNHP14, NHP_genomes_2$fields_pasted, ignore.case = TRUE),] 
# with exceptions: 
keep_NHP <- c("1290.286", "1292.160", "1304.1560", "1747.298", "1681.94", "1805478.3", "2109685.3", "2109687.3", "1841863.3", 
              "1805478.3","1870984.4", "1944646.3", "2109685.3", "2109687.3", "2109688.3", "2109690.3", "2109691.3", "2109692.3", "2545800.3","40216.95",
              "43675.878", "487.1231", "480.268") 

NHP_genomes_2 <- rbind(NHP_genomes_2, all_metadata_filt[which(all_metadata_filt$genome_id %in% keep_NHP), ])
nrow(NHP_genomes_2)

## ALSO if (verification phase):
# "host_disease_outcome" is NOT other than "not applicable", "none", "Missing", "Unknown", "colonised", "Expired" OR ""
table(NHP_genomes_2$host_disease_outcome)
table(NHP_genomes_2$host_health)
table(NHP_genomes_2$host_disease)
table(NHP_genomes_2$`host health state`)
table(NHP_genomes_2$host_disease_stage)
table(NHP_genomes_2$`host disease`)
table(NHP_genomes_2$host_health_ncbi)
table(NHP_genomes_2$`Host disease`)
table(NHP_genomes_2$disease_ncbi)
table(NHP_genomes_2$pathovar)
table(NHP_genomes_2$"observed biotic relationship")
table(NHP_genomes_2$"host status")
# ALL NHP
# if "Patient Sex" or "Patient Age" != "" then HP
table(NHP_genomes_2$"Patient Sex") # 0
table(NHP_genomes_2$"Patient Age") # 0
NHP_genomes <- NHP_genomes_2 # 522
## Mutants verification 
NHP_genomes_mutants <- NHP_genomes[grepl("mutant", NHP_genomes$fields_pasted, ignore.case = TRUE),c("genome_id","fields_pasted")] # 23
NHP_genomes_mutations <- NHP_genomes[grepl("mutations", NHP_genomes$fields_pasted, ignore.case = TRUE),c("genome_id","fields_pasted")] # 124
NHP_genomes_edit <- NHP_genomes[grepl("\\bedit\\w*\\b", NHP_genomes$fields_pasted, ignore.case = TRUE),c("genome_id","fields_pasted")] #1


## FINAL
HP_genomes$Label <- "HP" # 11330
NHP_genomes$Label <- "NHP" # 
## MERGE HP and NHP
Genomes_labeled <- rbind(HP_genomes, NHP_genomes)
HP_genomes <- HP_genomes[!(duplicated(HP_genomes$genome_id)),] 
NHP_genomes <- NHP_genomes[!(duplicated(NHP_genomes$genome_id)),] 

## not opportunistic in fields 
Genomes_labeled <- Genomes_labeled[!grepl("opportunistic", Genomes_labeled$fields_pasted, ignore.case = TRUE),] # 374; 309 (v2)
exclude_keywords <- "diabetes|cancer|hiv|aids|leukemia|carcinoma|guillain|AIDS|immunocompromised|interleukin-12 receptor deficiency"
Genomes_labeled <- Genomes_labeled[grepl(exclude_keywords, Genomes_labeled$fields_pasted, ignore.case = TRUE), ]

## Reassessments according to manual review guided by other terms or insights

inconclusive <- c("1091046.5","1280.17963","1280.17965" ,"1280.16544", "1280.16545", "1313.15800", "1352.3597", "1352.3598", "135487.47", "1390.556", "1396.2548", "1409.43", "158836.446", "1747.1400", "1747.868", 
                  "263.153", "263.154", "263.155", "2702.191", "28035.68", "28035.69", "28035.70", "28080.85", "287.7095", "303.689", "545.130", "546.1525", "546.1526", "550.1577", 
                  "550.3158", "562.79032", "573.15371", "573.15374", "573.15375", "573.41302", "582.117", "61647.65", "777.114", "817.2528", "83554.72", "83554.73", "83554.74", "991915.3", "562.80701"
)
genid_keep_FDA <- c("727.1050", "28450.828", "487.1303", "485.881", "1305.692", "676.375", "1338.30", "28035.29", "28035.30", "747.530", "106654.118", "28035.20", "50719.19", "34105.36", "1839798.3", "1302.125", "1303.167", "57706.19", "1117645.198", "490.6", "1311.1349", "1335.40", "1302.83", "1304.240", "666.3112", "485.879", "545.38", "13373.94", "87883.148", "210.4220", "32022.293", "2018067.3", "210.4150", "2545800.3", "44275.5", "673.18", "670.1110", "50719.21", "1660.163", "2545799.3", "823.3406", "47770.612", "37326.9", "1290.286", "1292.160", "1304.1560", "40216.95", "1717.335", "57975.30", "480.268", "487.1302", "487.1300", "487.1301", "1773.22837", "562.20517", "1038927.31", "1038927.40", "556499.3", "1353.23", "450.24", "827.94", "571.319", "29385.174", "1463165.70", "562.28386", "573.17909", "587.58", "29385.175", "1260.177")
Genomes_labeled <- Genomes_labeled[!Genomes_labeled$bioproject_accession %in% c("PRJNA689979", "PRJNA428178", "PRJNA393749", "PRJNA358390"),]
Genomes_labeled <- Genomes_labeled[!(Genomes_labeled$bioproject_accession == "PRJNA231221") | (Genomes_labeled$bioproject_accession == "PRJNA231221" & Genomes_labeled$genome_id %in% genid_keep_FDA), ]
Genomes_labeled <- Genomes_labeled[!Genomes_labeled$genome_id %in% inconclusive,  ] # -44, 6749
rm_carriage <- c("1280.21672", "1352.2737", "158836.977", "2587529.4", "486.36", "571.1455", "573.15360", "573.15361", "573.15364", "573.15366", "573.15369", "573.15370", "573.15372", "573.15376")
Genomes_labeled <- Genomes_labeled[!Genomes_labeled$genome_id %in% rm_carriage,]
rm_coloniz <- c("1134687.258","1134687.342", "1259973.7", "1280.31570", "1280.35619", "1282.2497", "1313.34904", "1352.3594", "1352.3595", "1352.3596", "2058152.231", "225992.16", "2608867.60", "299766.116","562.45436","562.45437","562.50775", "562.79034","562.96044","573.33881", "573.35024", "573.7226", "573.7245", "573.7354", "573.7356", "573.7359", "573.9645", "67824.39" )
Genomes_labeled <- Genomes_labeled[!Genomes_labeled$genome_id %in% rm_coloniz,]
Genomes_labeled <- Genomes_labeled[!Genomes_labeled$genome_id %in% "1134687.133",]


# ------------------------------------------------------------------------------
# Export final labeled dataset - columns to be further assessed in "Final_dbase_genomes_labeled"
# ------------------------------------------------------------------------------
write.csv(Genomes_labeled, file = "Genomes_labeled.csv", row.names = FALSE)

