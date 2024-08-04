# Load necessary libraries
library(dplyr)
library(readr)
library(stringr)

# Load updated PATRIC genome metadata
patric_genome_metadata <- read.delim("genome_metadata", dec = ".", colClasses = c("character", rep(NA, 65)))
rownames(patric_genome_metadata) <- patric_genome_metadata$genome_id

# Load BVBRC genome data with specific filters as detailed in publication
BVBRC_genome_nohf_compl <- read.delim("BVBRC_genome_COMPL_NOHF.txt", dec = ".", colClasses = c("character", rep(NA, 86)))
rownames(BVBRC_genome_nohf_compl) <- BVBRC_genome_nohf_compl$Genome.ID

# Create a host table and save it
host_table <- as.data.frame(table(BVBRC_genome_nohf_compl$Host.Name))
host_table <- host_table[order(host_table$Freq, decreasing = TRUE), ]
write.table(host_table, "host_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Filter for human hosts
BVBRC_genome_nohf_compl <- subset(BVBRC_genome_nohf_compl, Host.Name %in% c("Human", "Homo sapiens", "Human, Homo sapiens", 
                                                                                  "Homo sapiens sapiens", "Homo sapiens (adult)", "Homo sapiens (female)", 
                                                                                  "Homo sapiens (respiratory patient)", "Homo sapiens female", "Homo sapiens male",
                                                                                  "Homo sapiens subsp sapiens", "Homo sapiens; 55 year old", "Homo-sapiens",
                                                                                  "Human (adult)", "Human (Female)", "Human (Male)", "Human (Woman)", 
                                                                                  "Human body fluids", "Human female", "Human gut", "Human male", "Human skin", 
                                                                                  "Human woman"))

# Filter based on date inserted
BVBRC_genome_nohf_compl$Date.Inserted <- as.Date(BVBRC_genome_nohf_compl$Date.Inserted)
BVBRC_genome_filt <- BVBRC_genome_nohf_compl[BVBRC_genome_nohf_compl$Date.Inserted > as.Date("2017-01-01"), ]

# Ensure Genome.ID is character type
BVBRC_genome_filt$Genome.ID <- as.character(BVBRC_genome_filt$Genome.ID)

# Match PATRIC genome metadata with filtered BVBRC genome data
patric_genome_metadata_sel <- patric_genome_metadata[intersect(BVBRC_genome_filt$Genome.ID, patric_genome_metadata$genome_id), ]

# Divide data into groups for processing
num_groups <- 4
groups <- cut(seq_along(patric_genome_metadata_sel$biosample_accession), 
              breaks = seq(0, length(patric_genome_metadata_sel$biosample_accession), length.out = num_groups + 1), 
              labels = c("Group 1", "Group 2", "Group 3", "Group 4"), include.lowest = TRUE)

patric_genome_metadata_sel$groups <- groups
group1 <- patric_genome_metadata_sel[patric_genome_metadata_sel$groups=="Group 1",]
group2 <- patric_genome_metadata_sel[patric_genome_metadata_sel$groups=="Group 2",]
group3 <- patric_genome_metadata_sel[patric_genome_metadata_sel$groups=="Group 3",]
group4 <- patric_genome_metadata_sel[patric_genome_metadata_sel$groups=="Group 4",]
write.table(group1$biosample_accession, "bios_acc_patric_group1_v2.txt", quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(group2$biosample_accession, "bios_acc_patric_group2_v2.txt", quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(group3$biosample_accession, "bios_acc_patric_group3_v2.txt", quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(group4$biosample_accession, "bios_acc_patric_group4_v2.txt", quote = FALSE, row.names = FALSE, col.names = FALSE) 
