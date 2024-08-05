library(xml2)
library(dplyr)

# List of XML file paths  
file_paths <- c("biosample_result_PATRIC_group1.xml", 
                "biosample_result_PATRIC_group2.xml",
                "biosample_result_PATRIC_group3.xml",
                "biosample_result_PATRIC_group4.xml")

# Initialize an empty XML document to hold the merged content
merged_xml <- read_xml("<root></root>")

# Loop through each XML file and merge its content into the main document
for (file_path in file_paths) {
  if (!file.exists(file_path)) {
    warning(paste("File does not exist:", file_path))
    next  # Skip this iteration of the loop
  }
  # Proceed with reading the file
  xml_content <- read_xml(file_path)
  
  # Extract the root node of the current XML content
  root_node <- xml_root(xml_content)
  
  # Append the root node to the merged XML document
  xml_add_child(merged_xml, root_node)
}

# Save or print the merged XML document
#write_xml(merged_xml, "merged.xml")
print(merged_xml)


## Function to extract the attributes from xml
extractAttributeData <- function(xml_document) {
  # Parse the XML document
  #doc <- read_xml(xml_document)
  doc <- xml_document
  # Find all BioSample nodes
  biosample_nodes <- xml_find_all(doc, "//BioSample")
  
  # Initialize a list to store dataframes
  list_of_dataframes <- list()
  
  # Process each BioSample node and create a dataframe
  for (i in seq_along(biosample_nodes)) {
    biosample_id <- xml_text(xml_find_first(biosample_nodes[[i]], ".//Id[@db='BioSample' and @is_primary='1']"))
    
    attribute_nodes <- xml_find_all(biosample_nodes[[i]], ".//Attributes/Attribute")
    
    if (length(attribute_nodes) > 0) {
      attribute_data <- data.frame(
        BiosampleID = rep(biosample_id, length(attribute_nodes)),
        AttributeName = xml_attr(attribute_nodes, "attribute_name"),
        AttributeValue = xml_text(attribute_nodes),
        stringsAsFactors = FALSE
      )
      
      list_of_dataframes[[i]] <- attribute_data
    } else {
      list_of_dataframes[[i]] <- data.frame(
        BiosampleID = biosample_id,
        AttributeName = character(),
        AttributeValue = character(),
        stringsAsFactors = FALSE
      )
    }
  }
  
  return(list_of_dataframes)
}

# Apply the function to the XML document and get a list of dataframes - one per xml biosample file
list_of_dataframes <- extractAttributeData(merged_xml)

transposed_dataframes <- lapply(list_of_dataframes, function(df) {
  t_df <- t(df)
  t_df <- as.data.frame(t_df)
  col_names <- as.character(t_df[2, ])
  row_names <- as.character(t_df[1, 1])
  t_df <- t_df[-c(1,2), , drop = FALSE]
  colnames(t_df) <- col_names
  rownames(t_df) <- row_names
  as.data.frame(t_df)  # Convert the matrix back to a data frame
}) #11154

## Combine the columns using dplyr's bind_rows
library(dplyr)
#merged_df <- bind_rows(transposed_dataframes_1, .id = "source")
merged_df <- bind_rows(transposed_dataframes, .id = "source")
#merged_df_xml2 <- bind_rows(transposed_dataframes_xml2, .id = "source")

## Just to assess
env_biome <- merged_df[!is.na(merged_df$"env_biome"),"env_biome"]
env_biome_water <- merged_df[grepl("water", merged_df$env_biome),] # (WATER WASTE)

## Select relevant variables
names(merged_df)
keep_att_patsel <- c(3,4,5,7,11,17,19,24,25,30,36, 38,41,58,62,72,77,
                     85,94,100,109,111,118,125,129,130,143,167,168,180,192,210,217,
                     228,236,280,317,348,371,384,396,410,411,425)

names(merged_df)[keep_att_patsel]
keep_att_patsel <-c("env_broad_scale", "env_local_scale", "env_medium", "isol_growth_condt", "project_name", "health_disease_stat", 
                    "pathogenicity", "sample_type", "isolation_source", "host_disease", "host_health_state", "health_state", 
                    "host_disease_outcome", "host health state", "host_description", "host_disease_stage", "subsrc_note", "pathotype", 
                    "host_disease_stat", "isolation source", "subsource_note", "note", "description", "biotic_relationship", "env_biome", 
                    "env_feature", "host status", "environment (biome)", "environment (feature)", "pathovar", "isolation-source",
                    "host disease", "biome", "Host disease", "host_status", "host_health", "metagenome_source", "risk group", 
                    "observed biotic relationship", "Note on infection mode", "Isolation_source", "Patient Sex", "Patient Age", "disease")

## Subset df with these columns
merged_df_sel <- merged_df[,keep_att_patsel] #11154