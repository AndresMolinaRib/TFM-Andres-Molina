library(dplyr)
library(msdap)
library(limma)
library(Biostrings)
library(writexl)
library(openxlsx)
library(readxl)

#Import files
dataset <- import_dataset_diann(filename = "report.tsv")
dataset <- import_fasta(dataset, files = "archivo_combinado.fasta")
dataset <- import_sample_metadata(dataset, filename = "sample_metadata.xlsx")


#Create contrast
dataset = setup_contrasts(dataset,
                          contrast_list = list(c("MWd4_L","MSd4_L"),
                                               c("MWd4_Exo","MSd4_Exo"),
                                               c("MWd4_S","MSd4_S"),
                                               c("FW_L","FS_L"),
                                               c("FW_S","FS_S"),
                                               c("MW_L","MS_L"),
                                               c("MW_S","MS_S")
                          ))


#Run analysis
dataset = analysis_quickstart(dataset,
                              filter_min_detect = 1,
                              filter_min_quant = 1,
                              filter_fraction_detect = 0.0,
                              filter_fraction_quant = 0.0,
                              filter_min_peptide_per_prot = 1,
                              filter_by_contrast = TRUE,
                              norm_algorithm = c(
                                "vsn",
                                "modebetween_protein"
                              ),
                              dea_algorithm = c("deqms",
                                                "msempire", "msqrob"
                              ),
                              dea_qvalue_threshold = 0.05,
                              dea_log2foldchange_threshold = NA,
                              diffdetect_min_peptides_observed = 1,
                              diffdetect_min_samples_observed = 1,
                              diffdetect_min_fraction_observed = 0,
                              output_qc_report = TRUE,
                              output_abundance_tables = TRUE,
                              output_dir = "msdap_results",
                              output_within_timestamped_subdirectory = TRUE
)

## Create dataset for pipeline
protein_data = get_protein_matrix(dataset, intensity_column = "intensity_by_group")


protein_data_matrix <- as.data.frame(protein_data$matrix)

df_with_rownames <- rownames_to_column(protein_data_matrix, var = "protein_id")

result = inner_join(df_with_rownames, data_complete, by = "protein_id")


# Extract standard deviation
sd <- as.data.frame(dataset$de_proteins$standarddeviation)
sd$protein_id <- dataset$de_proteins$protein_id
sd$contrast <- dataset$de_proteins$contrast
sd$method <- dataset$de_proteins$dea_algorithm

df <- sd 

colnames(df)[1] <- "standarderror"

df$combinada <- paste(df$method, df$contrast, sep = "_")


df <- df %>%
  mutate(standarderror = as.numeric(standarderror))


# Delete duplications
df_sin_duplicados <- df %>%
  distinct(protein_id, combinada, .keep_all = TRUE)


df_wide <- df_sin_duplicados %>%
  pivot_wider(
    id_cols = protein_id,
    names_from = combinada,
    values_from = standarderror
  )

df_wide <- df_wide %>%
  rename_with(
    .fn = ~ paste0("sd_", .),
    .cols = starts_with("msempire")
  )

df_wide <- df_wide %>%
  rename_with(
    .fn = ~ paste0("sd_", .),
    .cols = starts_with("msqrob")
  )

df_wide <- df_wide %>%
  rename_with(
    .fn = ~ paste0("sd_", .),
    .cols = starts_with("deqms")
  )

resultado = inner_join(result, df_wide, by = "protein_id")

write_xlsx(resultado,"dea_final_prueba.xlsx" )
