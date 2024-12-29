#普通版
library(dplyr)
library(bigreadr)
main.folder <- "/public/home/biostat06/TWAS/FUSION_RESULT"
output_main_folder <- "/public/home/biostat06/TWAS/fusion_result"
subfolders <- list.dirs(main.folder,full.names = TRUE)
for (subfolder in subfolders){
  subfolder_name <- basename(subfolder)
  output_subfolder <- file.path(output_main_folder, subfolder_name)
  if (!dir.exists(output_subfolder)) {
    dir.create(output_subfolder)
  }
  files<-list.files(path = subfolder, pattern="*.dat", full.names = TRUE)
  combined_data <- data.frame()
  for (file in files){
    df <-fread2(file, stringsAsFactors = FALSE, header = T)
    combined_data <- rbind(combined_data, df)
    combined_data<-combined_data %>% 
      filter(TWAS.P > 0.05/(nrow(combined_data)-1))
  }
  output_file_path <- file.path(output_subfolder, paste0(subfolder_name, "_filtered.txt"))
  write.table(combined_data, file = output_file_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

#第一个表有
library(dplyr)
library(bigreadr)

main.folder <- "/public/home/biostat06/TWAS/FUSION_RESULT"
output_main_folder <- "/public/home/biostat06/TWAS/fusion_result"

subfolders <- list.dirs(main.folder, full.names = TRUE)

for (subfolder in subfolders) {
  subfolder_name <- basename(subfolder)
  output_subfolder <- file.path(output_main_folder, subfolder_name)
  if (!dir.exists(output_subfolder)) {
    dir.create(output_subfolder)
  }
  
  files <- list.files(path = subfolder, pattern = "*.dat", full.names = TRUE)
  combined_data <- data.frame()
  header_added <- FALSE
  for (file in files) {
     if (!header_added) {
      df <- fread2(file, stringsAsFactors = FALSE, header = TRUE)
      header_added <- TRUE
      } else {
        df <- fread2(file, skip = 1, stringsAsFactors = FALSE, header = FALSE)
        colnames(df) <- colnames(combined_data) 
            }
       combined_data <- rbind(combined_data, df)
  }
  combined_data <- combined_data %>%
    filter(TWAS.P > 0.05 / (nrow(combined_data) - 1))
  output_file_path <- file.path(output_subfolder, paste0(subfolder_name, "_filtered.txt"))
  write.table(combined_data, file = output_file_path, sep = "\t", quote = FALSE, row.names = FALSE)
}


#手动、
library(data.table)
library(dplyr)
directory <- "/public/home/biostat06/TWAS/FUSION_RESULT/GTExv8.EUR.Thyroid"
dat_files <- list.files(directory, pattern = "\\.dat$", full.names = TRUE) #也可以"*.dat"
dat_list <- lapply(dat_files, fread)
combined_data <- rbindlist(dat_list)
 #print(combined_data)
#output_file <- "/public/home/biostat06/TWAS/fusion_result/Adipose_Subcutaneous.csv"
#fwrite(combined_data, file = output_file, sep = ",")

#df<-fread2("/public/home/biostat06/TWAS/fusion_result/Adipose_Subcutaneous.csv", header = T)
#df_filter<-df %>% filter(TWAS.P > 0.05 / (nrow(combined_data) - 1))
combined_data <- combined_data %>% filter(TWAS.P < 0.05 / (nrow(combined_data) - 1))
output_file <- "/public/home/biostat06/TWAS/fusion_plot/thyroidfilter.csv"
fwrite(combined_data, file = output_file, sep = ",")
