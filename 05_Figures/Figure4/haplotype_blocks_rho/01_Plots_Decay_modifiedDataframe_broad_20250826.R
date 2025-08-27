###########################################################
#                                                         #
# Statistics recombination spots and Cla-element location #
#                   Script by Laura Pettrich              #
#                           REVISED                       #
#                         Juli 2024                       #
###########################################################

## Clean environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(dplyr)
#library(chromoMap)
#library(scales)
#library(ggplot2)

# Set directory
# setwd("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv")
setwd("~/Uni-Köln/PhD/CRIP/plink_haplotype_blocks/run02/")


getwd()

#######################################################################################
# READ IN AND MODIFY DATA
#######################################################################################
#-------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------
# Here you define vectors for every data table you have
# This info will be added a column later
# each population
pop <- c("MF_Chr1", "MF_Chr2", "MF_Chr3", "MF_Chr4", 
         "MG_Chr1", "MG_Chr2", "MG_Chr3", "MG_Chr4", 
         "NMF_Chr1", "NMF_Chr2", "NMF_Chr3", "NMF_Chr4", 
         "SI_Chr1", "SI_Chr2", "SI_Chr3", "SI_Chr4", 
         "SS_Chr1", "SS_Chr2", "SS_Chr3", "SS_Chr4")

# each chromosome
chr <- c("Chr1", "Chr2", "Chr3", "Chr4")

# each individual
ind_numbers <- 1:4
pop_prefix <- c("MF", "MG", "NMF", "SI", "SS")
ind_prefix <- paste(rep(pop_prefix, each = 4), 1:4, sep = "")
# Find indices where "MG" appears and change numbers
mg_indices <- grep("^MG", ind_prefix)
ind_prefix[mg_indices] <- paste("MG", 2:5, sep = "")
print(ind_prefix)

chr_suffixes <- paste("Chr", 1:4, sep = "")

ind <- paste(rep(ind_prefix, each = 4), chr_suffixes, sep = "_")
print(ind)

# Here you tell R the path where it should look for dataframes
# common_path = "C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv/02_less500bp"
common_path = "./"

# Read in bedtools closest results
# These dataframes will be used for analysis
# You add each population seperately
## MF
files_to_read_1 = list.files(
  path = common_path,        # directory to search within
  pattern = "^MF.*\\.bed$", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_1

data_lst_1 <- lapply(files_to_read_1, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})


# Read in bedtools closest results
## MG
files_to_read_2 = list.files(
  path = common_path,        # directory to search within
  pattern = "^MG.*\\.bed$", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_2

data_lst_2 <- lapply(files_to_read_2, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})

# Read in bedtools closest results
## NMF
files_to_read_3 = list.files(
  path = common_path,        # directory to search within
  pattern = "^NMF.*\\.bed$", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_3

data_lst_3 <- lapply(files_to_read_3, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})

# Read in bedtools closest results
## SI
files_to_read_4 = list.files(
  path = common_path,        # directory to search within
  pattern = "^SI.*\\.bed$", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_4

data_lst_4 <- lapply(files_to_read_4, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})

# Read in bedtools closest results
## SS
files_to_read_5 = list.files(
  path = common_path,        # directory to search within
  pattern = "^SS.*\\.bed$", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_5

data_lst_5 <- lapply(files_to_read_5, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})

#-------------------------------------------------------------------
# Rename columns
#-------------------------------------------------------------------
# Rename data list names based on files_to_read
names(data_lst_1) <- ind[1:16]
names(data_lst_2) <- ind[17:32]
names(data_lst_3) <- ind[33:48]
names(data_lst_4) <- ind[49:64]
names(data_lst_5) <- ind[65:80]

# Check if renaming was done correctly
print(names(data_lst_1))
print(names(data_lst_2))
print(names(data_lst_3))
print(names(data_lst_4))
print(names(data_lst_5))

# Drop dataframes with NULL enttries
for (i in 1:5) {
  obj_name <- paste0("data_lst_", i)
  assign(obj_name, Filter(function(df) !is.null(df), get(obj_name)))
}


# Create a list of all the data frames
data_lst <- list(data_lst_1, data_lst_2, data_lst_3, data_lst_4, data_lst_5)

# Add a "ind" column to each data frame with the corresponding 'ind' value
for (j in 1:5) {
  for (i in 1:length(data_lst[[j]])) {
    data_lst[[j]][[i]]$ind <- names(data_lst[[j]])[i]
  }
}

# Create separate data lists (data_lst_1, data_lst_2, etc.) from data_lst
for (j in 1:5) {
  assign(paste0("data_lst_", j), data_lst[[j]])
}

# Remove all dataframes with no Cla-Element in adjacent regions (value -1 in column v6 - v8)
# Alternatively, if you want to remove the entire dataframe from the list
data_lst_1 <- data_lst_1[sapply(data_lst_1, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]
data_lst_2 <- data_lst_2[sapply(data_lst_2, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]
data_lst_3 <- data_lst_3[sapply(data_lst_3, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]
data_lst_4 <- data_lst_4[sapply(data_lst_4, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]
data_lst_5 <- data_lst_5[sapply(data_lst_5, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]



# Mutate dataframe 
df1 <- bind_rows(data_lst_1)
df2 <- bind_rows(data_lst_2)
df3 <- bind_rows(data_lst_3)
df4 <- bind_rows(data_lst_4)
df5 <- bind_rows(data_lst_5)

# Add column with Pop name
df1$pop <- "MF"
df2$pop <- "MG"
df3$pop <- "NMF"
df4$pop <- "SI"
df5$pop <- "SS"

# Bind data frames together
dfm <- rbind(df1, df2, df3, df4, df5)
head(dfm)
tail(dfm)
# columns means:
# "chromosome", "start", "end", "rho", 
# "chromosome", "start", "end", "insertions", "score", "strand"
# "distance", "ind", "pop"

# Rename columns in data frame
new_colnames <- c("chromosome", "start", "end", "rho", "chromosome_ins", "start_hb", "end_hb",
                  "dist_to_hb", "ind_chr", "pop") # look up what the score means

dfm <- dfm %>%
  rename_with(~ new_colnames)
head(dfm)

write.table(dfm, "~/sciebo/RecombinationLandscape_CRIP/02_Revision_2/BMC_Genomics_CripRecLand/Revision_20250820/NewResults/bedtools_closest_ouput_1kb_rho_hb_broad.csv", sep = ",", col.names = TRUE)

head(dfm)

#-----------------------------------
# 100 kb broad


###########################################################
#                                                         #
# Statistics recombination spots and Cla-element location #
#                   Script by Laura Pettrich              #
#                           REVISED                       #
#                         Juli 2024                       #
###########################################################

## Clean environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(dplyr)
#library(chromoMap)
#library(scales)
#library(ggplot2)

# Set directory
# setwd("C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv")
setwd("~/Uni-Köln/PhD/CRIP/plink_haplotype_blocks/run03/")


getwd()

#######################################################################################
# READ IN AND MODIFY DATA
#######################################################################################
#-------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------
# Here you define vectors for every data table you have
# This info will be added a column later
# each population
pop <- c("MF_Chr1", "MF_Chr2", "MF_Chr3", "MF_Chr4", 
         "MG_Chr1", "MG_Chr2", "MG_Chr3", "MG_Chr4", 
         "NMF_Chr1", "NMF_Chr2", "NMF_Chr3", "NMF_Chr4", 
         "SI_Chr1", "SI_Chr2", "SI_Chr3", "SI_Chr4", 
         "SS_Chr1", "SS_Chr2", "SS_Chr3", "SS_Chr4")

# each chromosome
chr <- c("Chr1", "Chr2", "Chr3", "Chr4")

# each individual
ind_numbers <- 1:4
pop_prefix <- c("MF", "MG", "NMF", "SI", "SS")
ind_prefix <- paste(rep(pop_prefix, each = 4), 1:4, sep = "")
# Find indices where "MG" appears and change numbers
mg_indices <- grep("^MG", ind_prefix)
ind_prefix[mg_indices] <- paste("MG", 2:5, sep = "")
print(ind_prefix)

chr_suffixes <- paste("Chr", 1:4, sep = "")

ind <- paste(rep(ind_prefix, each = 4), chr_suffixes, sep = "_")
print(ind)

# Here you tell R the path where it should look for dataframes
# common_path = "C:/Users/laupe/Documents/sciebo/RecombinationLandscape_CRIP/Scripts/01_Cla-annotation_IndivAnalysis/run07_Indiv/02_less500bp"
common_path = "./"

# Read in bedtools closest results
# These dataframes will be used for analysis
# You add each population seperately
## MF
files_to_read_1 = list.files(
  path = common_path,        # directory to search within
  pattern = "^MF.*\\.bed$", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_1

data_lst_1 <- lapply(files_to_read_1, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})


# Read in bedtools closest results
## MG
files_to_read_2 = list.files(
  path = common_path,        # directory to search within
  pattern = "^MG.*\\.bed$", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_2

data_lst_2 <- lapply(files_to_read_2, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})

# Read in bedtools closest results
## NMF
files_to_read_3 = list.files(
  path = common_path,        # directory to search within
  pattern = "^NMF.*\\.bed$", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_3

data_lst_3 <- lapply(files_to_read_3, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})

# Read in bedtools closest results
## SI
files_to_read_4 = list.files(
  path = common_path,        # directory to search within
  pattern = "^SI.*\\.bed$", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_4

data_lst_4 <- lapply(files_to_read_4, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})

# Read in bedtools closest results
## SS
files_to_read_5 = list.files(
  path = common_path,        # directory to search within
  pattern = "^SS.*\\.bed$", # regex pattern, some explanation below
  recursive = FALSE,          # search subdirectories
  full.names = TRUE          # return the full path
)
files_to_read_5

data_lst_5 <- lapply(files_to_read_5, function(file) {
  if (file.size(file) > 0) {
    read.table(file, header = FALSE)
  } else {
    # Handle the case where the file is empty
    # You can choose to return NULL or some other value
    NULL
  }
})

#-------------------------------------------------------------------
# Rename columns
#-------------------------------------------------------------------
# Rename data list names based on files_to_read
names(data_lst_1) <- ind[1:16]
names(data_lst_2) <- ind[17:32]
names(data_lst_3) <- ind[33:48]
names(data_lst_4) <- ind[49:64]
names(data_lst_5) <- ind[65:80]

# Check if renaming was done correctly
print(names(data_lst_1))
print(names(data_lst_2))
print(names(data_lst_3))
print(names(data_lst_4))
print(names(data_lst_5))

# Drop dataframes with NULL enttries
for (i in 1:5) {
  obj_name <- paste0("data_lst_", i)
  assign(obj_name, Filter(function(df) !is.null(df), get(obj_name)))
}


# Create a list of all the data frames
data_lst <- list(data_lst_1, data_lst_2, data_lst_3, data_lst_4, data_lst_5)

# Add a "ind" column to each data frame with the corresponding 'ind' value
for (j in 1:5) {
  for (i in 1:length(data_lst[[j]])) {
    data_lst[[j]][[i]]$ind <- names(data_lst[[j]])[i]
  }
}

# Create separate data lists (data_lst_1, data_lst_2, etc.) from data_lst
for (j in 1:5) {
  assign(paste0("data_lst_", j), data_lst[[j]])
}

# Remove all dataframes with no Cla-Element in adjacent regions (value -1 in column v6 - v8)
# Alternatively, if you want to remove the entire dataframe from the list
data_lst_1 <- data_lst_1[sapply(data_lst_1, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]
data_lst_2 <- data_lst_2[sapply(data_lst_2, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]
data_lst_3 <- data_lst_3[sapply(data_lst_3, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]
data_lst_4 <- data_lst_4[sapply(data_lst_4, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]
data_lst_5 <- data_lst_5[sapply(data_lst_5, 
                                function(df) all(!(df$V6 == -1 & df$V7 == -1 & df$V8 == -1)))]



# Mutate dataframe 
df1 <- bind_rows(data_lst_1)
df2 <- bind_rows(data_lst_2)
df3 <- bind_rows(data_lst_3)
df4 <- bind_rows(data_lst_4)
df5 <- bind_rows(data_lst_5)

# Add column with Pop name
df1$pop <- "MF"
df2$pop <- "MG"
df3$pop <- "NMF"
df4$pop <- "SI"
df5$pop <- "SS"

# Bind data frames together
dfm <- rbind(df1, df2, df3, df4, df5)
head(dfm)
tail(dfm)
# columns means:
# "chromosome", "start", "end", "rho", 
# "chromosome", "start", "end", "insertions", "score", "strand"
# "distance", "ind", "pop"

# Rename columns in data frame
new_colnames <- c("chromosome", "start", "end", "rho", "chromosome_ins", "start_hb", "end_hb",
                  "dist_to_hb", "ind_chr", "pop") # look up what the score means

dfm <- dfm %>%
  rename_with(~ new_colnames)
head(dfm)

write.table(dfm, "~/sciebo/RecombinationLandscape_CRIP/02_Revision_2/BMC_Genomics_CripRecLand/Revision_20250820/NewResults/bedtools_closest_ouput_100kb_rho_hb_broad.csv", sep = ",", col.names = TRUE)

head(dfm)