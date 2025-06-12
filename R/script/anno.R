# Functional Annotation

# ---- Load Packages ----
library("tidyverse")
library("tidylog")
library("ggpubr")

setwd("./R")


# ---- Load Tables ----
# Reference
cat <- read.delim("functionalanno/ref/cat.hierarchy24.txt", header = FALSE, sep = "\t")
cog <- read.delim("functionalanno/ref/cog.to.cat24.txt", header = FALSE, sep = "\t")

# fix column names
colnames(cat) <- c("Letter", "Function", "Category")
colnames(cog) <- c("Letter", "COG")

# Read in files 
# Directory containing the .cog files
cog_dir <- "functionalanno/cogfiles"

# List all files ending in _4.cog
cog_files <- list.files(cog_dir, pattern = "_3\\.cog$", full.names = TRUE)

# Placeholder for combined results
cogannos <- list()

# Loop over each file
for (file in cog_files) {
  
  file_name <- basename(file)
  group <- str_remove(file_name, "_3\\.cog")
  
  
  tab <- read.delim(file, sep = "\t", header = FALSE)
  
  reformatted <- tab %>%
    mutate(
      sra_accession = str_remove_all(V1, "./"),
      sra_accession = str_remove_all(sra_accession, ".cog"),
      gene = str_remove_all(V3, "gene="),
      COG = str_remove_all(V2, "COG:"),
      group = group  # Add sample_id as a new column
    ) %>%
    select(group, sra_accession, gene, COG) %>%
    full_join(cog, by = "COG") %>%
    full_join(cat, by = "Letter") %>%
    group_by(gene) %>%
    mutate(Nstrains = n()) %>%
    ungroup() %>%
    select(group, gene, Nstrains, COG, Letter, Function, Category) %>%
    unique()
  
  # Append to the list
  cogannos[[file]] <- reformatted
}

# Combine all results into one data frame
allgroups <- bind_rows(cogannos)


# drop all NAs in group/gene columns
allgroups <- allgroups %>%
  mutate(group = str_replace_na(group, "drop")) %>%
  filter(group != "drop")

# fix the gene issue that happeend for COG2948

allgroups <- allgroups %>%
  filter(gene == "COG:COG2948") %>%
  mutate(gene = "virB10") %>%
  mutate(COG = "COG2948") %>%
  mutate(Letter = "U") %>%
  mutate(Function = "Intracellular trafficking, secretion, and vesicular transport") %>%
  mutate(Category = "CELLULAR PROCESSES AND SIGNALING") %>%
  unique() %>%
  rbind(allgroups %>%
          filter(gene != "COG:COG2948")) %>%
  filter(!str_detect(gene, "SRR"))

# acccount for combination categoruies 
allgroups <- allgroups %>%
  mutate(Function = str_replace_na(Function, "Combination")) %>%
  mutate(Category = str_replace_na(Category, "COMBINATION"))



# read in gene lists for reference as to what was annotated and what wasnt
# Directory containing the .cog files
ref_dir <- "functionalanno/genelists"

# List all files ending in _4.cog
ref_files <- list.files(ref_dir, pattern = ".txt", full.names = TRUE)

# Placeholder for table list
reflists <- list()

# Loop over each file
for (file in ref_files) {
  
  file_name <- basename(file)
  group <- str_remove(file_name, "genes.txt")
  tab <- read.delim(file, sep = "\t", header = FALSE)
  
  reformatted <- tab %>%
    mutate(gene = V1) %>%
    mutate(group = group) %>%
    select(-c(V1))
  
  # Append to the list
  reflists[[file]] <- reformatted
}

# do this but filtering to remove the "group####" genes 
# Combine all results into one data frame
allrefs <- bind_rows(reflists)


# Placeholder for table list
reffiltlists <- list()

# Loop over each file
for (file in ref_files) {
  
  file_name <- basename(file)
  group <- str_remove(file_name, "genes.txt")
  tab <- read.delim(file, sep = "\t", header = FALSE)
  
  reformatted <- tab %>%
    mutate(gene = V1) %>%
    mutate(group = group) %>%
    filter(!str_detect(gene, "group")) %>%
    select(-c(V1))
  
  # Append to the list
  reffiltlists[[file]] <- reformatted
}

# Combine all results into one data frame
allrefsfilt <- bind_rows(reffiltlists)


# ----- Check for missing genes and possible duplicates
dropped <- allrefs %>%
  filter(!gene %in% allgroups$gene)

spec <- dropped %>%
  filter(!str_detect(gene, "group"))


dupes <- allrefsfilt %>%
  # Create a new column with the first four characters of the target column
  mutate(prefix = substr(gene, 1, 4)) %>%
  # Group by this prefix
  group_by(prefix) %>%
  # Filter groups with more than one (for duplicates) or up to three (for triplicates) rows
  filter(n() >= 2 & n() <= 3) %>%
  # Remove the temporary 'prefix' column if not needed in the result
  ungroup() %>%
  select(-prefix)

dupes2 <- dupes %>%
  mutate(gene2 = str_split_fixed(gene, "_", 2)[,1]) %>%
  select(group, gene2) %>%
  unique()

allrefsfilt %>%
  filter(!gene %in% dupes$gene)

# Plot comparisons -----

allgroups <- allgroups %>%
  mutate(source = case_when(
    str_detect(group, "bov") ~ "Clinical Bovine",
    str_detect(group, "hum") ~ "Clinical Human",
    str_detect(group, "env") ~ "Environmental"
  )) %>%
  mutate(level = case_when(
    str_detect(group, "core") ~ "Core",
    str_detect(group, "cloud") ~ "Cloud",
    str_detect(group, "shell") ~ "Shell"
  ))

# Plot Figure 3D

ggplot(allgroups %>%
         group_by(gene) %>%
         filter(n() == 1) %>%
         ungroup() %>%
         filter(Function != "Combination") %>%
         #mutate(Function = str_replace_all(Function, " and ", " + ")) %>%
         mutate(Category = str_to_title(Category)) %>%
         mutate(Function = str_replace_all(Function, "and", "+")), aes(y = fct_infreq(Function), fill = level)) +
  geom_bar(color = "white")+
  scale_x_continuous(expand = c(0,0.5))+
  #scale_y_discrete(expand = c(0,0))+
  theme_classic() +
  theme(legend.position = "bottom",
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black", fill = "white"),
        strip.placement = "inside",
        strip.text.y = element_text(angle = 90, size = 10, color = "black"),
        strip.text.x.top = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 9, color = "black"),
        axis.title = element_text(size=9, color = "black"),
        legend.text = element_text(size = 2, color = "black"),
        legend.title = element_blank()) +
  scale_fill_viridis_d(option = "A", begin = 0.1, end = 0.7) +
  ylab("COG Function") +
  xlab("Number of Unique Genes") +
  facet_grid(Category~source, scales = "free_y", switch = "y", space = "free_y")
