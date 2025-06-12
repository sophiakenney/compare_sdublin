# Genetic and SNP Distances

# ---- Load Packages ----
library("tidyverse")
library("tidylog")
library("readxl")
library("ggpubr")
library("ape")
library("chorddiag")

setwd("./R")

sourcepal <- c("Clinical Bovine" = "#382A54FF",
               "Environmental" = "#3497A9FF",
               "Clinical Human" = '#395D9CFF')

# ---- Figure 4B,D ----

# Load tree for phylogenetic distances

tree <- ape::read.tree("phylo/iq.treefile")
cophenetic_matrix <- cophenetic(tree)
met <- read.delim("meta/finalset_meta.tsv", sep = "\t")

tab2 <- reshape2::melt(cophenetic_matrix)

filt2 <- tab2 %>%
  mutate(sra_accession = Var1,
         variable = Var2) %>%
  select(sra_accession, variable, value) %>%
  # filter duplicate combinations
  rowwise() %>%
  mutate(pair = paste(sort(c(sra_accession, variable)), collapse = "_")) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair) %>%
  # filter same strain comparisons
  filter(sra_accession != variable)


filt2 <- filt2 %>%
  #add strain 1 meta
  inner_join(met %>%
               mutate(source_1 = source) %>%
               select(sra_accession, source_1), by = "sra_accession") %>%
  # add strain 2 meta
  inner_join(met %>%
               mutate(source_2 = source,
                      variable = sra_accession) %>%
               select(variable, source_2), by = "variable") %>%
  mutate(source_1 = case_when(
    source_1 == "bovclin" ~ "Clinical Bovine",
    source_1 == "bovenv" ~ "Environmental",
    source_1 == "humall" ~ "Clinical Human",
  )) %>%
  mutate(source_2 = case_when(
    source_2 == "bovclin" ~ "Clinical Bovine",
    source_2 == "bovenv" ~ "Environmental",
    source_2 == "humall" ~ "Clinical Human",
  ))


# plot violin Figure 4B
ggplot(filt2 %>%
         filter(source_1 == source_2),
       aes(x=source_2, y = value, fill = source_2)) +
  geom_violin(alpha=0.7) +
  theme_classic() +
  scale_fill_manual(values = sourcepal) +
  theme(axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 14, color="black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.position = "none") +
  ylab("Genetic Distances") +
  xlab("")

# plot density Figure 4D

ggplot(filt2 %>%
         mutate(interintra=case_when(
           source_1 == source_2 ~ "Within",
           source_1 != source_2 ~ "Between"
         )),
       aes(x=value, color=interintra, fill=interintra, color=interintra)) +
  geom_density(alpha=.3) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c("Within" = "#D2693CFF",
                               "Between" = "gold")) +
  scale_color_manual(values = c("Within" = "#D2693CFF",
                                "Between" = "gold")) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color="black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  ylab("Density")+
  xlab("Genetic Distances")


# ---- Figure 4C ----

tab <- read.delim("phylo/distances.tab")
snp <- reshape2::melt(tab %>%
                         mutate(sra_accession = snp.dists.0.8.2) %>%
                         select(-c(snp.dists.0.8.2))) %>%
  inner_join(met, by = "sra_accession") %>%
  mutate(source_1 = source) %>%
  select(-c(source)) %>%
  inner_join(met %>%
               mutate(variable = sra_accession) %>%
               mutate(HHS_region_2 = HHS_region,
                      date_2 = date,
                      source_2 = source) %>%
               select(variable, source_2, HHS_region_2, date_2), by = "variable")

snp2 <- snp %>%
  # remove rows where the comparison is between the same single strain
  mutate(same = ifelse(sra_accession == variable, TRUE, FALSE)) %>%
  filter(same == FALSE) %>%
  select(-c(same)) %>%
  #remove rows with redunant combinations so as not to inflate snp difference distribution
  rowwise() %>%
  mutate(pair = paste(sort(c(sra_accession, variable)), collapse = "_")) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)  %>%
  mutate(snpbin = case_when(
    value == 0 ~ "0 SNPs",
    value > 0 & value <= 5 ~ "0 < SNPs ≤ 5",
    value > 5 & value <= 10 ~ "5 < SNPs ≤ 10",
    value > 10 & value <=20 ~ "10 < SNPs ≤ 20",
    value > 20 ~ "SNPS > 20"
  )) %>%
  mutate(freq = 1)



all <- snp2 %>%
  group_by(source_1, source_2) %>%
  mutate(freq_grouped = sum(freq)) %>%
  ungroup() %>%
  select(source_1, source_2, freq_grouped) %>%
  unique() %>%
  mutate(source_1 = case_when(
    source_1 == "bovclin" ~ "Clinical Bovine",
    source_1 == "bovenv" ~ "Environmental",
    source_1 == "humall" ~ "Clinical Human"
  )) %>%
  mutate(source_2 = case_when(
    source_2 == "bovclin" ~ "Clinical Bovine",
    source_2 == "bovenv" ~ "Environmental",
    source_2 == "humall" ~ "Clinical Human"
  )) %>%
  pivot_wider(id_cols = "source_1", names_from= "source_2", values_from = "freq_grouped")

all <- as.matrix(all[,2:ncol(all)])
rownames(all) <- colnames(all)

snp0wide <- snp2 %>%
  filter(value == 0) %>%
  group_by(source_1, source_2) %>%
  mutate(freq_grouped = sum(freq)) %>%
  ungroup() %>%
  select(source_1, source_2, freq_grouped) %>%
  unique() %>%
  mutate(source_1 = case_when(
    source_1 == "bovclin" ~ "Clinical Bovine",
    source_1 == "bovenv" ~ "Environmental",
    source_1 == "humall" ~ "Clinical Human"
  )) %>%
  mutate(source_2 = case_when(
    source_2 == "bovclin" ~ "Clinical Bovine",
    source_2 == "bovenv" ~ "Environmental",
    source_2 == "humall" ~ "Clinical Human"
  )) %>%
  pivot_wider(id_cols = "source_1", names_from= "source_2", values_from = "freq_grouped")

snp0wide[is.na(snp0wide)] <- 0
snp0wide <- as.matrix(snp0wide[,2:ncol(snp0wide)])
rownames(snp0wide) <- colnames(snp0wide)

snp5wide <- snp2 %>%
  filter(value > 0) %>%
  filter(value <= 5) %>%
  group_by(source_1, source_2) %>%
  mutate(freq_grouped = sum(freq)) %>%
  ungroup() %>%
  select(source_1, source_2, freq_grouped) %>%
  unique() %>%
  mutate(source_1 = case_when(
    source_1 == "bovclin" ~ "Clinical Bovine",
    source_1 == "bovenv" ~ "Environmental",
    source_1 == "humall" ~ "Clinical Human"
  )) %>%
  mutate(source_2 = case_when(
    source_2 == "bovclin" ~ "Clinical Bovine",
    source_2 == "bovenv" ~ "Environmental",
    source_2 == "humall" ~ "Clinical Human"
  )) %>%
  pivot_wider(id_cols = "source_1", names_from= "source_2", values_from = "freq_grouped")

snp5wide[is.na(snp5wide)] <- 0
snp5wide <- as.matrix(snp5wide[,2:ncol(snp5wide)])
rownames(snp5wide) <- colnames(snp5wide)

snp10wide <- snp2 %>%
  filter(value > 5) %>%
  filter(value <= 10) %>%
  group_by(source_1, source_2) %>%
  mutate(freq_grouped = sum(freq)) %>%
  ungroup() %>%
  select(source_1, source_2, freq_grouped) %>%
  unique() %>%
  mutate(source_1 = case_when(
    source_1 == "bovclin" ~ "Clinical Bovine",
    source_1 == "bovenv" ~ "Environmental",
    source_1 == "humall" ~ "Clinical Human"
  )) %>%
  mutate(source_2 = case_when(
    source_2 == "bovclin" ~ "Clinical Bovine",
    source_2 == "bovenv" ~ "Environmental",
    source_2 == "humall" ~ "Clinical Human"
  )) %>%
  pivot_wider(id_cols = "source_1", names_from= "source_2", values_from = "freq_grouped")

snp10wide[is.na(snp10wide)] <- 0
snp10wide <- as.matrix(snp10wide[,2:ncol(snp10wide)])
rownames(snp10wide) <- colnames(snp10wide)

snp20wide <- snp2 %>%
  filter(value > 10) %>%
  filter(value <= 20) %>%
  group_by(source_1, source_2) %>%
  mutate(freq_grouped = sum(freq)) %>%
  ungroup() %>%
  select(source_1, source_2, freq_grouped) %>%
  unique() %>%
  mutate(source_1 = case_when(
    source_1 == "bovclin" ~ "Clinical Bovine",
    source_1 == "bovenv" ~ "Environmental",
    source_1 == "humall" ~ "Clinical Human"
  )) %>%
  mutate(source_2 = case_when(
    source_2 == "bovclin" ~ "Clinical Bovine",
    source_2 == "bovenv" ~ "Environmental",
    source_2 == "humall" ~ "Clinical Human"
  )) %>%
  pivot_wider(id_cols = "source_1", names_from= "source_2", values_from = "freq_grouped")

snp20wide[is.na(snp20wide)] <- 0
snp20wide <- as.matrix(snp20wide[,2:ncol(snp20wide)])
rownames(snp20wide) <- colnames(snp20wide)

snpover20wide <- snp2 %>%
  filter(value > 20) %>%
  group_by(source_1, source_2) %>%
  mutate(freq_grouped = sum(freq)) %>%
  ungroup() %>%
  select(source_1, source_2, freq_grouped) %>%
  unique() %>%
  mutate(source_1 = case_when(
    source_1 == "bovclin" ~ "Clinical Bovine",
    source_1 == "bovenv" ~ "Environmental",
    source_1 == "humall" ~ "Clinical Human"
  )) %>%
  mutate(source_2 = case_when(
    source_2 == "bovclin" ~ "Clinical Bovine",
    source_2 == "bovenv" ~ "Environmental",
    source_2 == "humall" ~ "Clinical Human"
  )) %>%
  pivot_wider(id_cols = "source_1", names_from= "source_2", values_from = "freq_grouped")

snpover20wide[is.na(snpover20wide)] <- 0
snpover20wide <- as.matrix(snpover20wide[,2:ncol(snpover20wide)])
rownames(snpover20wide) <- colnames(snpover20wide)


# Plot figures - note that the directionality in these should be visualized as the modified figure in the manuscript
chorddiag(snp0wide, groupColors = c("#382A54FF",'#395D9CFF',"#3497A9FF"), 
          groupnamePadding =10,
          showTicks = FALSE,
          showTooltips = FALSE)
chorddiag(snp5wide, groupColors = c("#382A54FF","#3497A9FF", '#395D9CFF'), 
          groupnamePadding =10,
          showTicks = FALSE,
          showTooltips = FALSE)
chorddiag(snp10wide, groupColors = c("#382A54FF",'#395D9CFF',"#3497A9FF"), 
          groupnamePadding =10,
          showTicks = FALSE,
          showTooltips = FALSE)
chorddiag(snp20wide, groupColors = c("#382A54FF",'#395D9CFF',"#3497A9FF"), 
          groupnamePadding =10,
          showTicks = FALSE,
          showTooltips = FALSE)
f