# Abricate: Virulence Factors and Plasmids

# ---- Load Packages ----
library("tidyverse")
library("tidylog")
library("gt")

setwd("./R")

# ---- Load Tables ----
vf <- read.delim("abricate/vf-all_summary.txt", sep = "\t")
pld <- read.delim("abricate/plasmid-all_summary.txt", sep = "\t")
vfdb <- read.delim("abricate/vfdbmeta.txt", sep = "\t")  %>%
  mutate(group = GENE) %>%
  select (-c(GENE))
met <- read.delim("meta/finalset_meta.tsv", sep = "\t")

# format to binary and long

# VF
vf2 <- vf %>%
  mutate(ID=str_remove_all(X.FILE, "_assembly.fasta")) %>%
  select(-c(NUM_FOUND, X.FILE)) %>%
  relocate(ID) %>%
  filter(!str_detect(ID, "ARS")) %>%
  filter(!ID %in% c("NDSU2", "NDSU5"))

vf2[2:ncol(vf2)] <- replace(vf2[2:ncol(vf2)], vf2[2:ncol(vf2)] !=".", 1 ) #replace presence as 1
vf2[2:ncol(vf2)] <- replace(vf2[2:ncol(vf2)], vf2[2:ncol(vf2)] ==".", 0 ) #replace absence as 0

vf2 <- cbind(vf2$ID, mutate_all(vf2[2:ncol(vf2)], function(x) as.numeric(as.character(x))))
colnames(vf2)[1] <- "ID"

vf.m <- reshape2::melt(vf2)

# Plasmid
pld2 <- pld %>%
  mutate(ID=str_remove_all(X.FILE, "_assembly.fasta")) %>%
  select(-c(NUM_FOUND, X.FILE)) %>%
  relocate(ID) %>%
  filter(!str_detect(ID, "ARS")) %>%
  filter(!ID %in% c("NDSU2", "NDSU5"))


pld2[2:ncol(pld2)] <- replace(pld2[2:ncol(pld2)], pld2[2:ncol(pld2)] !=".", 1 ) #replace presence as 1
pld2[2:ncol(pld2)] <- replace(pld2[2:ncol(pld2)], pld2[2:ncol(pld2)] ==".", 0 ) #replace absence as 0

pld2 <- cbind(pld2$ID, mutate_all(pld2[2:ncol(pld2)], function(x) as.numeric(as.character(x))))
colnames(pld2)[1] <- "ID"

pld.m <- reshape2::melt(pld2)

# ---- Virulence Analysis ----

# summary counts
cts <- vf %>%
  mutate(ID=str_remove_all(X.FILE,"_assembly.fasta")) %>%
  mutate(vfct = NUM_FOUND) %>%
  select(ID, vfct) %>%
  filter(!str_detect(ID, "ARS")) %>%
  filter(!ID %in% c("NDSU2", "NDSU5"))

summary(cts$vfct)

# add database metadata
vfdb <- vfdb %>%
  filter(group %in% vf.m$variable)

vf.m <- left_join(vfdb %>%
                    select(-c(strain_assoc, protein_fxn)),
                  vf.m %>%
                    filter(variable %in% vfdb$group)%>%
                    mutate(group = variable) %>%
                    select(-c(variable)), by = "group") %>%
  unique()

vf.genect <- vf.m %>%
  filter(value!=0) %>%
  group_by(group) %>%
  mutate(n = n()) %>%
  select(category, subcategory, group, n) %>%
  unique() 

# Table 1
gt(vf.genect %>%
     mutate(Category = subcategory) %>%
   group_by(Category) %>%
  count() %>%
    ungroup())


# ---- Plasmid ----

pldct <-pld.m %>%
  filter(value!=0) %>%
  group_by(variable) %>%
  mutate(n = n()) %>%
  select(variable, n) %>%
  unique()

pldctsbysource <- pld.m %>%
  filter(value ==1) %>%
  mutate(sra_accession = ID) %>%
  inner_join(met %>% select(sra_accession, HHS_region, source), by = "sra_accession") %>%
  group_by(variable, source) %>%
  mutate(bysource = n()) %>%
  ungroup() %>%
  group_by(variable, HHS_region) %>%
  mutate(byhhs = n()) %>%
  ungroup() %>%
  select(variable, source, bysource, HHS_region, byhhs) %>%
  unique() %>%
  arrange(source, HHS_region) %>%
  select(variable, bysource, source) %>%
  mutate(ppn = case_when(
    source == "bovclin" ~ round(bysource/581, digits = 3)*100,
    source == "bovenv" ~ round(bysource/905, digits = 3)*100,
    source == "humall" ~ round(bysource/664, digits = 3)*100
  )) %>%
  mutate(Source = case_when(
    source == "bovclin" ~ "Clinical Bovine",
    source == "bovenv" ~ "Environmental",
    source == "humall" ~ "Clinical Human"
  )) %>%
  unique()

pldctsbysource2 <- pldctsbysource %>%
  filter(Source == "Clinical Bovine") %>%
  mutate(`Clinical Bovine` = ppn) %>%
  mutate(CB_ct = bysource) %>%
  select(variable, `Clinical Bovine`, CB_ct) %>%
  full_join(pldctsbysource %>%
              filter(Source == "Clinical Human") %>%
              mutate(CH_ct = bysource) %>%
              mutate(`Clinical Human` = ppn) %>%
              select(variable, `Clinical Human`, CH_ct), by = "variable") %>%
  full_join(pldctsbysource %>%
              filter(Source == "Environmental") %>%
              mutate(E_ct = bysource) %>%
              mutate(Environmental = ppn) %>%
              select(variable, Environmental, E_ct), by = "variable")

pldctsbysource2[is.na(pldctsbysource2)] <- 0


# Table 4 - without statistical comparisons included
gt(pldctsbysource2 %>%
     mutate(variable_clean = str_replace_all(variable, "\\.", "/")) %>%
     mutate(variable_clean = str_remove_all(variable_clean, "_1")) %>%
     mutate(variable_clean = case_when(
       variable_clean == "IncFII/S/" ~ "IncFII(S)",
       variable_clean == "IncFIB/AP001918/" ~ "IncFIB(AP001918)",
       variable_clean == "IncFIB/pB171/_pB171" ~ "IncFIB(pB171)",
       variable_clean == "IncHI1B/R27/_R27" ~ "IncHI1B(R27)",
       variable_clean == "IncI1_Alpha" ~ "IncI1(alpha)",
       variable_clean == "IncFIA/HI1/_HI1" ~ "IncFIA(HI1)",
       variable_clean == "IncI2_Delta" ~ "IncI2(delta)",
       variable_clean == "Col/BS512/" ~ "Col(BS512)",
       variable_clean == "Col/MG828/" ~ "Col(MG828)",
       variable_clean == "IncFII/pRSB107/_pRSB107" ~ "IncFII(pRSB107)",
       TRUE ~ variable_clean
     ))%>%
     select(-variable) %>%
     select(variable_clean, `Clinical Bovine`, CB_ct, `Clinical Human`, CH_ct, Environmental, E_ct) %>%
     group_by(variable_clean) %>%
     mutate(total = CB_ct + CH_ct + E_ct) %>%
     mutate(total_ppn = round((total/2150)*100, digits = 1)) %>%
     arrange(desc(total)) %>%
     ungroup() %>%
     mutate(Plasmid = variable_clean) %>%
     mutate(`Clinical Bovine (N)` = paste0(`Clinical Bovine`, "% (", CB_ct, ")")) %>%
     mutate(`Clinical Human (N)` = paste0(`Clinical Human`, "% (", CH_ct, ")")) %>%
     mutate(`Environmental (N)` = paste0(Environmental, "% (", E_ct, ")")) %>%
     mutate(`Total (N)` = paste0(total_ppn, "% (", total, ")")) %>%
     select(Plasmid, `Clinical Bovine (N)`, `Clinical Human (N)`, `Environmental (N)`,`Total (N)`))


# statistical analysis

# subset plasmids present in at least 5 strains of each group
pldctsbysource3 <- pldctsbysource2 %>% filter(CB_ct >= 5 & CH_ct >= 5, E_ct >=5) %>%
  mutate(CB_total = 581,
         CH_total = 664, 
         E_total = 905) %>%
  mutate(variable = factor(variable, levels=c("IncFII.S._1" ,"IncX1_1" ,"ColRNAI_1","IncA.C2_1"))) %>%
  select(variable, CB_ct, CB_total, CH_ct, CH_total, E_ct, E_total) %>%
  ungroup()


# Loop through each plasmid and perform pairwise prop test

results_list <- list()

# Loop through each row of `pldctsbysource3`
for (i in 1:nrow(pldctsbysource3)) {
  # Extract occurrences and sample sizes for each group
  occurrences <- as.numeric(pldctsbysource3[i, c("CB_ct", "CH_ct", "E_ct")])
  sample_sizes <- as.numeric(pldctsbysource3[i, c("CB_total", "CH_total", "E_total")])
  
  # Perform pairwise proportion test
  test_result <- pairwise.prop.test(x = occurrences, n = sample_sizes, p.adjust.method = "none")
  
  # Extract results as a data frame
  test_df <- data.frame(
    Variable = pldctsbysource3$variable[i],            # Name of the variable
    Comparison = rownames(test_result$p.value),       # Pairwise comparisons (e.g., CB vs CH)
    p.value = test_result$p.value                     # P-values for comparisons
  )
  
  # Add the data frame to the list
  results_list[[i]] <- test_df
}

# Combine all results into a single data frame
final_results <- do.call(rbind, results_list)


# add the actual variable 
final_results <- reshape2::melt(final_results %>%
                                  mutate(Comparison = case_when(
                                    Comparison == "2"~ "Clinical Human",
                                    Comparison == "3" ~ "Environmental"
                                  )) %>%
                                  mutate(`Clinical Bovine` = p.value.1) %>%
                                  mutate(`Clinical Human` = p.value.2) %>%
                                  select(-c("p.value.1", "p.value.2"))) %>%
  filter(value > 0)



final_results$bonferroni_adjusted <- p.adjust(final_results$value, method = "bonferroni")

# check for significant differences
final_results %>%
  filter(bonferroni_adjusted < 0.05) %>%
  arrange(Variable) %>%
  mutate(pval = round(bonferroni_adjusted, digits = 4))

