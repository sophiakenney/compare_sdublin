# Pangenome 

# ---- Load Packages ----
library("tidyverse")
library("tidylog")
library("vegan")
library("pairwiseAdonis")
library("pagoo")

setwd("./R")

# ---- Load Tables ---- 
met <- read.delim("meta/finalset_meta.tsv", sep = "\t")
full <- read.csv("pangenome/raw/gene_presence_absence.csv")
tab <- read.delim("pangenome/raw/gene_presence_absence.Rtab")
bov <- read.csv("pangenome/raw/bovclin_gene_presence_absence.csv")
env <- read.csv("pangenome/raw/bovenv_gene_presence_absence.csv")
hum <- read.csv("pangenome/raw/humall_gene_presence_absence.csv")

# set color palette for plots 

sourcepal <- c("Clinical Bovine" = "#382A54FF",
               "Environmental" = "#3497A9FF",
               "Clinical Human" = '#395D9CFF')

# ---- Pangenome PCA - Figure 3A ---- 

rownames(tab) <- tab$Gene
tab2 <- t(tab[,2:ncol(tab)])
colnames(tab2) <- rownames(tab)
rownames(tab2) <- colnames(tab[,2:ncol(tab)])

r2mat <- as.matrix(tab2)
roarymds<- metaMDS(r2mat, distance = "jaccard", autotransform = FALSE, k= 3, trymax = 50)


#plot
roarymds <- as.data.frame(scores(roarymds)$site)
roarymds$sra_accession <- rownames(roarymds)
roarymds <- merge(roarymds, met %>%
                    select(sra_accession, HHS_region, date, source), by = "sra_accession")

# permanova
set.seed(12345)
roarydist <- vegdist(r2mat, method = "jaccard")


#ADONIS test
# source
vegan::adonis2(roarydist ~ as.factor(roarymds$source), data = roarymds, permutations = 10000)

# pairwise for source
statdf <- pairwise.adonis(roarydist, as.factor(roarymds$source))
statdf

# plot figure 3A
ggplot(roarymds %>%
         mutate(Source = case_when(
           source == "bovclin" ~ "Clinical Bovine",
           source == "bovenv" ~ "Environmental",
           source == "humall" ~ "Clinical Human"
         )), aes(x=NMDS1,y=NMDS2, color = Source)) +
  theme_classic() +
  stat_ellipse(linewidth = 0.5)+
  scale_color_manual(values = sourcepal)+
  geom_point(size=1, alpha = 0.4) +
  theme(axis.title = element_text(size=10, color = "black"),
        axis.text = element_text(size=10, color = "black"),
        legend.title = element_text(size=10, color = "black"),
        legend.text = element_text(size = 10, color = "black"), 
        legend.position = "none")


# ---- Pangenome Permutation Analysis - Figure 3B ---- 

# ----- Permutations
  # create 100 data frames with randomly 30 selected strains
  random_bov <- data.frame()
  listofdfs_bov <- list()
  for (i in 1:100) {
    random_bov <- bov[ , sample(15:595, 30)]
    random_bov_bind <- cbind(bov[1:14], random_bov[1:15])
    listofdfs_bov[[i]] <- random_bov_bind
  }
  
  # create 100 data frames with randomly 30 selected strains
  random_env <- data.frame()
  listofdfs_env <- list()
  for (i in 1:100) {
    random_env <- env[ , sample(15:919, 30)]
    random_env_bind <- cbind(env[1:14], random_env[1:15])
    listofdfs_env[[i]] <- random_env_bind
  }
  
  # create 100 data frames with randomly 30 selected strains
  random_hum <- data.frame()
  listofdfs_hum <- list()
  for (i in 1:100) {
    random_hum <- hum[ , sample(15:678, 30)]
    random_hum_bind <- cbind(hum[1:14], random_hum[1:15])
    listofdfs_hum[[i]] <- random_hum_bind
  }
  
  # save dataframes as csv files
  for (i in 1:100) {
    write.csv(listofdfs_bov[[i]], paste0("pangenome/subsets/random_bov_bind", i, ".csv"), row.names = FALSE)
  }
  
  # save dataframes as csv files
  for (i in 1:100) {
    write.csv(listofdfs_env[[i]], paste0("pangenome/subsets/random_env_bind", i, ".csv"), row.names = FALSE)
  }
  
  # save dataframes as csv files
  for (i in 1:100) {
    write.csv(listofdfs_hum[[i]], paste0("pangenome/subsets/random_hum_bind", i, ".csv"), row.names = FALSE)
  }
  
  # run power law model on 100 gene presence/absence matrices with randomly selected cattle isolates
  pg_bov <- list()
  powerlaw_bov <- list()
  for (i in 1:100) {
    pg_bov[[i]] <- roary_2_pagoo(gene_presence_absence_csv = paste0("pangenome/subsets/random_bov_bind", i, ".csv"), sep = ",")
    powerlaw_bov[[i]] <- pg_bov[[i]]$pg_power_law_fit()
  }
  
  # run power law model on 100 gene presence/absence matrices with randomly selected environmental isolates
  pg_env <- list()
  powerlaw_env <- list()
  for (i in 1:100) {
    pg_env[[i]] <- roary_2_pagoo(gene_presence_absence_csv = paste0("pangenome/subsets/random_env_bind", i, ".csv"), sep = ",")
    powerlaw_env[[i]] <- pg_env[[i]]$pg_power_law_fit()
  }
  
  # run power law model on 100 gene presence/absence matrices with randomly selected human isolates
  pg_hum <- list()
  powerlaw_hum <- list()
  for (i in 1:100) {
    pg_hum[[i]] <- roary_2_pagoo(gene_presence_absence_csv = paste0("pangenome/subsets/random_hum_bind", i, ".csv"), sep = ",")
    powerlaw_hum[[i]] <- pg_hum[[i]]$pg_power_law_fit()
  }
  
  # create vector with alpha values from calculated power law models (bov isolates)
  random_alphas_bov <- vector()
  for (i in 1:100) {
    random_alphas_bov[[i]] <- attr(powerlaw_bov[[i]], which = "alpha")
  }
  hist(random_alphas_bov)
  qqnorm(random_alphas_bov)
  mean_bov <- mean(random_alphas_bov)
  sd_bov <- sd(random_alphas_bov)
  
  # create vector with alpha values from calculated power law models (env isolates)
  random_alphas_env <- vector()
  for (i in 1:100) {
    random_alphas_env[[i]] <- attr(powerlaw_env[[i]], which = "alpha")
  }
  hist(random_alphas_env)
  qqnorm(random_alphas_env)
  mean_env <- mean(random_alphas_env)
  sd_env <- sd(random_alphas_env)
  
  # create vector with alpha values from calculated power law models (human isolates)
  random_alphas_hum <- vector()
  for (i in 1:100) {
    random_alphas_hum[[i]] <- attr(powerlaw_hum[[i]], which = "alpha")
  }
  hist(random_alphas_hum)
  qqnorm(random_alphas_hum)
  mean_hum <- mean(random_alphas_hum)
  sd_hum <- sd(random_alphas_hum)
  
  # create a dataframe containing random alphas
  random_alphas_df <- data.frame(random_alphas_bov, random_alphas_env, random_alphas_hum)
  write.csv(random_alphas_df,"pangenome/subsets/random_alphas_source.csv")
  alphas2 <- reshape2::melt(random_alphas_df) %>%
    mutate(alpha = value,
           source = str_remove_all(variable, "random_alphas_")) %>%
    select(alpha, source)
  
  summary(random_alphas_df$random_alphas_bov)
  summary(random_alphas_df$random_alphas_env)
  summary(random_alphas_df$random_alphas_hum)
  
  # run anova on random alphas
  isolationsource_aov <- aov(alpha ~ source, alphas2)
  summary(isolationsource_aov)
  tukey_aov <- TukeyHSD(isolationsource_aov)
  tukey_aov
  tukey.plot.aov<-aov(alpha ~ source, data = alphas2)
  tukey.plot.test <- TukeyHSD(tukey.plot.aov)
  par(mar = c(5, 10, 4, 2)) 
  plot(tukey.plot.test, las = 1)
  
  
  # Convert Tukey HSD output to a data frame
  tukey_df <- as.data.frame(tukey.plot.test$source)
  tukey_df$comparison <- rownames(tukey_df)
  
  # Plot using ggplot2
  ggplot(tukey_df %>%
           mutate(comparison2 = case_when(
             comparison == "hum-env" ~ "Human vs Environmental", 
             comparison == "hum-bov" ~ "Human vs Bovine Clinical", 
             comparison == "env-bov" ~ "Environmental vs Bovine Clinical" 
           )), aes(x = comparison2, y = diff)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +  # Add error bars
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed", size = 1) + 
    coord_flip() +  # Flip coordinates for better readability
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) + 
    labs(
      title = "Tukey Test of Heap's Law Alphas", 
      x = "", 
      y = "Mean Difference"
    )
  
  # Plot Figure 3B
  set.seed(1243)
  ggplot(data = alphas2 %>%
           mutate(source = case_when(
             source == "bov" ~ "Clinical Bovine",
             source == "env" ~ "Environmental",
             source == "hum" ~ "Clinical Human"
           )), aes(x = source, y = alpha))+
    geom_boxplot(aes(fill = source), outliers = FALSE, alpha = 0.9) + # since plotted with jitter
    geom_jitter(color="black", size=1, alpha=0.5, width = 0.2) +
    scale_fill_viridis_d(option = "G", begin = 0.2, end = 0.6) + 
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 10, color = "black"),
          axis.text = element_text(size = 9, color = "black"), 
          legend.position = "none") +
    ylab("Alpha") +
    coord_flip()  
  
# ---- Pangenome Composition - Figure 3C ---- 

ppns <- rbind (bov %>%
                 mutate(ppn = No..isolates/581) %>%
                 mutate(cat = case_when(
                   ppn >= 0.95 ~ "core",
                   ppn <= 0.15 ~ "cloud", 
                   ppn < 0.95 & ppn > 0.15 ~ "shell" # to match pagoo core used to calc
                 )) %>%
                 select(ppn,cat) %>%
                 group_by(cat) %>%
                 count() %>%
                 mutate(source = "bovclin"),
               
               env %>%
                 mutate(ppn = No..isolates/581) %>%
                 mutate(cat = case_when(
                   ppn >= 0.95 ~ "core",
                   ppn <= 0.15 ~ "cloud", 
                   ppn < 0.95 & ppn > 0.15 ~ "shell" # to match pagoo core used to calc
                 )) %>%
                 select(ppn,cat) %>%
                 group_by(cat) %>%
                 count() %>%
                 mutate(source = "bovenv"),
               
               hum %>%
                 mutate(ppn = No..isolates/581) %>%
                 mutate(cat = case_when(
                   ppn >= 0.95 ~ "core",
                   ppn <= 0.15 ~ "cloud", 
                   ppn < 0.95 & ppn > 0.15 ~ "shell" # to match pagoo core used to calc
                 )) %>%
                 select(ppn,cat) %>%
                 group_by(cat) %>%
                 count() %>%
                 mutate(source = "humall")) %>%
  group_by(source) %>%
  mutate(perc = n/sum(n))

ggplot(ppns %>%
         mutate(Category = str_to_sentence(cat))%>%
         mutate(source = case_when(
           source == "bovclin" ~ "Clinical Bovine",
           source == "bovenv" ~ "Environmental",
           source == "humall" ~ "Clinical Human"
         )), aes(y = source, x = perc, fill = Category)) +
  geom_bar(stat = "identity", color="white") +
  #geom_text(size = 10, position = position_stack(vjust = 0.5), color = "white") +
  theme_classic() + 
  scale_fill_viridis_d(option = "A", begin = 0.1, end = 0.7) +
  scale_x_continuous(expand = c(0,0)) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.position = "bottom",
        plot.margin = margin(20,5,5,5)) +
  xlab("Proportion")

# ---- Perform Openness Test/Calc ----

# on all - use permutational analysis for by-source
p <- roary_2_pagoo("pangenome/raw/gene_presence_absence.csv")

colnames(met)[1] <- "org" # change to match name in p object 

# circle back to adding metadata
met <- arrange(met, org)
met$org
p$add_metadata(map="org", met)

#check
p$organisms
met$org


# fit power law function
p$pg_power_law_fit() # for pangenome

# ---- Extract List for Anno ----
p2 <- roary_2_pagoo("pangenome/raw/bovclin_gene_presence_absence.csv")
p3 <- roary_2_pagoo("pangenome/raw/bovenv_gene_presence_absence.csv")
p4 <- roary_2_pagoo("pangenome/raw/humall_gene_presence_absence.csv")

p2$add_metadata(map="org", met %>%
                  filter(org %in% colnames(bov)))
p3$add_metadata(map="org", met %>%
                  filter(org %in% colnames(env)))
p4$add_metadata(map="org", met %>%
                  filter(org %in% colnames(hum)))

#check
p2$organisms
p3$organisms
p4$organisms

# extract lists 

bovcore <- as.data.frame(p2$core_clusters)
envcore <- as.data.frame(p3$core_clusters)
humcore <- as.data.frame(p4$core_clusters)


bovshell <- as.data.frame(p2$shell_clusters)
envshell <- as.data.frame(p3$shell_clusters)
humshell <- as.data.frame(p4$shell_clusters)


bovcloud <- as.data.frame(p2$cloud_clusters)
envcloud <- as.data.frame(p3$cloud_clusters)
humcloud <- as.data.frame(p4$cloud_clusters)


# compare core
corediff <- list(bovonly = lapply(bovcore %>% filter(!cluster %in% humcore$cluster) %>% filter(!cluster %in% envcore$cluster), as.character),
                 humonly = lapply(humcore %>% filter(!cluster %in% bovcore$cluster) %>% filter(!cluster %in% envcore$cluster), as.character),
                 envonly = lapply(envcore %>% filter(!cluster %in% humcore$cluster) %>% filter(!cluster %in% bovcore$cluster), as.character))


# compare cloud
clouddiff <- list(
  bovonly = lapply(bovcloud %>% filter(!cluster %in% humcloud$cluster) %>% filter(!cluster %in% envcloud$cluster), as.character),
  humonly = lapply(humcloud %>% filter(!cluster %in% bovcloud$cluster) %>% filter(!cluster %in% envcloud$cluster), as.character),
  envonly = lapply(envcloud %>% filter(!cluster %in% humcloud$cluster) %>% filter(!cluster %in% bovcloud$cluster), as.character))

# compare shell
shelldiff <- list(bovonly = lapply(bovshell %>% filter(!cluster %in% humshell$cluster) %>% filter(!cluster %in% envshell$cluster), as.character),
                  humonly = lapply(humshell %>% filter(!cluster %in% bovshell$cluster) %>% filter(!cluster %in% envshell$cluster), as.character),
                  envonly = lapply(envshell %>% filter(!cluster %in% humshell$cluster) %>% filter(!cluster %in% bovshell$cluster), as.character))

bovshell <- shelldiff$bovonly$cluster
bovcloud <- clouddiff$bovonly$cluster
bovcore <- corediff$bovonly$cluster

humshell <- shelldiff$humonly$cluster
humcloud <- clouddiff$humonly$cluster
humcore <- corediff$humonly$cluster

envshell <- shelldiff$envonly$cluster
envcloud <- clouddiff$envonly$cluster
envcore <- corediff$envonly$cluster

write.table(bovshell, "functionalanno/genelists/bovshellgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(bovcloud, "functionalanno/genelists/bovcloudgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(bovcore, "functionalanno/genelists/bovcoregenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(humshell, "functionalanno/genelists/humshellgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(humcloud, "functionalanno/genelists/humcloudgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(humcore, "functionalanno/genelists/humcoregenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(envshell, "functionalanno/genelists/envshellgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(envcloud, "functionalanno/genelists/envcloudgenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(envcore, "functionalanno/genelists/envcoregenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# create gene (all not just unique) for these

write.table(bov %>% select(Gene), "functionalanno/genelists/allbovgenes.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(env %>% select(Gene), "functionalanno/genelists/allenvgenes.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(hum %>% select(Gene), "functionalanno/genelists/allhumgenes.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


# ---- Table 5 Summary Counts ----
p$summary_stats
p2$summary_stats
p3$summary_stats
p4$summary_stats
