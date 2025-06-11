# Assembly QC 

# ---- Load Packages ----
library("tidyverse")
library("tidylog")

setwd("./R")

# ---- Load Tables ----
quast <- read.delim("assembqc/report.tsv", sep = "\t")
checkm1 <- read.delim("assembqc/omics_checkm1.tsv", sep = "\t") 
checkm2 <- read.delim("assembqc/omics_checkm2.tsv", sep = "\t")
checkm3 <- read.delim("assembqc/omics_checkm3.tsv", sep = "\t")
ss <- read.delim("assembqc/seqseroentero_summary.tsv", sep = "\t")
sr <- read.delim("assembqc/sistrentero_output.txt", sep = "\t")
st <- read.delim("assembqc/mlst.tsv", sep = "\t", header = FALSE) # something here is messed up so circle back 


# ---- Format to combine ----
# CheckM
colnames(checkm1)

checkmall <- rbind(checkm1 %>%
                     mutate(ID=str_split_fixed(Bin.Id, "\\_", 2)[,1]) %>% # cut _assembly off name
                     select(ID, Completeness, Contamination, Strain.heterogeneity), # select most relevant columns,
                   checkm2 %>%
                     mutate(ID=str_split_fixed(Bin.Id, "\\_", 2)[,1]) %>% # cut _assembly off name
                     select(ID, Completeness, Contamination, Strain.heterogeneity), # select most relevant columns,
                   checkm3 %>%
                     mutate(ID=str_split_fixed(Bin.Id, "\\_", 2)[,1]) %>% # cut _assembly off name
                     select(ID, Completeness, Contamination, Strain.heterogeneity)) # select most relevant columns


# QUAST 
colnames(quast)
quast2<- data.table::transpose(quast)
quast2 <- as.data.frame(apply(quast2, 2, as.numeric))
colnames(quast2) <- quast$Assembly
quast2$ID <- colnames(quast)
quast2$ID <- str_remove_all(quast2$ID, "_assembly")
rownames(quast2) <- quast2$ID
quast2 <- quast2[-1,]

# SISTR
colnames(sr)

sr2 <- sr %>%
  mutate(ID=str_split_fixed(genome, "\\_", 2)[,1]) %>% # cut _assembly off name
  select(ID, o_antigen, h1, serogroup, serovar)


# SeqSero2 
colnames(ss)

ss2 <- ss %>%
  mutate(ID=Sample.name) %>%
  filter(!str_detect(ID, "alleles")) %>%
  select(ID, Predicted.antigenic.profile, Predicted.serotype)

# MLST
colnames(st)

st2 <- st %>%
  mutate(ID=str_split_fixed(V1, "\\_", 2)[,1]) %>%
  mutate(ST=V3) %>%
  select(ID, ST)

all <- checkmall %>%
  inner_join(quast2, by = "ID") %>%
  inner_join(sr2, by = "ID") %>%
  inner_join(ss2, by = "ID") %>%
  full_join(st2, by = "ID")

#save this 
#write.table(all, "assembly/qc/allcombinedentero.txt", sep = "\t")

# ----- Filter : Sero -----

# First by serotype
table(all$serovar)
table(all$Predicted.serotype)

notdub <- all %>%
  filter(serovar != "Dublin") # first drop the clearly not Dublin

all %>%
  filter(serovar == "Dublin") %>%
  group_by(Predicted.serotype) %>%
  summarise(count=n())  # SeqSero2 has issues with the O ag 

all %>%
  filter(serovar == "Dublin") %>%
  filter(Predicted.serotype != "Dublin") %>% # check those without O ag
  group_by(o_antigen) %>% # against sistr o ag
  summarise(count=n())  # sistr doesnt have any issues  so keep all 

all %>%
  filter(serovar == "Dublin") %>%
  group_by(ST) %>% # check ST
  summarise(count=n()) 

# sub set these and proceed 
dub <- all %>%
  filter(serovar == "Dublin")

# final dataset

# ----  Filter: Assemb Stats ----

summary(dub %>%
          select(Completeness, Contamination, Strain.heterogeneity)) # these are all fine

summary(dub %>%
          select(`Total length (>= 0 bp)`, `GC (%)`, N50, `# contigs (>= 0 bp)`)) # these are all fine

# no filtering is needed based on these criteria 

# Get details on dropped strains

#which non dublin seros were there?
notdub %>%
  group_by(serovar) %>%
  summarise(count=n())

# which ST types are represented by the dublin seros 
dub %>%
  group_by(ST) %>% 
  summarise(count=n())

