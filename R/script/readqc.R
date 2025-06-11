# Check Read QC 


# ---- Load Packages ----
library("tidyverse")
library("tidylog")
library("readxl")
library("ggpubr")

setwd("./R")

# ---- Load Tables ----
# read stats
raw <- read.delim("readqc/rawstats.txt", sep = "\t")
trim <- read.delim("readqc/trimstats.txt", sep = "\t")
filt <- read.delim("readqc/enterofiltstats.txt", sep = "\t")

# add full stats from separate study where additional strains were included in this dataset 
dub <- read.delim("readqc/readloss_throughk2filt.txt", sep = "\t")

# classification results
tab <- readxl::read_xlsx("readqc/combinedkreports.xlsx", sheet = "combinedkreports")
key <- read_xlsx("readqc/combinedkreports.xlsx", sheet = "IDs_cleaned", col_names = FALSE)

# ---- Check Classfication ----
# subset table 
tab_all <- tab %>%
  select(name, taxid, lvl_type, `#perc`, contains("all"))

# rename columns to strain ID based on key
# as a check before just overwriting names
namesall <- data.frame(original = colnames(tab_all)[6:ncol(tab_all)])

colnames(key) <- c("col1", "id")

names <- namesall %>%
  mutate(col1 = str_split_fixed(original, "_", 2)[,1]) %>%
  inner_join(key %>%
               mutate(col1 = as.character(col1)), by = "col1")

colnames(tab_all)[6:ncol(tab_all)] <- paste0(names$id, "_all")

# transpose
all.t <- data.frame(t(tab_all))
colnames(all.t) <- all.t[1,]
all.t2 <- all.t[6:nrow(all.t),]

sal2 <- all.t2[,1:9]

sal2 <- as.data.frame(apply(sal2, 2, function(x) gsub('\\s+', '', x)))

sal2[] <- lapply(sal2,as.numeric)

sal2 <- sal2 %>%
  mutate(total = root + unclassified) %>%
  mutate(ppnsal = Salmonella/total) %>%
  mutate(ppn_noclass = unclassified/total)

check <- sal2 %>%
  filter(ppnsal < 0.90)

check$id <- rownames(check)

check_all <- tab_all %>%
  filter(lvl_type == "G") %>%
  select(name, lvl_type, all_of(check$id)) 

check_all.m <- reshape2::melt(check_all)

sal2$variable <- rownames(sal2)

check_all.m <- merge(sal2 %>%
                       select(variable, total),
                     check_all.m, by = "variable")

unique(check_all.m$variable)

# subset by strain
bystrain <- list()

for (i in check_all.m$variable) {
  df <- check_all.m[check_all.m$variable == i, ]
  bystrain[[i]] <- df
}

bystrain[[i]] # check

top_ppns_list <- list() # make list

# loop to get the top genera in each by ppn and combine into one table
for (i in check_all.m$variable){
  result <- bystrain[[i]] %>%
    filter(value > 0) %>%
    mutate(ppn = value/total) %>%
    filter(ppn > 0.1) %>%
    unique()
  
  # add the result to the list
  top_ppns_list[[i]] <- result
}

# combine all results into a single data frame
top.ppns <- do.call(rbind, top_ppns_list) 

# check for contaminants
top.ppns %>%
  filter(!str_detect(name, "Salmonella")) 



 # ---- Check Read Loss and Final Average Coverage ---- 

# combine read stats from all
raw <- raw %>%
  mutate(raw_sumlen = sum_len) %>%
  mutate(raw_numseqs = num_seqs) %>%
  select(ID, read, raw_numseqs, raw_sumlen)

trim <- trim %>%
  mutate(ID = str_split_fixed(file, "\\_", 2)[,1]) %>%
  filter(str_detect(file, "P")) %>%
  mutate(read = case_when(
    str_detect(file, "_1P.fq.gz") ~ "forward",
    str_detect(file, "_2P.fq.gz") ~ "reverse"
  )) %>%
  mutate(trim_sumlen = sum_len) %>%
  mutate(trim_numseqs = num_seqs) %>%
  select(ID, read, trim_numseqs, trim_sumlen)

filt2 <- filt %>%
  mutate(ID = str_split_fixed(file, "\\.", 2)[,1]) %>%
  mutate(read = case_when(
    str_detect(file, "_1P.fq") ~ "forward",
    str_detect(file, "_2P.fq") ~ "reverse"
  )) %>%
  mutate(filt_sumlen = sum_len) %>%
  mutate(filt_numseqs = num_seqs) %>%
  select(ID, read, filt_numseqs, filt_sumlen)


all <- raw %>%
  inner_join(trim, by = c("ID", "read")) %>%
  inner_join(filt2, by = c("ID", "read"))

trim %>% filter(!ID %in% filt2$ID) # what was droppped?

all2 <- all %>%
  mutate(rawtotrim = (raw_numseqs-trim_numseqs)/raw_numseqs) %>%
  mutate(trimtofilt = (trim_numseqs-filt_numseqs)/trim_numseqs) %>%
  mutate(rawtofilt = (raw_numseqs-filt_numseqs)/raw_numseqs) %>%
  # add in strains from separate study
  rbind(dub %>%
          filter(str_detect(ID, "NDSU")) %>%
          filter(!ID %in% c("NDSU2", "NDSU5"))) 

# plot read loss distribution 
ggarrange(ggplot(all2 %>%
                   filter(read == "forward"), aes(x = rawtotrim))+
            theme_classic() +
            geom_histogram(bins = 100) +
            theme(text = element_text(size = 12)),
          ggplot(all2 %>%
                   filter(read == "forward"), aes(x = rawtofilt))+
            theme_classic() +
            geom_histogram(bins = 100) +
            theme(text = element_text(size = 12)), nrow = 2)


# check coverage - query only as other study strains meet inclusion criteria
cov <- filt %>%
  mutate(ID = str_split_fixed(file, "\\.", 2)[,1]) %>%
  mutate(read = case_when(
    str_detect(file, "_1P.fq") ~ "forward",
    str_detect(file, "_2P.fq") ~ "reverse"
  )) %>%
  mutate(totalseq = avg_len*num_seqs) %>%
  group_by(ID) %>%
  mutate(totalseq_bothreads = sum(totalseq)) %>%
  ungroup() %>%
  select(ID, totalseq_bothreads) %>%
  unique() %>%
  mutate(coverage = totalseq_bothreads/4857000)

# how many are under 30x coverage 
nrow(cov %>%
       filter(coverage <=30)) # 25 

# what strains are they 
cov %>%
  filter(coverage <=30) %>%
  select(ID) %>%
  unique()
