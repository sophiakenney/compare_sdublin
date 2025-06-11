# AMRFinder 

# ---- Load Packages ----
library("tidyverse")
library("tidylog")
library("gt")
library("vegan")
library("pairwiseAdonis")

setwd("./R")

# ---- Load Tables ----
tab <- read.delim("amrfinder/amrfinderentero.txt", sep = "\t")
met <- read.delim("meta/finalset_meta.tsv", sep = "\t")

# ---- Table 2 (not manually rearranged for manuscript) ----

gt(tab %>%
     filter(!Element.type %in% c("VIRULENCE")) %>%
     filter(Element.subtype != "ACID") %>%
     select(filename, Gene.symbol, Element.type, Class) %>%
     unique() %>%
     group_by(Gene.symbol, Class) %>%
     mutate(N_strains = n()) %>%
     mutate(ppn = round((N_strains/2150)*100, digits = 1)) %>%
     mutate(ppn = paste0(ppn, "%")) %>%
     ungroup() %>%
     mutate(Class = str_to_sentence(Class)) %>%
     mutate(Class = ifelse(str_detect(Class, "Quin"), "Quinolone", Class)) %>%
     select(Class, Gene.symbol, N_strains, ppn) %>%
     group_by(Class) %>%
     arrange(.by_group = TRUE) %>%
     unique()) %>%
  cols_label(
    N_strains = "Number of isolates") %>%
  cols_label(
    ppn = "Positive rate (%)") %>%
  cols_label(
    Gene.symbol = "Resistance Gene") %>%
  tab_style(
    style = cell_text(align = "left", weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = cell_text(align = "center", weight = "bold"),
    locations = cells_column_labels(columns = c(Gene.symbol, Class))
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_body(columns = c(N_strains, ppn))
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = c(Gene.symbol))
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman"),
    locations = cells_body()
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = cell_borders(sides = c("left", "right"), color = "lightgrey", weight = px(1)),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_borders(sides = c("left", "right"), color = "lightgrey", weight = px(1)),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_borders(sides= c( "top", "bottom"), color = "lightgrey", weight = px(1)),
    locations = cells_body(columns = everything())
  ) %>%
  opt_vertical_padding(scale = 0)

# ---- Table 3: Statistical Comparisons ----

# Pairwise proportion test for genes in >5% of all strains

tab2 <- tab %>%
  mutate(sra_accession = filename) %>%
  filter(Element.type == "AMR") %>%
  select(sra_accession, Contig.id, Strand, Stop, Start, Gene.symbol) %>%
  inner_join(met %>% select(source, sra_accession), by = "sra_accession") %>%
  group_by(source, Gene.symbol) %>%
  mutate(relabun = case_when(
    source == "bovclin" ~ n()/581,
    source == "bovenv" ~ n()/905,
    source == "humall" ~ n()/664
  )) %>%
  ungroup() %>%
  group_by(Gene.symbol) %>%
  mutate(totalrelabun = n()/2150) %>%
  ungroup() %>%
  filter(!Gene.symbol %in% c("mdsB", "mdsA")) %>%
  filter(totalrelabun >=0.05) %>%
  select(Gene.symbol, source, relabun, totalrelabun) %>%
  unique() %>%
  mutate(count = case_when(
    source == "bovclin" ~ relabun*581,
    source == "bovenv" ~ relabun*905,
    source == "humall" ~ relabun*664
  )) %>%
  mutate(total = case_when(
    source == "bovclin" ~ 581,
    source == "bovenv" ~ 905,
    source == "humall" ~ 664
  ))

# group by gene and perform chisq.test
global_chi2 <- tab2 %>%
  group_by(Gene.symbol) %>%
  summarise(
    pval_chi = {
      counts <- count
      totals <- total
      mat <- rbind(counts, totals - counts)
      test <- chisq.test(mat)
      test$p.value
    },
    .groups = "drop"
  ) %>%
  mutate(BH_adjusted_chi = p.adjust(pval_chi, method = "bonferroni"))


# generate all source pairs
source_pairs2 <- combn(unique(tab2$source), 2, simplify = FALSE)

# for each gene, run pairwise prop.test
pairwise_results2 <- tab2 %>%
  group_by(Gene.symbol) %>%
  filter(n() >= 2) %>%
  summarise(
    pairwise = list({
      gene_data <- cur_data()
      map_dfr(source_pairs2, function(pair) {
        g1 <- gene_data %>% filter(source == pair[1])
        g2 <- gene_data %>% filter(source == pair[2])
        
        if (nrow(g1) == 1 & nrow(g2) == 1) {
          r1 <- g1$count
          r2 <- g2$count
          t1 <- g1$total
          t2 <- g2$total
          
          test <- prop.test(c(r1, r2), c(t1, t2), alternative = "two.sided", correct = TRUE)
          
          tibble(
            source1 = pair[1],
            source2 = pair[2],
            pval = test$p.value
          )
        } else {
          tibble(source1 = pair[1], source2 = pair[2], pval = NA)
        }
      })
    }),
    .groups = "drop"
  ) %>%
  unnest(pairwise) %>%
  group_by(Gene.symbol) %>%
  mutate(BH_adjusted_pairwise = p.adjust(pval, method = "BH")) %>%
  ungroup()

pairwise_results2 %>%
  mutate(pval = round(pval, digits = 5),
         BH_adjusted_pairwise = round(BH_adjusted_pairwise, digits = 5)) 

# ---- Figure 2A ----
tab|>
  dplyr::summarise(n = dplyr::n(), .by = c(filename, Gene.symbol)) |>
  dplyr::filter(n > 1L) 

argmat <- tab %>%
  filter(Element.type == "AMR") %>%
  mutate(value = 1) %>%
  pivot_wider(id_cols = "filename", names_from= "Gene.symbol", values_from = "value")

argmat[is.na(argmat)] <- 0
rownames(argmat) <- argmat$filename

argmat2 <- as.matrix(argmat[2:ncol(argmat)])
rownames(argmat2) <- argmat$filename
armds<- metaMDS(argmat2, distance = "jaccard", autotransform = FALSE, k= 3, trymax = 50)

armds2 <- as.data.frame(scores(armds)$site)
armds2$sra_accession <- rownames(armds2)
armds2 <- merge(armds2, met %>%
                  select(sra_accession, HHS_region, date, source, sample_type_std, ST), by = "sra_accession")

sourcepal <- c("Clinical Bovine" = "#382A54FF",
               "Environmental" = "#3497A9FF",
               "Clinical Human" = '#395D9CFF')


ggplot(armds2 %>%
               mutate(Source = case_when(
                 source == "bovclin" ~ "Clinical Bovine",
                 source == "bovenv" ~ "Environmental",
                 source == "humall" ~ "Clinical Human"
               )), aes(x=NMDS1,y=NMDS2, color = Source)) +
  theme_classic() +
  stat_ellipse(linewidth = 0.75)+
  scale_color_manual(values=sourcepal) + 
  geom_point(size=2.5, alpha = 0.4) +
  theme(axis.title = element_text(size=10, color = "black"),
        axis.text = element_text(size=10, color = "black"),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 10), 
        legend.position = "right")


# ---- Figure 2B and 2C----


# plot arg for preso and copy this code to a different script later 
argbysource <- tab %>%
  mutate(sra_accession = filename) %>%
  select(-c(filename)) %>%
  inner_join(met, by = "sra_accession") %>%
  filter(Element.type == "AMR") %>%
  select(Gene.symbol, Class, sra_accession, source, Element.subtype) %>%
  group_by(Gene.symbol, source) %>%
  mutate(total = n()) %>% # total number of strains containing a given ARG by source
  mutate(ppn = case_when( # ppn of strains containing a give ARG by soruce
    source == "bovclin" ~ total/581,
    source == "bovenv" ~ total/905,
    source == "humall" ~ total/664
  )) %>%
  select(Gene.symbol, Class, source, total, ppn, Element.subtype) %>%
  unique()

CLASSpal <- c("Aminoglycoside" = "#A11A5BFF",
              "Beta-lactam" = "#D2204CFF",
              "Sulfonamide" = "#EE523FFF",
              "Phenicol" = "#F58B63FF",
              "Quinolone" = "#F6BD9AFF",
              "Tetracycline" = "#E98440FF",
              "Trimethoprim" = "#FAEBDDFF")

ggplot(argbysource %>%
                   mutate(Source = case_when(
                     source == "bovclin" ~ "Clinical Bovine",
                     source == "bovenv" ~ "Environmental",
                     source == "humall" ~ "Clinical Human"
                   )) %>%
                   filter(Class != "EFFLUX") %>%
                   filter(Element.subtype != "POINT") %>%
                   mutate(Class = str_to_sentence(Class)), aes(y = Gene.symbol, x = ppn, fill = Class)) +
            scale_x_continuous(expand = c(0,0)) + 
            geom_bar(stat = "identity") + 
            scale_fill_manual(values = CLASSpal) + 
            theme_classic() +
            theme(axis.text.x = element_text(size = 18),
                  axis.text.y = element_text(face = "italic", size = 18),
                  strip.text = element_text(size = 24),
                  axis.title.x = element_text(size = 18),
                  legend.position = "right",
                  legend.text = element_text(size = 18),
                  legend.title = element_text(size = 18)) +
            ylab("")+
            xlab("Proportion Strains with ARG") + 
            facet_wrap(~Source)
          
ggplot(argbysource %>%
         mutate(Source = case_when(
           source == "bovclin" ~ "Clinical Bovine",
           source == "bovenv" ~ "Environmental",
           source == "humall" ~ "Clinical Human"
         )) %>%
         filter(Class != "EFFLUX") %>%
         filter(Element.subtype == "POINT") %>%
         mutate(Class = str_to_sentence(Class)) %>%
         mutate(Class = if_else(str_detect(Class, "uino"), "Quinolone", Class)), aes(y = Gene.symbol, x = ppn, fill = Class)) +
  scale_x_continuous(expand = c(0,0), breaks = c(0, 0.1, 0.2)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("Quinolone" = "#F6BD9AFF",
                               "Multidrug" = "black")) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(face = "italic", size = 10, color = "black"),
        strip.text = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black"),
        legend.position = "right",
        legend.text = element_text(size = 9, color = "black"),
        legend.title = element_text(size = 10, color = "black")) +
  ylab("")+
  xlab("Proportion Strains with ARG") + 
  facet_wrap(~Source)



