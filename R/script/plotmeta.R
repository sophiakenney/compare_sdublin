# Plot metadata

# ---- Load Packages ----
library(tidylog)
library(tidyverse)
library(ggpubr)
library(maps)

setwd("./R")

# ---- Load Tables ----
met <- read.delim("meta/finalset_meta.tsv", sep = "\t")

# ----- Plot Map ---- 
nrow(met %>% filter(HHS_region == "Not Specified")) # 25 without specified HHS region

states_map <- map_data("state")
ggplot(states_map, aes(long, lat, group = region))+
  geom_polygon(fill="lightgray", colour = "white") +
  theme_void()+
  theme(legend.position = "none")

local <- as.data.frame(table(met$location))
hhs <- as.data.frame(table(met$HHS_region))
states_map$region %>%
  unique()


local2 <- local %>%
  filter(Var1 != "Not Specified") %>%
  mutate(region = case_when(
    Var1 == "USA:AZ" ~ "arizona",
    Var1 == "USA:CA" ~ "california",
    Var1 == "USA:CO" ~ " colorado",
    Var1 == "USA:CT" ~ "connecticut",
    Var1 == "USA:DE" ~ "delaware",
    Var1 == "USA:FL" ~ "florida",
    Var1 == "USA:GA" ~ "georgia",
    Var1 == "USA:IA" ~ "iowa",
    Var1 == "USA:ID" ~ "idaho",
    Var1 == "USA:IL" ~ "illinois",
    Var1 == "USA:IN" ~ "indiana",
    Var1 == "USA:KS" ~ "kansas",
    Var1 == "USA:KY" ~ "kentucky",
    Var1 == "USA:MA" ~ "maine",
    Var1 == "USA:MD" ~ "maryland",
    Var1 == "USA:MI" ~ "michigan",
    Var1 == "USA:MN" ~ "minnesota",
    Var1 == "USA:MO" ~ "missouri", 
    Var1 == "USA:MT" ~ "montana",
    Var1 == "USA:NC" ~ "north carolina",
    Var1 == "USA:ND" ~ "north dakota",
    Var1 == "USA:NE" ~ "nebraska" ,
    Var1 == "USA:NJ" ~ "new jersey",
    Var1 == "USA:NM" ~ "new mexico",
    Var1 == "USA:NV" ~ "nevada", 
    Var1 == "USA:NY" ~ "new york",
    Var1 == "USA:OH" ~ "ohio" ,
    Var1 == "USA:OK" ~ "oklahoma",
    Var1 == "USA:OR" ~ "oregon",
    Var1 == "USA:PA" ~ "pennsylvania",
    Var1 == "USA:SC" ~ "south carolina",
    Var1 == "USA:SD" ~ "south dakota",
    Var1 == "USA:TN" ~ "tennessee",
    Var1 == "USA:TX" ~ "texas",
    Var1 == "USA:UT" ~ "utah",
    Var1 == "USA:WA" ~ "washington",
    Var1 == "USA:WI" ~ "wisconsin",
    Var1 == "USA:WY" ~ "wyoming"
  ))

colnames(local2) = c("location", "Freq", "region")

local2<-merge(local2, met %>%
                select(location, HHS_region), by = "location") %>%
  unique()

all_map <- left_join(states_map, local2, by = "region")

all_map2 <- all_map %>% mutate(Freq = case_when(
  HHS_region == "1" ~ 32,
  HHS_region == "2" ~ 154,
  HHS_region == "3" ~ 104,
  HHS_region == "4" ~ 92,
  HHS_region == "5" ~ 563,
  HHS_region == "6" ~ 242,
  HHS_region == "7" ~ 153,
  HHS_region == "8" ~ 123,
  HHS_region == "9" ~ 398,
  HHS_region == "10" ~ 264))

ggplot(all_map2, aes(long, lat, group = group))+
  geom_polygon(aes(fill = Freq), color = "black") + 
  theme_void()+
  scale_fill_gradient(high = "#702058FF", low =  "#EFE6EC", na.value = "white") + 
  theme(legend.position = "right", 
        text = element_text(size = 15),
        plot.margin = margin(10, 40, 10, 10)) + 
  labs(caption = "25 strains not assigned an HHS region")

# --- Plot Date ----
unique(met$date)

ggplot(met %>%
         # mutate(date = ifelse(date == "TBD", "2016", date))%>% # recode to 2016 for now
         mutate(source = case_when(
           source == "bovclin" ~ "Clinical Bovine",
           source == "bovenv" ~ "Environmental",
           source == "humall" ~ "Clinical Human"
         )), aes(x=date, fill=source))+
  geom_bar(color="white") +
  theme_classic() +
  scale_y_continuous(expand = c(0,1)) + 
  scale_x_continuous(breaks = c(2002:2023), expand = c(0,0)) + 
  scale_fill_viridis_d(option = "G", begin = 0.2, end = 0.6) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        #text = element_text(colour = "black", size = 14),
        legend.position = "bottom") +
  ylab("Count")


# ----- Plot pie charts ---- 
sourcehhs <- met %>%
  group_by(HHS_region) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  group_by(HHS_region, source) %>%
  mutate(nsource = n()) %>%
  mutate(ppn = nsource/total) %>%
  ungroup() %>%
  select(HHS_region, total, nsource, ppn, source) %>%
  unique() %>%
  mutate(source = case_when(
    source == "bovclin" ~ "Clinical Bovine",
    source == "bovenv" ~ "Environmental",
    source == "humall" ~ "Clinical Human"
  ))


ggplot(sourcehhs %>%
         filter(!HHS_region %in% c("Not Specified")) %>%
         mutate(HHS_region = factor(HHS_region, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))), aes(x = "", y = ppn, fill = source, width = 1)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_viridis_d(option = "G", begin = 0.2, end = 0.6) + 
  theme_void() + 
  theme(strip.text = element_text(size = 30),
        legend.position = "none") + 
  facet_wrap(~HHS_region, ncol = 2)