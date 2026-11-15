#!/usr/bin/env Rscript
# Stage the 2026-11 literature intake (data_raw/literature_to_add/_intake/) into
# a candidate analysis table WITHOUT overwriting the current one:
#   data_clean/lag_size_long_next.csv, dino_tree_next.nwk, species_to_tip_next.csv
# Changes relative to lag_size_long.csv:
#  1. Allosaurus: drop the 13 Myh_Allosaurus_* series (Myhrvold 2013 S6). Those
#     are Bybee et al. 2006's LAGs POOLED across bones per element and repeated
#     under alternative age assignments (fc1/fc2 = the same 38 femur LAGs), so
#     they are neither individuals nor independent. Add Bybee's native Table 1
#     per-bone series (17 bones not already held via D'Emic 2023).
#     age = preserved ring index and resorption = "present", as for the other
#     theropods; Bybee's retrocalculated ring count is kept in
#     data_clean/bybee_retrocalc.csv as an external check on the model's offset.
#     Periosteal-surface rows are not rings and are left out (as in the D'Emic series).
#     Bybee gives genus only; rows are labelled Allosaurus_fragilis to pool with
#     the same Cleveland-Lloyd femora already held under that name.
#  2. Diplocynodon hantoniensis (Hoffman et al. 2025): 6 femora, half
#     circumference (um) -> circumference_mm (x2/1000). No ages in the paper:
#     age = growth-mark index, resorption = "present". New tip, sister to Alligator.
suppressMessages({ library(dplyr); library(readr); library(ape) })
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean"); I <- file.path(root,"data_raw/literature_to_add/_intake")
base_cols <- c("clade","family","species","genus","specimen","element","dimension",
               "zone","lag_index","age","size","is_efs","source","resorption","efs_present","maturity")
cur <- read_csv(file.path(clean,"lag_size_long.csv"), show_col_types=FALSE, col_types=cols(lag_index=col_character()))

b <- read_csv(file.path(I,"bybee2006_allosaurus/lag_series.csv"), show_col_types=FALSE)
write_csv(b %>% filter(zone!="periosteal") %>% group_by(specimen, element) %>%
            summarise(n_preserved=n(), retrocalc=min(age - lag_index), .groups="drop"),
          file.path(clean,"bybee_retrocalc.csv"))
bn <- b %>% filter(zone!="periosteal", !grepl("^TRUE", already_held)) %>%
  mutate(species="Allosaurus_fragilis", age=lag_index, lag_index=as.character(lag_index), maturity="subadult") %>%
  select(all_of(base_cols))

h <- read_csv(file.path(I,"hoffman2025_diplocynodon/lag_series.csv"), show_col_types=FALSE) %>%
  mutate(specimen=gsub(" ","_",specimen), size=2*size/1000, dimension="circumference_mm",
         age=lag_index, lag_index=as.character(lag_index), is_efs=FALSE, family="Diplocynodontidae",
         maturity=ifelse(efs_present,"adult","subadult")) %>% select(all_of(base_cols))

nxt <- bind_rows(cur %>% filter(!grepl("^Myh_Allosaurus", specimen)), bn, h)
write_csv(nxt, file.path(clean,"lag_size_long_next.csv"))

tree <- read.tree(file.path(clean,"dino_tree.nwk"))
tree2 <- read.tree(text=sub("\\(Alligator,", "((Alligator,Diplocynodon),", write.tree(tree)))
write.tree(tree2, file.path(clean,"dino_tree_next.nwk"))
tm <- read_csv(file.path(clean,"species_to_tip.csv"), show_col_types=FALSE) %>%
  bind_rows(data.frame(species="Diplocynodon_hantoniensis", phylo_tip="Diplocynodon"))
write_csv(tm, file.path(clean,"species_to_tip_next.csv"))
stopifnot(all(unique(nxt$species) %in% tm$species), all(tm$phylo_tip %in% tree2$tip.label))

cat(sprintf("current: %d rows | %d species | %d specimens\n", nrow(cur), n_distinct(cur$species), n_distinct(cur$specimen)))
cat(sprintf("staged : %d rows | %d species | %d specimens | %d tips\n", nrow(nxt), n_distinct(nxt$species), n_distinct(nxt$specimen), Ntip(tree2)))
print(nxt %>% group_by(dimension) %>% summarise(rows=n(), species=n_distinct(species), specimens=n_distinct(specimen)) %>% as.data.frame())
