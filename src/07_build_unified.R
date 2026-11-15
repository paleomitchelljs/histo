#!/usr/bin/env Rscript
# =====================================================================
# 07_build_unified.R
# Merge all clades into one long table + categorical histology metadata,
# and build a genus-level cross-Dinosauria topology.
#
# Sources:
#   data_clean/lag_femur_long.csv            theropod femur circumference (mm)
#   data_clean/maiasaura_lag_long.csv        Maiasaura tibia circumference (mm)
#   data_raw/external/griebeler/...csv       sauropodomorph % final bone length
#
# Output:
#   data_clean/lag_size_long.csv   unified; response = `size`, with `dimension`
#   data_clean/dino_tree.nwk       12-tip genus topology
#   data_clean/species_to_tip.csv  species -> tree tip
#
# Categorical histology metadata (per specimen) -- for the offset model:
#   resorption  : absent / present / unknown  (drives missing-LAG offset prior)
#   efs_present : did growth cease (outer rings = EFS / near asymptote)?
#   maturity    : immature / subadult / adult
# =====================================================================
suppressMessages({ library(dplyr); library(readr); library(ape); library(stringr) })
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root, "data_clean")

base_cols <- c("clade","family","species","genus","specimen","element","dimension",
               "zone","lag_index","age","size","is_efs","source")

ther <- read_csv(file.path(clean,"lag_femur_long.csv"), show_col_types=FALSE) %>%
  transmute(clade, family, species, genus, specimen, element,
            dimension="circumference_mm", zone, lag_index=as.character(lag_index),
            age, size=circ_mm, is_efs, source="D'Emic et al. 2023 (Science)")

maia <- read_csv(file.path(clean,"maiasaura_lag_long.csv"), show_col_types=FALSE) %>%
  transmute(clade, family, species, genus, specimen, element,
            dimension="circumference_mm", zone, lag_index=as.character(lag_index),
            age, size=circ_mm, is_efs, source)

saur <- read_csv(file.path(root,"data_raw/external/griebeler/griebeler_sauropod_long.csv"),
                 show_col_types=FALSE) %>%
  transmute(clade, family, species, genus, specimen, element,
            dimension="pct_bone_length", zone=NA_character_,
            lag_index=as.character(annulus), age, size=pct_bone_length,
            is_efs=FALSE, source)

dat <- bind_rows(ther, maia, saur)

# ---- per-specimen categorical histology metadata -------------------
# resorption: theropods & sauropods preserve only outer rings (inner resorbed);
#   Maiasaura series are anchored at the perinate -> no missing inner record.
sp_meta <- dat %>% group_by(specimen) %>%
  summarise(clade1=first(clade), dimension=first(dimension),
            any_efs=any(is_efs), max_size=max(size), .groups="drop")
sp_meta <- sp_meta %>%
  group_by(species_dim = dimension) %>% mutate() %>% ungroup()

# species-level max for an EFS heuristic on Maiasaura/sauropods
sp_max <- dat %>% group_by(species) %>% summarise(sp_max=max(size), .groups="drop")
specmax <- dat %>% group_by(specimen, species) %>% summarise(sm=max(size), .groups="drop") %>%
  left_join(sp_max, by="species") %>% mutate(near_asym = sm >= 0.95*sp_max)

meta <- dat %>% distinct(specimen, species, clade) %>%
  left_join(dat %>% group_by(specimen) %>% summarise(any_efs=any(is_efs), .groups="drop"), by="specimen") %>%
  left_join(specmax %>% select(specimen, near_asym), by="specimen") %>%
  mutate(
    resorption = case_when(
      str_detect(species, "Maiasaura") ~ "absent",      # perinate-anchored
      clade %in% c("Sauropodomorpha","Sauropoda","Macronaria") ~ "present",
      TRUE ~ "present"),                                  # theropods: truncated
    efs_present = any_efs | (str_detect(species,"Maiasaura") & near_asym) |
                  (clade %in% c("Sauropodomorpha","Sauropoda","Macronaria") & near_asym),
    maturity = if_else(efs_present, "adult", "subadult")) %>%
  select(specimen, resorption, efs_present, maturity)

dat <- dat %>% left_join(meta, by="specimen")
write_csv(dat, file.path(clean,"lag_size_long.csv"))

# ---- genus-level cross-Dinosauria topology ------------------------
nwk <- paste0("(Maiasaura,((Plateosaurus,(Mamenchisaurus,((Apatosaurus,Diplodocus),",
              "Camarasaurus))),((Allosaurus,(Acrocanthosaurus,Meraxes)),",
              "(Rativates,(Gorgosaurus,Tyrannosaurus)))));")
tree <- read.tree(text=nwk)
write.tree(tree, file.path(clean,"dino_tree.nwk"))

tip_map <- tribble(~species, ~phylo_tip,
  "Tyrannosaurus_rex","Tyrannosaurus", "cf.Tyrannosaurus_rex","Tyrannosaurus",
  "Gorgosaurus_libratus","Gorgosaurus", "Rativates_evadens","Rativates",
  "Campanas_Carcharodontosaurid","Meraxes", "Acrocanthosaurus_atokensis","Acrocanthosaurus",
  "Allosaurus_fragilis","Allosaurus", "Maiasaura_peeblesorum","Maiasaura",
  "Plateosaurus_engelhardti","Plateosaurus", "Mamenchisauridae_indet","Mamenchisaurus",
  "Apatosaurus_sp","Apatosaurus", "Diplodocidae_indet","Diplodocus",
  "Camarasaurus_sp","Camarasaurus")
write_csv(tip_map, file.path(clean,"species_to_tip.csv"))

# ---- report -------------------------------------------------------
cat(sprintf("\nUNIFIED: %d rows | %d species | %d specimens | %d clades\n",
    nrow(dat), n_distinct(dat$species), n_distinct(dat$specimen), n_distinct(dat$clade)))
cat("dimensions:", paste(unique(dat$dimension), collapse=", "), "\n\n")
print(dat %>% group_by(clade, species) %>%
        summarise(specimens=n_distinct(specimen), points=n(),
                  age=sprintf("%g-%g", min(age), max(age)),
                  .groups="drop") %>% as.data.frame(), row.names=FALSE)
cat("\n-- specimen metadata tally --\n")
print(dat %>% distinct(specimen, resorption, efs_present) %>% count(resorption, efs_present))
stopifnot(all(tip_map$phylo_tip %in% tree$tip.label))
cat("\nall species map to tree tips:", all(unique(dat$species) %in% tip_map$species), "\n")
cat("tree tips:", paste(tree$tip.label, collapse=", "), "\n")
