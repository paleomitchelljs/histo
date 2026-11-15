#!/usr/bin/env Rscript
# =====================================================================
# 12_expand_dataset.R  -- low-hanging dataset expansion from in-hand data:
#   * Myhrvold 2013 Table S6  (pone.0081917.s018.xlsx): 31 digitized
#     age-vs-size series, ~13 taxa incl. NEW clades (ceratopsian, dromaeosaur,
#     coelophysoid, extra hadrosaur/sauropods) + multi-element Allosaurus.
#   * Alligator (alligators.csv, Lee & O'Connor 2013): first EXTANT outgroup.
#
# Caveats baked into metadata: Myhrvold data are digitized + rescaled + mix
# back-calc methods (dimension = "myhrvold_rescaled"); treated as a linear size
# with per-specimen asymptote, age ~ true age (resorption "absent").
#
# Writes EXPANDED data_clean files (backs up the 13-taxon versions as *_v1).
# =====================================================================
suppressMessages({ library(dplyr); library(readr); library(readxl); library(ape); library(stringr); library(tidyr) })
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean")

base_cols <- c("clade","family","species","genus","specimen","element","dimension",
               "zone","lag_index","age","size","is_efs","source","resorption","efs_present","maturity")

# ---- 1. existing 13-taxon unified table ----------------------------
cur <- read_csv(file.path(clean,"lag_size_long.csv"), show_col_types=FALSE) %>%
  mutate(lag_index=as.character(lag_index))

# ---- 2. Myhrvold S6 -------------------------------------------------
m <- suppressMessages(read_excel(file.path(root,"data_raw/supplements/pone.0081917.s018.xlsx"),
        col_names=c("label","age","size"), skip=1)) %>% as.data.frame()
m$label <- zoo::na.locf(m$label, na.rm=FALSE)
m <- m %>% filter(!is.na(age), !is.na(size)) %>% mutate(age=as.numeric(age), size=as.numeric(size))

# label -> taxonomy + element. genus = first token; element from Allosaurus codes.
elem_of <- function(lab){
  code <- str_extract(lab, "(?<= )[a-z]{2}[0-9]?$")
  if (is.na(code)) return("femur")               # Erickson-style default (femoral circ)
  p <- substr(code,1,2)
  c(fc="femur",fl="femur",hc="humerus",hl="humerus",uc="ulna",ul="ulna",tc="tibia")[p] %>% unname() %||% "femur"
}
`%||%` <- function(a,b) if (is.null(a)||is.na(a)) b else a
mk <- function(lab){
  g <- str_extract(lab, "^[A-Za-z]+")
  info <- dplyr::case_when(
    g=="Tyrannosaurus" ~ "Coelurosauria|Tyrannosauridae|Tyrannosaurus_rex|Tyrannosaurus",
    g=="Gorgosaurus"   ~ "Coelurosauria|Tyrannosauridae|Gorgosaurus_libratus|Gorgosaurus",
    g=="Albertosaurus" ~ "Coelurosauria|Tyrannosauridae|Albertosaurus_sarcophagus|Albertosaurus",
    g=="Saurornitholestes" ~ "Coelurosauria|Dromaeosauridae|Saurornitholestes_langstoni|Saurornitholestes",
    g=="Syntarsus"     ~ "Coelophysoidea|Coelophysidae|Syntarsus_rhodesiensis|Syntarsus",
    g=="Allosaurus"    ~ "Allosauroidea|Allosauridae|Allosaurus_fragilis|Allosaurus",
    g=="Psittacosaurus"& str_detect(lab,"m") ~ "Ceratopsia|Psittacosauridae|Psittacosaurus_mongoliensis|Psittacosaurus",
    g=="Psittacosaurus" ~ "Ceratopsia|Psittacosauridae|Psittacosaurus_lujiatunensis|Psittacosaurus",
    g=="Hypacrosaurus" ~ "Ornithischia|Hadrosauridae|Hypacrosaurus_sp|Hypacrosaurus",
    g=="Apatosaurus"   ~ "Sauropoda|Diplodocidae|Apatosaurus_sp|Apatosaurus",
    g=="Alamosaurus"   ~ "Sauropoda|Titanosauria|Alamosaurus_sanjuanensis|Alamosaurus",
    g=="Janenschia"    ~ "Sauropoda|Turiasauria|Janenschia_sp|Janenschia",
    g=="Northampton"   ~ "Sauropoda|Eusauropoda|Northampton_sauropod|Northampton",
    g=="Massospondylus" ~ "Sauropodomorpha|Massospondylidae|Massospondylus_carinatus|Massospondylus",
    TRUE ~ NA_character_)
  if (is.na(info)) return(NULL)
  parts <- str_split(info,"\\|")[[1]]
  data.frame(clade=parts[1], family=parts[2], species=parts[3], genus=parts[4],
             specimen=paste0("Myh_", gsub("[ .]","_",lab)), element=elem_of(lab))
}
myh <- bind_rows(lapply(unique(m$label), function(lab){
  meta <- mk(lab); if (is.null(meta)) return(NULL)
  rows <- m[m$label==lab,]
  cbind(meta[rep(1,nrow(rows)),], dimension="myhrvold_rescaled", zone=NA_character_,
        lag_index=as.character(seq_len(nrow(rows))), age=rows$age, size=rows$size,
        is_efs=FALSE, source="Myhrvold 2013 (PLOS ONE e0081917, Table S6; digitized/rescaled)",
        resorption="absent", efs_present=FALSE, maturity="subadult")
})) %>% select(all_of(base_cols))

# ---- 3. Alligator (reshape paired-CGM -> indexed series) ------------
al <- suppressMessages(read_csv(file.path(clean,"..","data_raw/compiled/alligators.csv"),
        skip=2, col_names=FALSE, show_col_types=FALSE))
# cols: specimen, smaller CGM, larger CGM (first 3 used)
al <- al[,1:3]; names(al) <- c("spec","small","large")
al <- al %>% filter(!is.na(spec), !is.na(small)) %>%
  mutate(spec=as.character(spec), small=as.numeric(small), large=as.numeric(large))
allig <- bind_rows(lapply(unique(al$spec), function(s){
  r <- al[al$spec==s,]; seq_circ <- unique(round(c(r$small, r$large),2))
  seq_circ <- sort(seq_circ)
  data.frame(clade="Crocodylia", family="Alligatoridae", species="Alligator_mississippiensis",
             genus="Alligator", specimen=paste0("Allig_",s), element="femur",
             dimension="circumference_mm", zone=NA_character_,
             lag_index=as.character(seq_along(seq_circ)), age=seq_along(seq_circ), size=seq_circ,
             is_efs=FALSE, source="Lee & O'Connor 2013 (Alligator; extant outgroup)",
             resorption="present", efs_present=FALSE, maturity="subadult")
})) %>% select(all_of(base_cols))

# ---- 4. combine + write --------------------------------------------
exp <- bind_rows(cur %>% select(all_of(base_cols)), myh, allig)
file.copy(file.path(clean,"lag_size_long.csv"), file.path(clean,"lag_size_long_v1.csv"), overwrite=TRUE)
write_csv(exp, file.path(clean,"lag_size_long.csv"))

# ---- 5. expanded tree (topology-only; approximate consensus) -------
nwk <- paste0("(Alligator,((Psittacosaurus,(Maiasaura,Hypacrosaurus)),",
  "(((Plateosaurus,Massospondylus),(Mamenchisaurus,(Northampton,",
  "((Apatosaurus,(Diplodocus,Janenschia)),(Camarasaurus,Alamosaurus))))),",
  "(Syntarsus,((Allosaurus,(Acrocanthosaurus,Meraxes)),",
  "(Rativates,(Saurornitholestes,(Albertosaurus,(Gorgosaurus,Tyrannosaurus)))))))));")
tree <- read.tree(text=nwk)
file.copy(file.path(clean,"dino_tree.nwk"), file.path(clean,"dino_tree_v1.nwk"), overwrite=TRUE)
write.tree(tree, file.path(clean,"dino_tree.nwk"))

tip_map <- exp %>% distinct(species, genus) %>%
  mutate(phylo_tip = case_when(
    species=="cf.Tyrannosaurus_rex" ~ "Tyrannosaurus",
    species=="Campanas_Carcharodontosaurid" ~ "Meraxes",
    species=="Mamenchisauridae_indet" ~ "Mamenchisaurus",
    species=="Diplodocidae_indet" ~ "Diplodocus",
    TRUE ~ genus))
file.copy(file.path(clean,"species_to_tip.csv"), file.path(clean,"species_to_tip_v1.csv"), overwrite=TRUE)
write_csv(tip_map %>% select(species, phylo_tip), file.path(clean,"species_to_tip.csv"))

# ---- report --------------------------------------------------------
stopifnot(all(tip_map$phylo_tip %in% tree$tip.label))
cat(sprintf("EXPANDED: %d rows | %d species | %d specimens | %d clades | %d tree tips\n",
    nrow(exp), n_distinct(exp$species), n_distinct(exp$specimen), n_distinct(exp$clade),
    length(tree$tip.label)))
cat("dimensions:", paste(unique(exp$dimension), collapse=", "), "\n")
cat("\nNEW species added:\n")
print(setdiff(unique(exp$species), unique(cur$species)))
cat("\ntips missing data:", paste(setdiff(tree$tip.label, tip_map$phylo_tip), collapse=", "), "\n")
cat("\nper-clade species/specimens:\n")
print(exp %>% group_by(clade) %>% summarise(species=n_distinct(species), specimens=n_distinct(specimen), points=n()) %>% as.data.frame())
