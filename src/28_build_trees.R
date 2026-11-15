#!/usr/bin/env Rscript
# Topologies for the staged data and for the wider intake taxa. Hand-assembled
# genus-level topologies, no branch lengths. Writes:
#   data_clean/dino_tree_next.nwk       model tree for lag_size_long_next.csv (23 tips)
#   data_clean/amniote_tree_draft.nwk   all taxa with ANY intake data (ring series,
#                                       age-class means, feature tables), for the
#                                       feature matrix / later threshold characters
#   data_clean/species_to_tip_draft.csv species label -> draft-tree tip
# Changes to the dinosaur part relative to dino_tree.nwk:
#   * Diplocynodon added as sister to Alligator.
#   * Saurornitholestes moved from sister-of-Tyrannosauridae to sister of Rativates
#     (Maniraptoriformes), as in D'Emic et al. 2023's tree_theropodGrowth.nex.
#   * Janenschia placed in Turiasauria: outside Neosauropoda, as sister to it
#     (was sister to Diplodocus). JM 2026-11-13.
#   * Tyrannosauridae: ((Albertosaurus,Gorgosaurus),Tyrannosaurus). JM 2026-11-13.
suppressMessages({ library(ape); library(readr); library(dplyr) })
root <- "/Users/jmitchell/Library/CloudStorage/Dropbox/Research/00_current_Histology"; clean <- file.path(root,"data_clean")

orn   <- "(Psittacosaurus,(Maiasaura,Hypacrosaurus))"
sauro <- "((Plateosaurus,Massospondylus),(Mamenchisaurus,(Northampton,(Janenschia,((Apatosaurus,Diplodocus),(Camarasaurus,Alamosaurus))))))"
allo  <- "(Allosaurus,(Acrocanthosaurus,Meraxes))"
tyr   <- "((Albertosaurus,Gorgosaurus),Tyrannosaurus)"
thero <- function(coelo, mani) sprintf("(%s,(%s,(%s,%s)))", coelo, allo, tyr, mani)
dino  <- function(o, t) sprintf("(%s,(%s,%s))", o, sauro, t)

model_tree <- sprintf("((Alligator,Diplocynodon),%s);", dino(orn, thero("Syntarsus", "(Rativates,Saurornitholestes)")))
write.tree(read.tree(text=model_tree), file.path(clean,"dino_tree_next.nwk"))

# ---- wider draft tree -------------------------------------------------
birds <- "((Struthio,Dromaius),(Anas,(Streptopelia,(Turdus,Aptenodytes))))"
mani  <- sprintf("(Rativates,((Troodon,Saurornitholestes),%s))", birds)
dinoW <- dino(sprintf("(Lesothosaurus,%s)", orn), thero("(Syntarsus,Coelophysis)", mani))
# Pseudosuchia: phytosaur (Rutiodon) placed as the basal pseudosuchian; some
# analyses put phytosaurs outside Archosauria instead.
pseudo <- "(Rutiodon,((Typothorax,Calyptosuchus),(Postosuchus,(Crocodylus,(Alligator,Diplocynodon)))))"
arch   <- sprintf("(%s,%s)", pseudo, dinoW)
turt   <- "(Chelodina,(Pelodiscus,Trachemys))"
lacert <- "(Psammodromus,(Ophisops,(Lacerta,Podarcis,Zootoca)))"          # Lacertini unresolved
iguania<- "((Iguana,Phymaturus),(Agama,Phrynocephalus))"
squam  <- sprintf("(%s,(Varanus,%s))", lacert, iguania)
mamm   <- "(Loxodonta,((Microcebus,Macaca),(Cavia,Mus)))"
draft  <- sprintf("(%s,(%s,(%s,%s)));", mamm, squam, turt, arch)
dt <- read.tree(text=draft); write.tree(dt, file.path(clean,"amniote_tree_draft.nwk"))

tm <- read_csv(file.path(clean,"species_to_tip_next.csv"), show_col_types=FALSE)
extra <- tribble(~species, ~phylo_tip,
  "Diplocynodon_hantoniensis","Diplocynodon", "Crocodylus_niloticus","Crocodylus",
  "Struthio_camelus","Struthio", "Dromaius_novaehollandiae","Dromaius", "Anas_platyrhynchos","Anas",
  "Streptopelia_decaocto","Streptopelia", "Turdus_merula","Turdus", "Aptenodytes_patagonicus","Aptenodytes",
  "Postosuchus","Postosuchus", "Typothorax","Typothorax", "Calyptosuchus","Calyptosuchus", "Rutiodon","Rutiodon",
  "Lesothosaurus","Lesothosaurus", "Coelophysis","Coelophysis", "Troodon","Troodon",
  "Trachemys_scripta","Trachemys", "Pelodiscus_sinensis","Pelodiscus", "Chelodina_siebenrocki","Chelodina",
  "Lacerta_vivipara","Zootoca", "Podarcis_muralis","Podarcis", "Lacerta_agilis_boemica","Lacerta", "Lacerta_strigata","Lacerta",
  "Psammodromus_algirus","Psammodromus", "Ophisops_elegans","Ophisops",
  "Varanus_exanthematicus","Varanus", "Varanus_niloticus","Varanus",
  "Iguana_iguana","Iguana", "Phymaturus_aguanegra","Phymaturus", "Phymaturus_zapalensis","Phymaturus",
  "Agama_impalearis","Agama", "Phrynocephalus_vlangalii","Phrynocephalus",
  "Loxodonta_africana","Loxodonta", "Microcebus_murinus","Microcebus", "Macaca_arctoides","Macaca",
  "Cavia_porcellus","Cavia", "Mus_musculus","Mus")
tmd <- bind_rows(tm, extra) %>% distinct(species, .keep_all=TRUE)
write_csv(tmd, file.path(clean,"species_to_tip_draft.csv"))
mt <- read.tree(file.path(clean,"dino_tree_next.nwk"))
stopifnot(all(tm$phylo_tip %in% mt$tip.label), all(mt$tip.label %in% tm$phylo_tip),
          all(tmd$phylo_tip %in% dt$tip.label), all(dt$tip.label %in% tmd$phylo_tip), !anyDuplicated(dt$tip.label))
cat(sprintf("model tree: %d tips | draft tree: %d tips, %d species mapped\n", Ntip(mt), Ntip(dt), nrow(tmd)))
