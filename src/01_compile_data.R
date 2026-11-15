#!/usr/bin/env Rscript
# =====================================================================
# 01_compile_data.R
# Assemble analysis-ready LAG growth data + a matching topology.
#
# Inputs  (data_raw/)
#   compiled/dinosaur_growth_theropod_femora.txt  -- clean long-format spine
#   compiled/measures.csv                         -- wide circ+mass compilation
#   external/.../tree_theropodGrowth.nex          -- D'Emic 2023 genus tree
#
# Outputs (data_clean/)
#   lag_femur_long.csv        one row / LAG / specimen (femur circumference)
#   measures_long.csv         reshaped wide table (broad taxa; circ + DME mass)
#   theropod_femur_tree.nwk   genus-level topology pruned to the femur taxa
#   data_dictionary.md        column definitions + provenance
#
# Design choice (see docs/dino_growth_hierarchical_plan.md sec 1.2):
# we model the *natively measured* diaphyseal circumference, not a
# back-calculated mass. age is the LAG index; the unknown resorption
# offset (first preserved ring != birth) is absorbed downstream by a
# per-specimen t0.
# =====================================================================

suppressMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ape)
})

root      <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
raw       <- file.path(root, "data_raw")
clean     <- file.path(root, "data_clean")
dir.create(clean, showWarnings = FALSE)

# ---------------------------------------------------------------------
# 1. Femur LAG spine  -------------------------------------------------
# ---------------------------------------------------------------------
femur_path <- file.path(raw, "compiled", "dinosaur_growth_theropod_femora.txt")
femur <- read.table(femur_path, header = TRUE, stringsAsFactors = FALSE)

femur_long <- femur %>%
  transmute(
    clade        = Clade1,
    family       = Clade2,
    species      = Species,
    genus        = str_replace(species, "^cf\\.", "") %>% str_extract("^[^_]+"),
    specimen     = Specimen,
    element      = "femur",
    zone         = Location,                 # PBC = cortex; EFS = adult plateau
    lag_index    = as.integer(LAG),
    age          = as.integer(LAG),          # age proxy = LAG count (see header)
    circ_mm      = as.numeric(Circum_mm),
    is_efs       = Location == "EFS"
  ) %>%
  arrange(species, specimen, lag_index)

write_csv(femur_long, file.path(clean, "lag_femur_long.csv"))

# quick per-specimen summary to the console
spec_summary <- femur_long %>%
  group_by(species, specimen) %>%
  summarise(n_lag = n(),
            has_efs = any(is_efs),
            circ_min = min(circ_mm), circ_max = max(circ_mm),
            .groups = "drop")
cat("\n--- femur LAG spine: per-specimen ---\n")
print(as.data.frame(spec_summary))
cat(sprintf("\nTotal: %d LAG obs | %d specimens | %d species\n",
            nrow(femur_long),
            n_distinct(femur_long$specimen),
            n_distinct(femur_long$species)))

# ---------------------------------------------------------------------
# 2. Reshape the wide measures.csv (broad taxa; circ + DME mass) ------
#    Uses the author's existing python reshaper, unchanged.
# ---------------------------------------------------------------------
measures_in  <- file.path(raw, "compiled", "measures.csv")
measures_out <- file.path(clean, "measures_long.csv")
py <- Sys.which("python3")
if (nzchar(py) && file.exists(measures_in)) {
  reshaper <- file.path(root, "src", "reshape_measures.py")
  # the script hard-codes relative paths; run it from a scratch dir with copies
  tmp <- tempfile("reshape"); dir.create(tmp)
  file.copy(measures_in, file.path(tmp, "measures.csv"), overwrite = TRUE)
  file.copy(reshaper,     file.path(tmp, "reshape_measures.py"), overwrite = TRUE)
  old <- setwd(tmp); on.exit(setwd(old), add = TRUE)
  system2(py, "reshape_measures.py")
  if (file.exists(file.path(tmp, "measures_long.csv")))
    file.copy(file.path(tmp, "measures_long.csv"), measures_out, overwrite = TRUE)
  setwd(old)
  cat(sprintf("\nWrote %s\n", measures_out))
} else {
  cat("\n[skip] python3 or measures.csv not found; measures_long.csv not built\n")
}

# ---------------------------------------------------------------------
# 3. Genus-level topology pruned to the femur taxa --------------------
#    Topology-only (plan sec 1.4 Option A): branch lengths discarded.
# ---------------------------------------------------------------------
nex <- list.files(file.path(raw, "external"), pattern = "tree_theropodGrowth\\.nex$",
                  recursive = TRUE, full.names = TRUE)[1]
tree <- read.nexus(nex)

# map each femur species -> a tip in the D'Emic genus tree
tip_map <- c(
  Tyrannosaurus_rex            = "Tyrannosaurus",
  cf.Tyrannosaurus_rex         = "Tyrannosaurus",
  Gorgosaurus_libratus         = "Gorgosaurus",
  Rativates_evadens            = "Rativates",
  Acrocanthosaurus_atokensis   = "Acrocanthosaurus",
  Allosaurus_fragilis          = "Allosaurus",
  Campanas_Carcharodontosaurid = "Mapusaurus"   # indet. carcharodontosaurid proxy
)
keep <- unique(tip_map[unique(femur_long$species)])
keep <- keep[keep %in% tree$tip.label]
pruned <- keep.tip(tree, keep)
write.tree(pruned, file.path(clean, "theropod_femur_tree.nwk"))

# also persist the species->tip lookup so the model script can build `gr(cov=)`
write_csv(tibble(species = names(tip_map), phylo_tip = unname(tip_map)),
          file.path(clean, "species_to_tip.csv"))

cat("\n--- pruned topology tips ---\n"); print(pruned$tip.label)

# ---------------------------------------------------------------------
# 4. Data dictionary --------------------------------------------------
# ---------------------------------------------------------------------
dd <- c(
"# Data dictionary -- data_clean/",
"",
"## lag_femur_long.csv  (PRIMARY analysis table)",
"One row per LAG (growth ring) per specimen. Femur diaphyseal circumference.",
"Source: data_raw/compiled/dinosaur_growth_theropod_femora.txt (D'Emic et al. 2023-derived).",
"",
"| column | type | meaning |",
"|---|---|---|",
"| clade | chr | higher clade (Coelurosauria / Allosauroidea) |",
"| family | chr | family (Tyrannosauridae, Carcharodontosauridae, ...) |",
"| species | chr | taxon label (matches series; `cf.` retained) |",
"| genus | chr | genus parsed from species |",
"| specimen | chr | museum/specimen number (the individual) |",
"| element | chr | skeletal element (all `femur` here) |",
"| zone | chr | PBC = primary bone cortex; EFS = external fundamental system (growth ~ceased) |",
"| lag_index | int | growth ring number, innermost preserved = 1 |",
"| age | int | age proxy = lag_index (resorption offset absorbed by per-specimen t0) |",
"| circ_mm | num | diaphyseal circumference at that LAG (mm) -- the modeled response |",
"| is_efs | lgl | TRUE if the ring is in the EFS (asymptote is real for these individuals) |",
"",
"## measures_long.csv  (SECONDARY; broader taxa)",
"Reshaped from the wide measures.csv. Columns: species, specimen, bone, lag_index,",
"circumference_mm, est_mass_kg, source, notes. Mixes histological (circumference)",
"and composite (mass-only) blocks; many taxa are single-element / single-specimen.",
"NOTE: measures.csv and alligators.csv are identical to copies inside the D'Emic",
"adc8714 package -- they are that supplement's data, not independent sources.",
"",
"## theropod_femur_tree.nwk",
"Genus-level topology (branch lengths dropped) pruned from D'Emic et al. 2023",
"tree_theropodGrowth.nex to the genera present in lag_femur_long.csv.",
"species_to_tip.csv maps each species label to its tree tip (the Campanas",
"carcharodontosaurid is proxied by Mapusaurus).",
""
)
writeLines(dd, file.path(clean, "data_dictionary.md"))
cat("\nWrote data_dictionary.md\n")
cat("\n=== compile complete ===\n")
