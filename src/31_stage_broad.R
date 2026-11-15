#!/usr/bin/env Rscript
# Stage the BROAD dataset for model 32 (individual age at death + element level +
# count-and-size module + three-state EFS). Nothing existing is overwritten. Writes
#   data_clean/broad_rings.csv        one row per sized growth mark (ring-series elements)
#   data_clean/broad_elements.csv     one row per element (ring-series AND count-and-size)
#   data_clean/broad_individuals.csv  one row per individual (age prior type)
#   data_clean/broad_tree.nwk, broad_species_to_tip.csv
#   data_clean/broad_excluded.csv     what was left out, and why
#
# OBSERVATION TYPES
#  ring  an element with a series of sized growth marks. age0 = mark index (or the
#        source's assigned age); the model adds the element's lost inner marks.
#  cs    count-and-size: an element with a mark COUNT and/or a SIZE AT DEATH. One
#        (age, size) point with uncertain age. Also used for a body-size measure of an
#        individual (skull length, SVL), for age-class means (nrep = n animals) and for
#        population cross-sections (one pseudo-individual per point).
# INDIVIDUAL AGE
#  age_est NA   age = (marks counted) + (marks lost), lost >= 0             [count-based]
#  age_est set  age ~ normal(age_est, age_sd): known ages, or the source's own age estimate
#  t_add        years added to the individual's age for a size-at-death observation:
#               0.5 when the age is a mark count (death falls between marks), 0 when the
#               age is the age at death itself.
# EFS (three states)  1 present | 0 absent (source says none) | 2 unknown (not stated,
#               or outer cortex eroded). eroded = 1 additionally lets OUTER marks be lost.
#
# JUDGEMENT CALLS MADE HERE (all flagged in histoPlan.md for JM):
#  * Myhrvold 2013 S6 series are split by `myh_class` below: cross-sections become
#    cs points with age ~ normal(assigned age, max(0.5, 0.15 age)); consecutive-age
#    monotone series stay ring series. Psittacosaurus: version l1 only; Tyrannosaurus: _1 only.
#  * Count ranges ("5 to 6") enter as the midpoint; "at least n" enters as n.
#  * Lystrosaurus size = % of maximum skull length where printed, otherwise the
#    midpoint of the authors' size class with size_logsd set from the class width.
#  * EFS = absent only where the source says so; otherwise unknown.
#  * Species with fewer than MIN_OBS sized observations are left out.
suppressMessages({ library(dplyr); library(readr); library(ape); library(tibble) })
# Usage: Rscript 31_stage_broad.R [MIN_OBS=2] [prefix=broad] [nomyh]
#   e.g. `Rscript 31_stage_broad.R 9 core` writes data_clean/core_*.csv: a narrower test set
#   without the weakest species (JM 2026-11-14: narrowing the data for testing is fine).
#   A third argument `nomyh` leaves out every Myhrvold 2013 S6 series (ring and cross-section);
#   added 2026-11-14 for model 35 (docs/model_revision_plan.md, data handling). Default: kept.
# Added 2026-11-14 (used by model 35, ignored by model 32):
#   elements: count_lo, count_hi  interval for the growth-mark count (hi = Inf for "at least");
#             parsed from the printed count text; equal to `count` where the count is exact.
#   individuals: sex  "female" / "male" / NA, where the source gives it.
a <- commandArgs(trailingOnly=TRUE)
MIN_OBS <- ifelse(is.na(a[1]), 2, as.numeric(a[1])); PRE <- ifelse(is.na(a[2]), "broad", a[2])
DROP_MYH <- !is.na(a[3]) && a[3] == "nomyh"
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean"); I <- file.path(root,"data_raw/literature_to_add/_intake")
rd <- function(...) read_csv(file.path(I, ...), show_col_types=FALSE, col_types=cols(.default="c"))
num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.]", "", x)))
# count parser: all integers in the string -> midpoint of their range
pcount <- function(x) vapply(x, function(s){ if (is.na(s)) return(NA_real_)
  v <- as.numeric(regmatches(s, gregexpr("[0-9]+\\.?[0-9]*", s))[[1]]); if (!length(v)) NA_real_ else mean(range(v)) }, 0)
eltype_of <- function(e){ e <- tolower(e)
  ifelse(grepl("femur",e),"femur", ifelse(grepl("humer",e),"humerus", ifelse(grepl("ulna",e),"ulna",
  ifelse(grepl("radius",e),"radius", ifelse(grepl("tibia",e),"tibia", ifelse(grepl("fibula",e),"fibula","other")))))) }

rings <- list(); els <- list(); inds <- list(); excl <- list()
add_ring <- function(d) rings[[length(rings)+1]] <<- d   # species individual element dimension age0 size source
add_el <- function(d){                                   # one row per element
  def <- list(kind="cs", size=NA_real_, count=NA_real_, count_min=FALSE, efs=2L, resorb=1L, eroded=0L,
              nrep=1, t_add=0.5, size_logsd=0, note=NA_character_)
  for (n in names(def)) if (!n %in% names(d)) d[[n]] <- def[[n]]
  els[[length(els)+1]] <<- d }
add_ind <- function(d){ if (!"age_est" %in% names(d)) { d$age_est <- NA_real_; d$age_sd <- NA_real_ }
  inds[[length(inds)+1]] <<- d }
skip <- function(source, what, why) excl[[length(excl)+1]] <<- tibble(source=source, what=what, why=why)

# ---- reader for the harmonised `model_inputs.csv` (one row per element x size measure) ----
# elem: named list  size_measure regex -> list(dim=, f=transform)   size of the sectioned element
# body: named list  size_measure       -> list(dim=, el=)           one body-size element per animal
# Counts: ranges enter as the midpoint, "at least" as given (count_min). EFS as coded by the
# agent (absent only where the source says so). Individuals are keyed on specimen number, plus
# the species where one number holds two species.
mi <- function(folder, src, elem=list(), body=list(), keep=function(d) d){
  d <- rd(folder, "model_inputs.csv") %>% mutate(species=gsub("[ .]+", "_", sub("\\.$", "", species))) %>% keep()
  d <- d %>% group_by(individual) %>% mutate(nsp=n_distinct(species)) %>% ungroup() %>% mutate(id=ifelse(nsp > 1, paste0(individual, " (", species, ")"), individual)) %>%
    mutate(cnt=pcount(count), cmin=!is.na(count_is_minimum) & count_is_minimum=="TRUE", val=num(size_value),
           efs=unname(c(present=1L, absent=0L, unknown=2L)[efs_state]), ero=as.integer(!is.na(outer_cortex_eroded) & outer_cortex_eroded=="TRUE"))
  first <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else NA }
  for (m in names(elem)) { e <- elem[[m]]; f <- if (is.null(e$f)) identity else e$f
    g <- d %>% group_by(species, id, element) %>%
      summarise(size=first(f(val[!is.na(size_measure) & grepl(m, size_measure)])), count=first(cnt), count_min=any(cmin),
                efs=ifelse(any(efs==1, na.rm=TRUE), 1L, ifelse(any(efs==0, na.rm=TRUE), 0L, 2L)), eroded=max(ero), raw=first(count), .groups="drop") %>%
      filter(!is.na(size) | !is.na(count))
    add_el(g %>% transmute(species, individual=id, element, dimension=e$dim, source=src, size, count, count_min, efs, eroded,
                           note=ifelse(is.na(raw), NA, paste0("count as printed: ", raw)))) }
  for (m in names(body)) { b <- body[[m]]
    g <- d %>% filter(size_measure==m) %>% group_by(species, id) %>% summarise(size=first(val), .groups="drop") %>% filter(!is.na(size))
    add_el(g %>% transmute(species, individual=id, element=b$el, dimension=b$dim, source=src, size, note=m)) }
  add_ind(d %>% distinct(species, individual=id)); invisible(d) }

# ============ 1. existing staged ring table ============================================
cur <- read_csv(file.path(clean,"lag_size_long_next.csv"), show_col_types=FALSE) %>% filter(!is.na(size), !is.na(age))
myh_class <- tribble(~specimen, ~class, ~basis,
  "Myh_Alamosaurus","ring","consecutive ages, monotone sizes; single-bone record (unverified against Myhrvold S6)",
  "Myh_Apatosaurus","ring","consecutive ages, monotone sizes (unverified)",
  "Myh_Janenschia","ring","consecutive ages, monotone sizes (unverified)",
  "Myh_Northampton","ring","consecutive ages, monotone sizes (unverified)",
  "Myh_Hypacrosaurus_fc","ring","consecutive ages, monotone sizes (unverified; fc and tc may be one animal, kept separate)",
  "Myh_Hypacrosaurus_tc","ring","consecutive ages, monotone sizes (unverified)",
  "Myh_Albertosaurus","cross","non-monotone sizes; one point per individual",
  "Myh_Gorgosaurus","cross","repeated ages; one point per individual",
  "Myh_Massospondylus","cross","repeated ages, non-monotone sizes",
  "Myh_Psittacosaurus_l1","cross","repeated ages; population sample",
  "Myh_Psittacosaurus_m1","cross","treated as population sample (ages consecutive; unverified)",
  "Myh_Saurornitholestes","cross","non-consecutive ages (unverified)",
  "Myh_Syntarsus","cross","repeated ages",
  "Myh_Tyrannosaurus_1","cross","non-monotone sizes; one point per individual. Version 2 adds two femora (78.8, 128) that match 'Jane' and a MOR specimen; Jane is held from Woodward 2020, so version 1 is used")
# 2026-11-14 intake of the Myhrvold supplements (S5-S7, Text S7; _intake/myhrvold2013_growth_rates/):
# the files do not say what a row is for any series, and the main paper is not held, so the
# classes above remain PATTERN-BASED. Confirmed: l2 = l1 without duplicates, l3 a subset, l4 = l3
# without duplicates; the values are native measurements (Allosaurus = Bybee circumferences in cm),
# not rescaled to an asymptote, so the dimension label `myhrvold_rescaled` is a misnomer kept only
# for continuity; `element = femur` is a default of src/12, not a source statement.
if (DROP_MYH) { skip("Myhrvold 2013 S6", "all series", "left out (nomyh): what a row is remains unknown; transcription errors against Bybee; S6/S7 ages disagree")
  cur <- cur %>% filter(!grepl("^Myh_", specimen)) }
myh <- cur %>% filter(grepl("^Myh_", specimen))
for (s in setdiff(unique(myh$specimen), myh_class$specimen))
  skip("Myhrvold 2013 S6", s, "alternative version of a pooled series; one version kept (l1 / _1)")
cr <- cur %>% filter(!specimen %in% myh_class$specimen[myh_class$class=="cross"],
                     !(grepl("^Myh_", specimen) & !specimen %in% myh_class$specimen))
# EFS: present where coded TRUE; absent only where the source codes EFS explicitly
# (D'Emic 2023, Bybee 2006, Woodward 2015, Hoffman 2025); otherwise unknown.
explicit <- "D'Emic|Bybee|Woodward|Hoffman"
ce <- cr %>% group_by(species, individual=specimen, element, dimension, source) %>%
  summarise(efs=ifelse(any(efs_present), 1L, ifelse(grepl(explicit, source[1]), 0L, 2L)),
            resorb=as.integer(resorption[1] != "absent"), .groups="drop")
add_ring(cr %>% transmute(species, individual=specimen, element, dimension, age0=age, size, source))
add_el(ce %>% mutate(kind="ring"))
add_ind(ce %>% distinct(species, individual))
# Myhrvold cross-sections: one pseudo-individual per point; identical points collapsed (nrep)
mx <- myh %>% filter(specimen %in% myh_class$specimen[myh_class$class=="cross"]) %>%
  group_by(species, specimen, element, dimension, source, age, size) %>% summarise(nrep=n(), .groups="drop") %>%
  group_by(specimen) %>% mutate(individual=paste0(specimen, "_p", row_number())) %>% ungroup()
if (nrow(mx)) add_el(mx %>% transmute(species, individual, element, dimension, source, size, nrep, t_add=0, resorb=0L,
                        note="population cross-section point (Myhrvold 2013 S6)"))
if (nrow(mx)) add_ind(mx %>% transmute(species, individual, age_est=age, age_sd=pmax(0.5, 0.15*age)))

# ============ 2. new ring series =======================================================
## Ceratosaurus (Sombathy et al. 2025): femur BYU 881-12893 (8 marks, all in the EFS), tibia MWC 1.
so <- rd("sombathy2025_ceratosaurus/lag_series.csv") %>% filter(dimension=="circumference_mm", zone %in% c("EFS","LAG")) %>%
  transmute(species, individual=specimen, element, dimension, age0=as.numeric(lag_index), size=as.numeric(size), source="Sombathy et al. 2025 (J Anat)")
add_ring(so); add_el(so %>% distinct(species, individual, element, dimension, source) %>% mutate(kind="ring", efs=1L))
add_ind(so %>% distinct(species, individual))
add_el(tibble(species="Ceratosaurus_sp", individual="BYU 725-5133", element="femur", dimension="circumference_mm",
              source="Sombathy et al. 2025 (J Anat)", size=217, count=0, efs=0L, note="no LAGs visible; Table 2 circumference"))
add_ind(tibble(species="Ceratosaurus_sp", individual="BYU 725-5133"))
skip("Sombathy 2025","rib and osteoderm counts; unsectioned femora and tibiae","no homologous size with a count / no count")

## Orthosuchus (Weiss et al. 2025): 6 limb elements of SAM-PK-K409 + femur of BP/1/4242.
## EFS per element from individuals.csv notes: humerus and ulna outer cortex missing
## (unknown, eroded); radius "No OCL discernible" (absent); K409 femur "possible OCL"
## (unknown); tibia, fibula, BP/1/4242 femur OCL (present).
we <- rd("weiss2025_orthosuchus/lag_series.csv") %>%
  transmute(species="Orthosuchus_stormbergi", individual=specimen, element, dimension, age0=as.numeric(lag_index),
            size=as.numeric(size), source="Weiss et al. 2025 (J Anat)")
wst <- tribble(~individual, ~element, ~efs, ~eroded,
  "SAM-PK-K409","humerus",2L,1L, "SAM-PK-K409","radius",0L,0L, "SAM-PK-K409","ulna",2L,1L, "SAM-PK-K409","femur",2L,0L,
  "SAM-PK-K409","tibia",1L,0L, "SAM-PK-K409","fibula",1L,0L, "BP/1/4242","femur",1L,0L)
add_ring(we); add_el(we %>% distinct(species, individual, element, dimension, source) %>% left_join(wst, by=c("individual","element")) %>% mutate(kind="ring"))
add_el(tibble(species="Orthosuchus_stormbergi", individual=c("SAM-PK-K409","BP/1/4242"), element="femur_length", dimension="bone_length_mm",
              source="Weiss et al. 2025 (J Anat)", size=c(61.6, 63.8), note="femur length at death"))
add_ind(tibble(species="Orthosuchus_stormbergi", individual=c("SAM-PK-K409","BP/1/4242")))
skip("Weiss 2025","SAM-PK-K409 radiale (1 mark) and dorsal rib (5-7)","no size series; rib remodelled")

## Moschorhinus (Huttenlocker & Botha-Brink 2013): sized marks per element; skull length per animal.
mo <- rd("huttenlocker2013_moschorhinus/lag_series.csv") %>%
  mutate(individual=sub("[ab]$", "", specimen), nmark=ave(as.numeric(lag_index), specimen, element, FUN=max))
mos <- "Huttenlocker & Botha-Brink 2013 (Paleobiology)"
mr <- mo %>% filter(!is.na(size)) %>% transmute(species="Moschorhinus_kitchingi", individual, element, dimension, age0=as.numeric(lag_index), size=as.numeric(size), source=mos)
add_ring(mr)
# count = all marks the element preserves (sized or not), so unsized outer marks still count
add_el(mo %>% filter(!is.na(size)) %>% distinct(individual, element, dimension, nmark) %>%
         transmute(species="Moschorhinus_kitchingi", individual, element, dimension, source=mos, count=nmark,
                   kind="ring", efs=ifelse(individual=="NMQR 3939", 0L, 2L)))
# elements with a count but no sized mark, in animals that have other data
mc <- mo %>% group_by(individual, element) %>% summarise(n_sized=sum(!is.na(size)), count=max(as.numeric(lag_index)), .groups="drop") %>%
  filter(n_sized==0, individual %in% c(mr$individual, "SAM-PK-K118"))
add_el(mc %>% transmute(species="Moschorhinus_kitchingi", individual, element, dimension="circumference_mm", source=mos, count, note="count only"))
add_el(tibble(species="Moschorhinus_kitchingi", individual="SAM-PK-K118", element=c("humerus","radius"), dimension="circumference_mm",
              source=mos, count=0, note="count only (0 marks, Table 2)"))
bsl <- tribble(~individual, ~size, "NMQR 3351",240, "NMQR 3939",170, "BP/1/4227",200, "SAM-PK-K118",157, "SAM-PK-K9953",155, "NMQR 48",163)
add_el(bsl %>% mutate(species="Moschorhinus_kitchingi", element="skull", dimension="bsl_mm", source=mos, note="basal skull length (Table 1)"))
add_ind(tibble(species="Moschorhinus_kitchingi", individual=union(mr$individual, "SAM-PK-K118")))
skip("Huttenlocker 2013","NMQR 1640 (femur 6, tibia 4), NMQR 3684 (femur 4)","counts with no measured size (size only as a regression estimate)")

## Crocodylus niloticus (Chinsamy & Pereyra 2026): 4 captive animals of known age, total
## length at three injection dates. Hatching "March-April 2011" taken as 2011-04-01.
## The measurement at death was taken without skin and is left out.
cn <- rd("chinsamy2026_crocodylus_niloticus/body_size_series.csv") %>% filter(!grepl("death", event)) %>%
  transmute(species="Crocodylus_niloticus", individual=paste("Cnil", specimen), element="body", dimension="total_length_mm",
            age0=round(as.numeric(as.Date(date) - as.Date("2011-04-01"))/365.25, 2), size=10*as.numeric(total_length_cm),
            source="Chinsamy & Pereyra 2026")
add_ring(cn); add_el(cn %>% distinct(species, individual, element, dimension, source) %>% mutate(kind="ring", resorb=0L, efs=0L, note="known-age length series"))
add_ind(cn %>% group_by(species, individual) %>% summarise(age_est=max(age0), age_sd=0.1, .groups="drop"))
skip("Chinsamy & Pereyra 2026 / Buso 2026","growth-mark counts of known-age Nile crocodiles","counts exceed years (non-annual marks); held for the marks-per-year calibration (R7)")

# ============ 3. count-and-size sources ================================================
## Tenontosaurus (Werning 2012): major-axis diameter and LAG count per element.
w <- rd("werning2012_tenontosaurus/individuals.csv") %>% filter(element != "rib")
ws <- "Werning 2012 (PLOS ONE)"
add_el(w %>% transmute(species="Tenontosaurus_tilletti", individual=specimen, element, dimension="bone_diameter_mm", source=ws,
        size=num(ifelse(MajA_mm=="?", NA, MajA_mm)), count=pcount(n_LAGs_preserved_as_printed),
        count_min=grepl("at least", n_LAGs_preserved_as_printed), efs=ifelse(efs_present=="TRUE", 1L, 0L),
        note=paste0("LAGs as printed: ", n_LAGs_preserved_as_printed)))
add_ind(w %>% distinct(individual=specimen) %>% mutate(species="Tenontosaurus_tilletti"))

## Jeholosaurus (Han et al. 2020 + Supplemental Tables S1-S2): radius from the section centroid to
## each LAG, 4 animals x 2 bones -> ring series (radius_mm; tibia and fibula use different
## reference points, which the element effect absorbs). The radius at the periphery is not a
## mark: it is used only for the animal with no LAGs (count 0). Element lengths enter as sizes at death.
hs <- "Han et al. 2020 (JVP) + Suppl. Tables S1-S2"
hl <- rd("han2020_jeholosaurus/lag_series.csv") %>% filter(zone=="LAG") %>%
  transmute(species="Jeholosaurus_shangyuanensis", individual=specimen, element=sub("right ", "", element), dimension="radius_mm",
            age0=as.numeric(lag_index), size=as.numeric(size), source=hs)
hm <- rd("han2020_jeholosaurus/model_inputs.csv") %>% mutate(element=sub("right ", "", element), val=num(size_value))
he <- hm %>% group_by(individual, element) %>% summarise(efs=ifelse(any(efs_state=="absent"), 0L, 2L), count=pcount(count[1]), .groups="drop")
add_ring(hl); add_el(hl %>% distinct(species, individual, element, dimension, source) %>% left_join(he, by=c("individual","element")) %>% mutate(kind="ring"))
add_el(hm %>% filter(size_measure=="periphery_radius_mm", pcount(count)==0) %>%
         transmute(species="Jeholosaurus_shangyuanensis", individual, element, dimension="radius_mm", source=hs, size=val, count=0, efs=2L, note="no LAGs; radius at the periphery"))
add_el(hm %>% filter(size_measure=="element_length_mm") %>%
         transmute(species="Jeholosaurus_shangyuanensis", individual, element=paste0(element, "_length"), dimension="bone_length_mm", source=hs, size=val, note="element length at death"))
add_ind(hm %>% distinct(individual) %>% mutate(species="Jeholosaurus_shangyuanensis"))

## Orodromeus, Dryosaurus (Horner et al. 2009): element length and LAG count.
## Tenontosaurus rows are the same specimens as Werning 2012 and are left out.
ho <- rd("horner2009_hypsilophodontids/individuals.csv") %>% filter(!grepl("Tenontosaurus", species)) %>%
  mutate(size=ifelse(grepl("cm", element_length_as_printed), 10*num(element_length_as_printed), NA))
add_el(ho %>% transmute(species, individual=specimen, element, dimension="bone_length_mm", source="Horner et al. 2009",
        size, count=pcount(n_LAGs_preserved), efs=ifelse(grepl("^TRUE", efs_present), 1L, 0L)))
add_ind(ho %>% distinct(species, individual=specimen))

## Coelophysis (Barta et al. 2022, Supp. Table 1 / Data 1): tibia circumference with the tibial
## growth-mark count, fibular count (no fibula size is tabulated), femur length (measured or
## estimated by the authors; not sectioned). Zone widths (Supp. Table 2) are not used yet.
bs <- "Barta et al. 2022 (Sci Rep) Supp. Table 1"
mi("barta2022_coelophysis", bs, elem=list(tibia_circumference_mm=list(dim="circumference_mm")), keep=function(d) d %>% filter(element != "femur"))
mi("barta2022_coelophysis", bs, elem=list(femur_length_mm=list(dim="bone_length_mm")), keep=function(d) d %>% filter(element == "femur") %>% mutate(count=NA))

## Tyrannosaurus juveniles (Woodward et al. 2020): femur length, CGM counts ("at least").
add_el(tibble(species="Tyrannosaurus_rex", individual=rep(c("BMRP 2002.4.1","BMRP 2006.4.4"), each=2), element=rep(c("femur","tibia"), 2),
              dimension="bone_length_mm", source="Woodward et al. 2020", size=c(688, NA, 774, NA), count=c(13, 10, 15, 15.5),
              count_min=c(TRUE, FALSE, TRUE, FALSE), efs=0L, note=c("length estimated", "count only", NA, "count only; 13 to 18")))
add_ind(tibble(species="Tyrannosaurus_rex", individual=c("BMRP 2002.4.1","BMRP 2006.4.4")))

## Asilisaurus (Griffin & Nesbitt 2016): four sectioned femora, estimated length, no LAGs.
add_el(tibble(species="Asilisaurus_kongwe", individual=paste0("NMT RB", c(211,212,213,226)), element="femur", dimension="bone_length_mm",
              source="Griffin & Nesbitt 2016", size=c(108.7, 113.8, 121.8, 156.6), count=0, efs=0L, note="length estimated by the authors"))
add_ind(tibble(species="Asilisaurus_kongwe", individual=paste0("NMT RB", c(211,212,213,226))))

## Varanopids (Huttenlocker & Shelton 2020): femur length and growth marks from the text.
vs <- "Huttenlocker & Shelton 2020"
add_el(tibble(species=c("Varanops_brevirostris","Mycterosaurus_sp","Mycterosaurinae_indet"),
              individual=c("OMNH 73758","UWBM 89464","Mycterosaurinae_group_mean"), element="femur", dimension="bone_length_mm", source=vs,
              size=c(88, 65, mean(c(49,50,51,57))), count=c(5.5, 1, 3.5), count_min=c(FALSE, TRUE, FALSE), nrep=c(1, 1, 4),
              efs=c(2L, 0L, 0L), note=c("5 to 6 marks; specimen assignment rests on a probable typo (open decision)",
                                         "at least one LAG", "group statement: on average three to four marks, 4 femora")))
add_ind(tibble(species=c("Varanops_brevirostris","Mycterosaurus_sp","Mycterosaurinae_indet"),
               individual=c("OMNH 73758","UWBM 89464","Mycterosaurinae_group_mean")))

## Lystrosaurus (Botha 2020, Supplemental Table S1): 77 elements, 38 animals. Element size =
## midshaft width; body size = % of the species' maximum skull length (Table S1 value).
mi("botha2020_lystrosaurus", "Botha 2020 (PeerJ) Table S1", elem=list(midshaft_width_mm=list(dim="bone_diameter_mm")),
   body=list(pct_max_BSL=list(dim="pct_max_size", el="skull")), keep=function(d) d %>% filter(element != "rib"))

## Therocephalians (Huttenlocker & Botha-Brink 2014 + Supporting Data): counts for 57 elements.
## Element size = midshaft cross-sectional area, entered as the circumference of the circle of
## equal area, sqrt(4 pi CSA), so that it is a linear measure like every other dimension.
## Body size = % of the largest skull of the genus. Moschorhinus rows are left out (held from 2013).
mi("huttenlocker2014_therocephalians", "Huttenlocker & Botha-Brink 2014 (PeerJ) + Supporting Data",
   elem=list(midshaft_csa_mm2=list(dim="circumference_mm", f=function(x) sqrt(4*pi*x))),
   body=list(pct_largest_BSL=list(dim="pct_max_size", el="skull")), keep=function(d) d %>% filter(!grepl("Moschorhinus", species)))

## Gorgonopsians and cynodonts (Botha 2025, Data S2 Table S1): section circumference and growth-mark
## count per element. Left out: specimens not identified to genus ("Gorgonopsia indeterminate",
## "possible new species": 34 elements of unknown and clearly different body sizes cannot share one
## species curve), ribs, and the comparison therocephalians already held from other papers.
mi("botha2025_gorgonopsians", "Botha 2025 (J Anat) Data S2", elem=list("^circumference_mm$"=list(dim="circumference_mm")),
   keep=function(d) d %>% filter(!grepl("indeterminate|new_species|Moschorhinus|scylacosaurid|Lycosuchus", species), element != "rib"))
skip("Botha 2025 Data S2","34 elements of Gorgonopsia indeterminate / possible new species","not identified to genus; different body sizes cannot share one species curve")
skip("Botha 2025 Data S2","zone-width series (214 rows, 55 elements)","radial increments with an unknown inner start; needs a cumulative-radius or increment likelihood (queued)")

## Psittacosaurus lujiatunensis (Zhao et al. 2013, Nat Commun, supplement): limb-bone lengths for 16
## specimens and preserved LAG counts for 42 sections. Zhao 2019's table conflicts with it for
## several specimens and has no sizes; only the 2013 table is used. "IVPP V12617" (Table S1) and
## "IVPP V12716" (Table S2) are kept as printed, i.e. as two animals.
mi("zhao2013_psittacosaurus", "Zhao et al. 2013 (Nat Commun) suppl.", elem=list("_length_mm$"=list(dim="bone_length_mm")))

## Frogs (Laking et al. 2025): SVL and growth-mark count (visual, mean of 4 readers).
fr <- rd("laking2025_oreophrynella_pristimantis/individuals.csv")
add_el(fr %>% transmute(species, individual=paste(genus, specimen), element="femur", dimension="svl_mm", source="Laking et al. 2025",
        size=as.numeric(svl_mm), count=as.numeric(gm_count_visual_mean_of_4_researchers)))
add_ind(fr %>% transmute(species, individual=paste(genus, specimen)))

## Andrias japonicus (Klein et al. 2024): known age and body length. Marks are not used.
an <- rd("klein2024_andrias/individuals.csv") %>% filter(species=="Andrias_japonicus") %>%
  mutate(a=regmatches(age_at_death_as_printed, gregexpr("[0-9]+", age_at_death_as_printed)),
         age=vapply(a, function(v) as.numeric(v[1]) + as.numeric(v[2])/12, 0))
add_el(an %>% transmute(species, individual=specimen, element="body", dimension="total_length_mm", source="Klein et al. 2024",
        size=as.numeric(body_length_mm), t_add=0, resorb=0L, note="known age"))
add_ind(an %>% transmute(species, individual=specimen, age_est=age, age_sd=0.1))

## Lizard skeletochronology: age-class mean SVL. One pseudo-individual per class x sex.
ac <- function(folder, var, src, sp, t_add=0.5, age_sd=0.5){
  d <- rd(folder, "age_class_means.csv") %>% filter(variable==var, !is.na(age), !is.na(mean), as.numeric(n) > 0)
  d <- d %>% mutate(individual=paste0(sp, "_age", age, "_", ifelse(is.na(sex), "both", sex)))
  add_el(d %>% transmute(species=sp, individual, element="body", dimension="svl_mm", source=src, size=as.numeric(mean),
                         nrep=as.numeric(n), t_add=t_add, resorb=0L, note="age-class mean SVL"))
  add_ind(d %>% transmute(species=sp, individual, age_est=as.numeric(age), age_sd=age_sd)) }
ac("gharzi2013_ophisops_elegans", "SVL_mean", "Gharzi & Yari 2013", "Ophisops_elegans")
ac("ma2022_phrynocephalus_vlangalii", "SVL", "Ma et al. 2022", "Phrynocephalus_vlangalii")
## Phymaturus (Cabezas-Cartes et al. 2018): group means; age = the authors' estimated age
## (LAGs + back-calculated resorbed marks), taken as the age at death.
ph <- rd("cabezas2018_phymaturus_aguanegra_zapalensis/age_class_means.csv") %>% mutate(mean=as.numeric(mean)) %>%
  mutate(v=ifelse(variable=="SVL", "svl", ifelse(grepl("^estimated_age", variable), "age", "lags"))) %>%
  select(species, group, sex, n, v, mean) %>% tidyr::pivot_wider(names_from=v, values_from=mean) %>%
  mutate(species=gsub(" ", "_", species), individual=paste(species, group, sex, sep="_"))
add_el(ph %>% transmute(species, individual, element="body", dimension="svl_mm", source="Cabezas-Cartes et al. 2018", size=svl,
                        nrep=as.numeric(n), t_add=0, resorb=0L, note="group mean SVL and mean estimated age"))
add_ind(ph %>% transmute(species, individual, age_est=age, age_sd=1))

## Caiman latirostris (Mascarenhas-Junior et al. 2021): 10 farm-reared females of known age (whole
## years, 1-6). Humerus length at a known age. Growth-mark counts are NOT used: they fall short of
## the age (R7, marks are not annual); held for the marks-per-year calibration.
cm <- rd("mascarenhasjunior2021_caiman_latirostris/model_inputs.csv") %>% filter(size_measure=="humerus_length_mm")
add_el(cm %>% transmute(species="Caiman_latirostris", individual=paste("Clat", individual), element="humerus", dimension="bone_length_mm",
        source="Mascarenhas-Junior et al. 2021", size=as.numeric(size_value), t_add=0, resorb=0L, efs=0L, note="known age (year class)"))
add_ind(cm %>% transmute(species="Caiman_latirostris", individual=paste("Clat", individual), age_est=as.numeric(age_known_yr), age_sd=0.3))

## Sauromalus ater (Petermann & Gauthier 2020, Table 2): SVL and the authors' age (LAG count; 13 of
## 31 are starred "first LAG lost"). Entered as source-estimated ages: sd 0.5 yr, 1 yr if starred.
sa <- rd("petermann2020_sauromalus/model_inputs.csv") %>% mutate(age=as.numeric(sub(".*table2_age_years=([0-9.]+).*", "\\1", note)), star=is.na(count))
add_el(sa %>% transmute(species="Sauromalus_ater", individual, element="body", dimension="svl_mm", source="Petermann & Gauthier 2020",
        size=as.numeric(size_value), resorb=0L, note=ifelse(star, "age from Table 2; first LAG lost", "age from Table 2")))
add_ind(sa %>% transmute(species="Sauromalus_ater", individual, age_est=age, age_sd=ifelse(star, 1, 0.5)))
skip("Schlief et al. 2025 (Eublepharis)","13 captive geckos, 5 of known age, no LAGs in any element","no body or bone size per individual in the tables; relevant to R7 (no marks at up to 3.9 yr)")
skip("Barta 2022 Supp. Tables 2-4","zone widths: 21 Coelophysis tibiae, 10 Alligator femora, 26 Maiasaura tibiae","radial increments from an unspecified origin; needs an increment likelihood (queued)")

## Left out, with reasons
skip("Zhao et al. 2019 (Psittacosaurus)","43 element LAG counts","superseded by Zhao 2013 supplement (sizes + counts); the two tables conflict for several specimens")
skip("Bybee 2006 individuals","bone length/circumference at death of the 20 sectioned bones","already ring series; unsectioned bones have no count")
skip("Zug & Rand 1987 (Iguana)","SVL by age of first-time nesting females","conditioned on first nesting: not a growth cross-section")
skip("Roitberg & Smirina 2006; El Mouden 1999","LAG1 diameters; mark-recapture growth rates","no size-at-age table (figures need digitizing)")
skip("Ong 2022 (ratites); Buso 2026 (Nile crocodile adult)","known-age animals","no size measure extracted")
skip("Goldsmith 2025; Guo 2026; Horner 2009 Camptosaurus","single sized specimens","fewer than MIN_OBS sized observations for the taxon")
skip("Navarro 2025; Pereyra 2025; Maillo 2025; Toefy 2025; Weiss 2022; Skutschas 2015; Whitney 2021, 2022; Prondvai 2018","counts per element","no size comparable across specimens in the tables")

# ============ assemble ================================================================
rings <- bind_rows(rings) %>% mutate(el_id=paste(individual, element, dimension, sep=" | "))
els <- bind_rows(els) %>% mutate(el_id=paste(individual, element, dimension, sep=" | "), eltype=eltype_of(element))
inds <- bind_rows(inds) %>% distinct(species, individual, .keep_all=TRUE)
stopifnot(!anyDuplicated(els$el_id), !anyDuplicated(inds$individual), all(els$individual %in% inds$individual),
          all(rings$el_id %in% els$el_id[els$kind=="ring"]), all(els$el_id[els$kind=="ring"] %in% rings$el_id))
# count intervals from the printed text: "5 to 6", "2 or 3", "13 to 18" -> [lo, hi]; "at least" -> hi = Inf
cint <- function(note, count, cmin){
  lo <- count; hi <- count
  for (r in seq_along(note)) { if (is.na(count[r]) || is.na(note[r])) next
    txt <- if (grepl("as printed: ", note[r])) sub(".*as printed: ", "", note[r]) else note[r]
    mt <- regmatches(txt, regexpr("[0-9]+(\\.[0-9]+)? *(to|or) *[0-9]+(\\.[0-9]+)?", txt))
    if (length(mt)) { v <- as.numeric(regmatches(mt, gregexpr("[0-9]+(\\.[0-9]+)?", mt))[[1]]); lo[r] <- min(v); hi[r] <- max(v) }
    if (grepl("at least", txt)) hi[r] <- Inf }
  hi[!is.na(count) & cmin] <- Inf
  list(lo=lo, hi=hi) }
ci <- cint(els$note, els$count, els$count_min); els$count_lo <- ci$lo; els$count_hi <- ci$hi
# sex, where the source gives it
sx <- bind_rows(
  tibble(individual=inds$individual, sex=sub(".*_(male|female)$", "\\1", inds$individual)) %>% filter(grepl("_(male|female)$", individual)),
  rd("laking2025_oreophrynella_pristimantis/individuals.csv") %>% transmute(individual=paste(genus, specimen), sex=tolower(sex)),
  rd("klein2024_andrias/individuals.csv") %>% transmute(individual=specimen, sex=ifelse(sex %in% c("female","male"), sex, NA)),
  rd("chinsamy2026_crocodylus_niloticus/individuals.csv") %>% transmute(individual=paste("Cnil", specimen), sex=tolower(sex)),
  rd("petermann2020_sauromalus/model_inputs.csv") %>% transmute(individual, sex=c(F="female", M="male")[sub("^sex ([FM]);.*", "\\1", note)]),
  tibble(individual=inds$individual[grepl("^Clat ", inds$individual)], sex="female")) %>%
  filter(sex %in% c("female","male")) %>% distinct(individual, .keep_all=TRUE)
inds$sex <- sx$sex[match(inds$individual, sx$individual)]
# drop cs elements carrying neither a size nor a count, then individuals with no sized observation
els <- els %>% filter(kind=="ring" | !is.na(size) | !is.na(count))
sized <- union(rings$individual, els$individual[!is.na(els$size)])
skip("(all)", paste(setdiff(inds$individual, sized), collapse="; "), "individual has no sized observation")
els <- els %>% filter(individual %in% sized); inds <- inds %>% filter(individual %in% sized)
# individuals with a count-based age need at least one count or ring series
okage <- inds$individual[!is.na(inds$age_est)]
okcnt <- union(rings$individual, els$individual[!is.na(els$count)])
bad <- setdiff(inds$individual, union(okage, okcnt))
if (length(bad)) { skip("(all)", paste(bad, collapse="; "), "individual has sizes but no growth-mark count, ring series or age")
  els <- els %>% filter(!individual %in% bad); inds <- inds %>% filter(!individual %in% bad); rings <- rings %>% filter(!individual %in% bad) }
# species with too few sized observations
nobs <- bind_rows(rings %>% count(species), els %>% filter(!is.na(size)) %>% count(species)) %>% group_by(species) %>% summarise(n=sum(n))
drop_sp <- nobs$species[nobs$n < MIN_OBS]
if (length(drop_sp)) skip("(all)", paste(drop_sp, collapse="; "), sprintf("fewer than %d sized observations", MIN_OBS))
rings <- rings %>% filter(!species %in% drop_sp); els <- els %>% filter(!species %in% drop_sp); inds <- inds %>% filter(!species %in% drop_sp)

# ============ tree (genus-level topology, hand-assembled; no branch lengths) ===========
# Dinosaur and crocodylian part = data_clean/dino_tree_next.nwk (src/28_build_trees.R) plus:
#  Coelophysis sister to Syntarsus (Coelophysoidea); Ceratosaurus sister to Tetanurae.
#  Ornithischia: Jeholosaurus and Orodromeus as successive outgroups to Cerapoda (after
#   Boyd 2015; CONTESTED: some analyses place both in Ornithopoda); Tenontosaurus and
#   Dryosaurus as successive outgroups to Hadrosauridae.
#  Asilisaurus (Silesauridae) sister to Dinosauria (Nesbitt et al. 2010; CONTESTED: silesaurids
#   as early ornithischians, Mueller & Garcia 2020).
#  Orthosuchus (Protosuchia) sister to Crocodylia; Crocodylus sister to Alligatoroidea.
#  Lepidosauria: Ophisops (Lacertidae) sister to Iguania (Phymaturus, Phrynocephalus).
#  Synapsida: varanopids (Mycterosaurus, Varanops) as basal synapsids (CONTESTED: Ford &
#   Benson 2020 recover Varanopidae among sauropsids); Lystrosaurus (Anomodontia) sister to
#   Moschorhinus (Therocephalia) within Therapsida.
#  Lissamphibia: Andrias sister to Anura (Oreophrynella, Pristimantis).
orn   <- "(Jeholosaurus,(Orodromeus,(Psittacosaurus,(Tenontosaurus,(Dryosaurus,(Maiasaura,Hypacrosaurus))))))"
sauro <- "((Plateosaurus,Massospondylus),(Mamenchisaurus,(Northampton,(Janenschia,((Apatosaurus,Diplodocus),(Camarasaurus,Alamosaurus))))))"
tet   <- "((Allosaurus,(Acrocanthosaurus,Meraxes)),(((Albertosaurus,Gorgosaurus),Tyrannosaurus),(Rativates,Saurornitholestes)))"
thero <- sprintf("((Syntarsus,Coelophysis),(Ceratosaurus,%s))", tet)
dino  <- sprintf("(%s,(%s,%s))", orn, sauro, thero)
#  Caiman sister to Alligator (Alligatoridae), Diplocynodon outside them (basal alligatoroid). Sauromalus (Iguanidae) with
#   Phymaturus (Liolaemidae) in Pleurodonta.
arch  <- sprintf("((Orthosuchus,(Crocodylus,(Diplocynodon,(Alligator,Caiman)))),(Asilisaurus,%s))", dino)
saurop<- sprintf("((Ophisops,((Phymaturus,Sauromalus),Phrynocephalus)),%s)", arch)
#  Therapsida (added 2026-11-14, provisional): Anomodontia (Lystrosaurus) sister to Theriodontia;
#   Gorgonopsia sister to Eutheriodontia (Therocephalia + Cynodontia). Therocephalia after
#   Huttenlocker & Botha-Brink 2014's framework: Lycosuchus, then Scylacosauridae (+ Glanosuchus),
#   then Eutherocephalia = (Akidnognathidae (Moschorhinus, Olivierosuchus) + Hofmeyriidae
#   (Hofmeyria, Mirotenthes)) and (Theriognathus, (Ictidosuchoides, Baurioidea (Tetracynodon,
#   (Scaloposaurus, Microgomphodon)))). Gorgonopsian interrelationships are not resolved here.
thero_s <- "(Lycosuchus,((Scylacosauridae,Glanosuchus),(((Moschorhinus,Olivierosuchus),(Hofmeyria,Mirotenthes)),(Theriognathus,(Ictidosuchoides,(Tetracynodon,(Scaloposaurus,Microgomphodon)))))))"
gorg  <- "(Cyonosaurus,Gorgonops,Aelurognathus,Tigricephalus,Arctops,Inostrancevia)"
synap <- sprintf("((Mycterosaurus,Varanops),(Lystrosaurus,(%s,(%s,(Cynognathus,Diademodon)))))", gorg, thero_s)
tree  <- read.tree(text=sprintf("((Andrias,(Oreophrynella,Pristimantis)),(%s,%s));", synap, saurop))
tm <- read_csv(file.path(clean,"species_to_tip_next.csv"), show_col_types=FALSE) %>% bind_rows(tribble(~species, ~phylo_tip,
  "Ceratosaurus_nasicornis","Ceratosaurus", "Ceratosaurus_magnicornis","Ceratosaurus", "Ceratosaurus_sp","Ceratosaurus",
  "Coelophysis_bauri","Coelophysis", "Asilisaurus_kongwe","Asilisaurus", "Orthosuchus_stormbergi","Orthosuchus",
  "Crocodylus_niloticus","Crocodylus", "Tenontosaurus_tilletti","Tenontosaurus", "Jeholosaurus_shangyuanensis","Jeholosaurus",
  "Orodromeus_makelai","Orodromeus", "Dryosaurus_altus","Dryosaurus", "Moschorhinus_kitchingi","Moschorhinus",
  "Lystrosaurus_maccaigi","Lystrosaurus", "Lystrosaurus_curvatus","Lystrosaurus", "Lystrosaurus_murrayi","Lystrosaurus",
  "Lystrosaurus_declivis","Lystrosaurus", "Varanops_brevirostris","Varanops", "Mycterosaurus_sp","Mycterosaurus",
  "Mycterosaurinae_indet","Mycterosaurus", "Ophisops_elegans","Ophisops", "Phrynocephalus_vlangalii","Phrynocephalus",
  "Phymaturus_aguanegra","Phymaturus", "Phymaturus_zapalensis","Phymaturus", "Andrias_japonicus","Andrias",
  "Oreophrynella_quelchii","Oreophrynella", "Pristimantis_aureoventris","Pristimantis",
  "Lycosuchus_vanderrieti","Lycosuchus", "Scylacosauridae_indet","Scylacosauridae", "Glanosuchus_macrops","Glanosuchus",
  "Olivierosuchus_parringtoni","Olivierosuchus", "Hofmeyria_atavus","Hofmeyria", "Mirotenthes_digitipes","Mirotenthes",
  "Theriognathus_microps","Theriognathus", "Ictidosuchoides_longiceps","Ictidosuchoides", "Tetracynodon_darti","Tetracynodon",
  "Scaloposaurus_constrictus","Scaloposaurus", "Microgomphodon_oligocynus","Microgomphodon",
  "Cyonosaurus","Cyonosaurus", "Gorgonops_torvus","Gorgonops", "Aelurognathus_tigriceps","Aelurognathus",
  "Tigricephalus_kingwilli","Tigricephalus", "Arctops_willistoni","Arctops", "Inostrancevia_africana","Inostrancevia",
  "Cynognathus_crateronotus","Cynognathus", "Diademodon_tetragonus","Diademodon",
  "Caiman_latirostris","Caiman", "Sauromalus_ater","Sauromalus")) %>% distinct(species, .keep_all=TRUE)
spp <- sort(unique(inds$species))
if (!all(spp %in% tm$species)) stop("species without a tip: ", paste(setdiff(spp, tm$species), collapse=", "))
tm <- tm %>% filter(species %in% spp)
tree <- keep.tip(tree, intersect(tree$tip.label, tm$phylo_tip))
stopifnot(all(tm$phylo_tip %in% tree$tip.label), all(tree$tip.label %in% tm$phylo_tip))

write_csv(rings, file.path(clean,paste0(PRE,"_rings.csv"))); write_csv(els, file.path(clean,paste0(PRE,"_elements.csv")))
write_csv(inds, file.path(clean,paste0(PRE,"_individuals.csv"))); write_csv(bind_rows(excl), file.path(clean,paste0(PRE,"_excluded.csv")))
write_csv(myh_class, file.path(clean,paste0(PRE,"_myhrvold_classes.csv")))
write.tree(tree, file.path(clean,paste0(PRE,"_tree.nwk"))); write_csv(tm, file.path(clean,paste0(PRE,"_species_to_tip.csv")))

cat(sprintf("species %d | tips %d | individuals %d (%d with an age estimate) | elements %d (%d ring, %d count-and-size) | sized marks %d\n",
            length(spp), Ntip(tree), nrow(inds), sum(!is.na(inds$age_est)), nrow(els), sum(els$kind=="ring"), sum(els$kind=="cs"), nrow(rings)))
ne <- els %>% count(individual); cat(sprintf("individuals with more than one element: %d\n", sum(ne$n > 1)))
cat("EFS states (elements):"); print(table(factor(els$efs, 0:2, c("absent","present","unknown")), els$kind))
options(width=200)
print(as.data.frame(els %>% group_by(species) %>% summarise(inds=n_distinct(individual), ring_el=sum(kind=="ring"), cs_sized=sum(kind=="cs" & !is.na(size)),
      cs_count_only=sum(kind=="cs" & is.na(size)), dims=paste(unique(dimension), collapse="|")) %>%
      left_join(rings %>% count(species, name="marks"), by="species")))
