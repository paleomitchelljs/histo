#!/usr/bin/env Rscript
# Fetch PhyloPic silhouettes for tips of data_clean/broad_tree.nwk that are not cached yet.
# Appends to data_clean/phylopic/phylopic_uuids.csv; existing files are left alone.
suppressMessages({ library(rphylopic); library(ape) })
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) }); dir <- file.path(root, "data_clean/phylopic")
tips <- read.tree(file.path(root, "data_clean/broad_tree.nwk"))$tip.label
tips <- tips[!file.exists(file.path(dir, paste0(tips, ".rds")))]
fb <- list(Diplocynodon=c("Alligator"), Orthosuchus=c("Protosuchus","Crocodyliformes"), Moschorhinus=c("Therocephalia","Theriodontia"),
           Asilisaurus=c("Silesaurus","Silesauridae"), Orodromeus=c("Thescelosaurus","Hypsilophodon"), Jeholosaurus=c("Hypsilophodon","Thescelosaurus"),
           Oreophrynella=c("Bufonidae","Bufo"), Pristimantis=c("Eleutherodactylus","Anura"), Ophisops=c("Lacertidae","Lacerta"),
           Phrynocephalus=c("Agamidae","Agama"), Phymaturus=c("Liolaemus","Liolaemidae","Iguania"), Andrias=c("Cryptobranchidae"),
           Dryosaurus=c("Dryosauridae","Iguanodontia"), Tenontosaurus=c("Iguanodontia"), Lystrosaurus=c("Dicynodontia"),
           # 2026-11-14 v3 additions: genera with no PhyloPic image borrow a clade image
           Cyonosaurus=c("Gorgonopsia","Gorgonops","Inostrancevia"), Tigricephalus=c("Therocephalia","Theriodontia"),
           Lycosuchus=c("Therocephalia","Theriodontia"), Glanosuchus=c("Therocephalia","Theriodontia"),
           Olivierosuchus=c("Therocephalia","Theriodontia"), Hofmeyria=c("Therocephalia","Theriodontia"),
           Mirotenthes=c("Therocephalia","Theriodontia"), Theriognathus=c("Therocephalia","Theriodontia"),
           Ictidosuchoides=c("Therocephalia","Theriodontia"), Tetracynodon=c("Therocephalia","Theriodontia"),
           Diademodon=c("Cynognathus","Cynodontia"), Sauromalus=c("Iguanidae","Ctenosaura","Iguana","Iguania"))
rows <- list()
for (g in tips) { uuid <- NA; used <- NA
  for (nm in c(g, fb[[g]])) { uuid <- tryCatch(get_uuid(name=nm, n=1), error=function(e) NA); if (!is.na(uuid)) { used <- nm; break } }
  ok <- FALSE
  if (!is.na(uuid)) { img <- tryCatch(get_phylopic(uuid=uuid, format="raster"), error=function(e) NULL)
    if (!is.null(img)) { saveRDS(img, file.path(dir, paste0(g, ".rds"))); ok <- TRUE } }
  rows[[g]] <- data.frame(genus=g, search=used, uuid=uuid, cached=ok)
  cat(sprintf("%-16s -> %-20s %s\n", g, ifelse(is.na(used), "(none)", used), ifelse(ok, "cached", "FAILED"))) }
if (length(rows)) { f <- file.path(dir, "phylopic_uuids.csv"); old <- read.csv(f)
  write.csv(rbind(old[!old$genus %in% names(rows),], do.call(rbind, rows)), f, row.names=FALSE) }
