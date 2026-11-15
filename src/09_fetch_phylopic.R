#!/usr/bin/env Rscript
# Fetch + cache PhyloPic silhouettes for every tip in data_clean/dino_tree.nwk.
# Saves each as a raster RDS in data_clean/phylopic/ + a uuid table.
suppressMessages({ library(rphylopic); library(ape) })
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
dir  <- file.path(root, "data_clean/phylopic"); dir.create(dir, showWarnings = FALSE, recursive = TRUE)
tips <- read.tree(file.path(root, "data_clean/dino_tree.nwk"))$tip.label

# fallbacks for tips PhyloPic lacks at genus level
fb <- list(Meraxes=c("Giganotosaurus","Carcharodontosauridae"),
           Rativates=c("Ornithomimus","Ornithomimidae"),
           Northampton=c("Cetiosaurus","Sauropoda"), Janenschia=c("Turiasaurus","Sauropoda"),
           Alamosaurus=c("Saltasaurus","Titanosauria","Sauropoda"),
           Syntarsus=c("Coelophysis","Coelophysoidea"),
           Saurornitholestes=c("Velociraptor","Dromaeosauridae"),
           Mamenchisaurus=c("Mamenchisaurus","Mamenchisauridae"))

rows <- list()
for (g in tips) {
  cands <- c(g, fb[[g]]); uuid <- NA; used <- NA
  for (nm in cands) { uuid <- tryCatch(get_uuid(name=nm, n=1), error=function(e) NA)
    if (!is.na(uuid)) { used <- nm; break } }
  ok <- FALSE
  if (!is.na(uuid)) { img <- tryCatch(get_phylopic(uuid=uuid, format="raster"), error=function(e) NULL)
    if (!is.null(img)) { saveRDS(img, file.path(dir, paste0(g,".rds"))); ok <- TRUE } }
  rows[[g]] <- data.frame(genus=g, search=ifelse(is.na(used),NA,used), uuid=ifelse(is.na(uuid),NA,uuid), cached=ok)
  cat(sprintf("%-18s -> %-22s %s\n", g, ifelse(is.na(used),"(none)",used), ifelse(ok,"cached","FAILED")))
}
ut <- do.call(rbind, rows); write.csv(ut, file.path(dir,"phylopic_uuids.csv"), row.names=FALSE)
cat(sprintf("\n%d/%d silhouettes cached\n", sum(ut$cached), nrow(ut)))
