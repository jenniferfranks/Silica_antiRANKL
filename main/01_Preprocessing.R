# -------------
# Preprocessing Silica_antiRANKL Data
#
# - Jennifer Franks
# - April 15, 2026
#
# - Read in Nextseq data, basic QC filtering
#
# --------------

library(dplyr)
library(monocle3)
library(Matrix)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(viridis)
library(cowplot)
library(ggplotify)
library(msigdbr)


main_dir <- "/data/hps/assoc/private/franks_lab/user/jfra11/projects/Silica/antiRANKL/"

setwd(main_dir)

cds <- readRDS("data/Jennifer_Franks/Project_159/Franks_RNA3-116-a_nextseq_data/monocle3_cds/all_Mouse_cds.RDS")
samplemap <- read.csv("data/Jennifer_Franks/Project_159/Franks_RNA3-116-a_nextseq_data/Franks_SampleIDMap.csv")

# format BAT sample IDs consistently
samplemap$BAT.Sample.ID_chr <- sprintf("%.4f", as.numeric(samplemap$BAT.Sample.ID))

# named lookup vector
converter <- setNames(
  as.character(samplemap$Sample.ID),
  samplemap$BAT.Sample.ID_chr
)

# map cds sample values back to Sample.ID
pData(cds)$BAT.Sample.ID <- sprintf("%.4f", as.numeric(as.character(pData(cds)$sample)))
pData(cds)$Sample.ID <- converter[pData(cds)$BAT.Sample.ID]

# initialize new columns
pData(cds)$mouse_id  <- NA_character_
pData(cds)$timepoint <- NA_character_
pData(cds)$treatment <- NA_character_
pData(cds)$exposure  <- NA_character_
pData(cds)$condition <- NA_character_

# identify controls
control_idx <- grepl("control$", pData(cds)$Sample.ID)

# controls: WT, no silica, no treatment
pData(cds)$mouse_id[control_idx]  <- sub("_control$", "", pData(cds)$Sample.ID[control_idx])
pData(cds)$timepoint[control_idx] <- "control"
pData(cds)$treatment[control_idx] <- "none"
pData(cds)$exposure[control_idx]  <- "WT_no_silica"
pData(cds)$condition[control_idx] <- "control"

# non-controls
noncontrol_idx <- !control_idx & !is.na(pData(cds)$Sample.ID)

tmp <- strcapture(
  "^(\\d+)_(d\\d+)_(IgG|antiRANKL)$",
  pData(cds)$Sample.ID[noncontrol_idx],
  proto = list(
    mouse_id = character(),
    timepoint = character(),
    treatment = character()
  )
)

pData(cds)$mouse_id[noncontrol_idx]  <- tmp$mouse_id
pData(cds)$timepoint[noncontrol_idx] <- tmp$timepoint
pData(cds)$treatment[noncontrol_idx] <- tmp$treatment
pData(cds)$exposure[noncontrol_idx]  <- "silica"
pData(cds)$condition[noncontrol_idx] <- paste(tmp$timepoint, tmp$treatment, sep = "_")


# ------------- EVERYTHING ---------------- 

# QUALITY CONTROL #
# parameters ---
PCs_to_use <- 100 
#UMImax <- min(quantile(pData(cds)$n.umi, 0.99), quantile(pData(cds)$n.umi, 0.75)+3*iqr(pData(cds)$n.umi))
UMImin <- 500 
mitoMax <- 10. 
cluster.k <- 50 # full dataset
scrublet.score.max <- 0.2 
set.seed(8248959)
# ---

# Doublet filtering
cds <- cds[,which(pData(cds)$scrublet_score < scrublet.score.max)]

# UMI filtering 
cds <- cds[,which(pData(cds)$n.umi > UMImin)]
cds <- cds[,which(pData(cds)$n.umi < UMImax)]

# Mitochondrial reads filtering 
cds <- cds[,pData(cds)$perc_mitochondrial_umis < mitoMax]

# --------------------

cds <- cds %>% 
  estimate_size_factors() %>%
  preprocess_cds(num_dim=PCs_to_use) %>% 
  detect_genes() %>%
  reduce_dimension() %>%
  cluster_cells(k=cluster.k, random_seed = 3752)

pData(cds)$cluster <- clusters(cds)

plot_cells(cds, color_cells_by = "cluster", label_cell_groups = TRUE, group_label_size=10) + 
  theme_void()

plot_cells(cds, color_cells_by = "sample", label_cell_groups = FALSE) + 
  theme_void()

plot_cells(cds, color_cells_by = "exposure", label_cell_groups = FALSE) + 
  theme_void()

plot_cells(cds, color_cells_by = "treatment", label_cell_groups = FALSE) + 
  theme_void()

plot_cells(cds, color_cells_by = "timepoint", label_cell_groups = FALSE) + 
  theme_void()

plot_cells(cds, color_cells_by = "scrublet_score", label_cell_groups = FALSE) + 
  theme_void()

table(pData(cds)$treatment, pData(cds)$timepoint)

plot_cells(cds, 
  genes=c("Sftpc", "Siglecf", "Pecam1", "Ptprc", "Cd68", "Col1a1"), 
  label_cell_groups = FALSE, scale_to_range = FALSE) 


cds2 <- cds[,which(pData(cds)$cluster %in% c(3,9))]
cds2 <- cds2 %>% 
  preprocess_cds(num_dim=PCs_to_use) %>% 
  reduce_dimension() %>%
  cluster_cells(k=10, random_seed = 3752)
pData(cds2)$myeloidcluster <- clusters(cds2)

plot_cells(cds2, color_cells_by = "myeloidcluster", label_cell_groups = FALSE, cell_size=1.5) + 
  theme_void() 
plot_cells(cds2, color_cells_by = "timepoint", label_cell_groups = FALSE, cell_size=1.5) + 
  theme_void() 
plot_cells(cds2, color_cells_by = "treatment", label_cell_groups = FALSE, cell_size=1.5) + 
  theme_void() 
plot_cells(cds2, color_cells_by = "exposure", label_cell_groups = FALSE, cell_size=1.5) + 
  theme_void() 


plot_cells(cds2, 
  genes=c("Cd68", "Siglecf", "Marco","Adgre4", "Itgam", "Ccr2", "Cd74"), 
  label_cell_groups = FALSE, scale_to_range = FALSE, cell_size=1.5)

plot_genes_by_group(cds2, markers=c("Cd68", "Siglecf","Marco", "Adgre4", "Itgam", "Ccr2", "Cd74"), 
  group_cells_by="myeloidcluster", ordering_type="maximal_on_diag")


table(pData(cds2)$myeloidcluster, pData(cds2)$timepoint)
table(pData(cds2)$myeloidcluster, pData(cds2)$condition)
