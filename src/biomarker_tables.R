library(tidyverse)
library(readr)
library(magrittr)

expression_dir <- # INSERT path to folder with general expression data
response_dir <- # INSERT path to folder with project specific results

# load the data (public versions)
GE <- data.table::fread(paste0(expression_dir, "/ge.csv"))
CNA <- data.table::fread(paste0(expression_dir, "/cna.csv"))
MET <- data.table::fread(paste0(expression_dir, "/met.csv"))
miRNA <- data.table::fread(paste0(expression_dir, "/mirna.csv"))
MUT <- data.table::fread(paste0(expression_dir, "/mut.csv"))
PROT <- data.table::fread(paste0(expression_dir, "/prot.csv"))
XPR <- data.table::fread(paste0(expression_dir, "/xpr.csv"))
LIN <- data.table::fread(paste0(expression_dir, "/lin.csv"))
shRNA <- data.table::fread(paste0(expression_dir, "/shrna.csv"))
repurposing <- data.table::fread(paste0(expression_dir, "/rep.csv"))


# data table mapping depmap to arxspan id
name_map <- LIN %>%
  dplyr::select(DepMap_ID, CCLE_Name)
  dplyr::rename(ccle_name = CCLE_Name, broad_id = DepMap_ID) %>%
  dplyr::distinct(ccle_name, broad_id)
cell_lines <- data.table::fread(paste0(response_dir, "/logMFI.csv")) %>%
  dplyr::distinct(ccle_name) %>%
  dplyr::left_join(name_map)


# for each table:
# 1) rename columns (add feature abbreviation to front)
# 2) rename rows by joining with name tables

# shRNA
shRNA <- t(shRNA)
colnames(shRNA) <- paste("shRNA", word(colnames(shRNA), 1, sep = " "), sep = "_")
shRNA <- as_tibble(shRNA, rownames = "broad_id", .name_repair = make.names) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(cell_lines) %>%
  dplyr::select(-broad_id) %>%
  dplyr::distinct() %>%
  column_to_rownames("ccle_name")

# miRNA
colnames(miRNA) <- paste("miRNA", colnames(miRNA), sep = "_")
miRNA <- as_tibble(miRNA, rownames = "broad_id") %>%
  dplyr::distinct() %>%
  dplyr::inner_join(cell_lines) %>%
  dplyr::select(-broad_id) %>%
  dplyr::distinct() %>%
  column_to_rownames("ccle_name")

# metabolomics
colnames(MET) <- paste("MET", colnames(MET), sep = "_")
MET <- as_tibble(MET, rownames = "broad_id") %>%
  dplyr::distinct() %>%
  dplyr::inner_join(cell_lines) %>%
  dplyr::select(-broad_id) %>%
  dplyr::distinct() %>%
  column_to_rownames("ccle_name")

# protein expression
colnames(RPPA) <- paste("RPPA", colnames(RPPA), sep = "_")
RPPA <- as_tibble(RPPA, rownames = "broad_id") %>%
  dplyr::distinct() %>%
  dplyr::inner_join(cell_lines) %>%
  dplyr::select(-broad_id) %>%
  dplyr::distinct() %>%
  column_to_rownames("ccle_name")

# lineage table (convert to matrix from long)
LIN %<>% dplyr::rename(ccle_name = `CCLE Name`) %>%
  dplyr::filter(!grepl("MERGED", ccle_name))
lin <- LIN %>% dplyr::distinct(ccle_name, lineage) %>%
  dplyr::filter(lineage != "")
sub <- LIN %>% dplyr::distinct(ccle_name, lineage_subtype) %>%
  dplyr::rename(lineage = lineage_subtype) %>%
  dplyr::filter(lineage != "")
sub2 <- LIN %>% dplyr::distinct(ccle_name, lineage_sub_subtype) %>%
  dplyr::rename(lineage = lineage_sub_subtype) %>%
  dplyr::filter(lineage != "")
LIN_table <- dplyr::bind_rows(lin, sub, sub2) %>%
  dplyr::distinct() %>%
  dplyr::mutate(from = 1)
LIN <- reshape2::acast(LIN_table, ccle_name ~ lineage,
                       value.var = "from", fill = 0)
colnames(LIN) <- paste("LIN", colnames(LIN), sep = "_")


colnames(GE) <- paste("GE", word(colnames(GE), 1, sep = " "), sep = "_")
GE <- as_tibble(GE, rownames = "broad_id", .name_repair = make.names) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(cell_lines) %>%
  dplyr::select(-broad_id) %>%
  dplyr::distinct() %>%
  column_to_rownames("ccle_name")

# gene effect (CRISPR)
colnames(XPR) <- paste("XPR", word(colnames(XPR), 1, sep = " "), sep = "_")
XPR <- as_tibble(XPR, rownames = "broad_id", .name_repair = make.names) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(cell_lines) %>%
  dplyr::select(-broad_id) %>%
  dplyr::distinct() %>%
  column_to_rownames("ccle_name")

# copy number
colnames(CNA) <- paste("CNA", word(colnames(CNA), 1, sep = " "), sep = "_")
CNA <- as_tibble(CNA, rownames = "broad_id") %>%
  dplyr::distinct() %>%
  dplyr::inner_join(cell_lines) %>%
  dplyr::select(-broad_id) %>%
  dplyr::distinct() %>%
  column_to_rownames("ccle_name")

# mutations (damaging, hotspot, other)
colnames(damMUT) <- paste("MUT_dam", word(colnames(damMUT), 1, sep = " "), sep = "_")
damMUT <- as_tibble(damMUT, rownames = "broad_id", .name_repair = make.names) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(cell_lines) %>%
  dplyr::select(-broad_id) %>%
  dplyr::distinct()
colnames(hsMUT) <- paste("MUT_hs", word(colnames(hsMUT), 1, sep = " "), sep = "_")
hsMUT <- as_tibble(hsMUT, rownames = "broad_id", .name_repair = make.names) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(cell_lines) %>%
  dplyr::select(-broad_id) %>%
  dplyr::distinct()
colnames(otherMUT) <- paste("MUT_other", word(colnames(otherMUT), 1, sep = " "), sep = "_")
otherMUT <- as_tibble(otherMUT, rownames = "broad_id", .name_repair = make.names) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(cell_lines) %>%
  dplyr::select(-broad_id) %>%
  dplyr::distinct()
MUT <- dplyr::inner_join(damMUT, hsMUT) %>%
  dplyr::inner_join(otherMUT) %>%
  column_to_rownames("ccle_name")

# repurposing
rep_info %<>% column_to_rownames(var = "column_name")
colnames(REP) <- paste("REP", rep_info[colnames(REP), ]$name, sep = "_")
REP <- as_tibble(REP, rownames = "broad_id", .name_repair = make.names) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(cell_lines) %>%
  dplyr::select(-broad_id) %>%
  dplyr::distinct() %>%
  column_to_rownames("ccle_name")

# total proteome
prots <- PROT$Uniprot
PROT %<>% dplyr::select(-c(1:4), -c(6:48))
colnames(PROT) <- word(colnames(PROT), start = 1, end = -2, sep = "_")
PROT <- t(PROT[,-1])
colnames(PROT) <- paste("PROT", prots, sep = "_")
PROT <- as_tibble(PROT, rownames = "ccle_name", .name_repair = make.names) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(cell_lines) %>%
  dplyr::select(-broad_id) %>%
  dplyr::distinct() %>%
  column_to_rownames("ccle_name")


# filter out variance
MUT <- MUT[, apply(MUT, 2, function(x) var(x, na.rm = TRUE)) > 0]
GE <- GE[, apply(GE, 2, function(x) var(x, na.rm = TRUE)) > 0]
shRNA <- shRNA[, apply(shRNA, 2, function(x) var(x, na.rm = TRUE)) > 0]
CNA <- CNA[, apply(CNA, 2, function(x) var(x, na.rm = TRUE)) > 0]
XPR <- XPR[, apply(XPR, 2, function(x) var(x, na.rm = TRUE)) > 0]
LIN <- LIN[, apply(LIN, 2, function(x) var(x, na.rm = TRUE)) > 0]
miRNA <- miRNA[, apply(miRNA, 2, function(x) var(x, na.rm = TRUE)) > 0]
RPPA <- RPPA[, apply(RPPA, 2, function(x) var(x, na.rm = TRUE)) > 0]
MET <- MET[, apply(MET, 2, function(x) var(x, na.rm = TRUE)) > 0]
REP <- REP[, apply(REP, 2, function(x) var(x, na.rm = TRUE)) > 0]
PROT <- PROT[, apply(PROT, 2, function(x) var(x, na.rm = TRUE)) > 0]


# write
write.csv(shRNA, paste0(expression_dir, "/shrna.csv"))
write.csv(CNA, paste0(expression_dir, "/cna.csv"))
write.csv(MUT, paste0(expression_dir, "/mut.csv"))
write.csv(GE, paste0(expression_dir, "/ge.csv"))
write.csv(LIN, paste0(expression_dir, "/lin.csv"))
write.csv(MET, paste0(expression_dir, "/met.csv"))
write.csv(miRNA, paste0(expression_dir, "/mirna.csv"))
write.csv(XPR, paste0(expression_dir, "/xpr.csv"))
write.csv(REP, paste0(expression_dir, "/rep.csv"))
write.csv(PROT, paste0(expression_dir, "/prot.csv"))


# combined tables for RF and ENet
X.ccle <- GE %>% rownames_to_column(var = "ccle_name") %>%
  dplyr::inner_join(as.tibble(LIN, rownames = "ccle_name")) %>%
  dplyr::inner_join(MUT %>% rownames_to_column(var = "ccle_name")) %>%
  dplyr::inner_join(CNA %>% rownames_to_column(var = "ccle_name")) %>%
  column_to_rownames("ccle_name")
X.all <- X.ccle %>% rownames_to_column(var = "ccle_name") %>%
  dplyr::inner_join(XPR %>% rownames_to_column(var = "ccle_name")) %>%
  dplyr::inner_join(PROT %>% rownames_to_column(var = "ccle_name")) %>%
  dplyr::inner_join(miRNA %>% rownames_to_column(var = "ccle_name")) %>%
  dplyr::inner_join(MET %>% rownames_to_column(var = "ccle_name")) %>%
  column_to_rownames("ccle_name")

X.ccle <- X.ccle[, apply(X.ccle, 2, function(x) !any(is.na(x)))]
X.all <- X.all[, apply(X.all, 2, function(x) !any(is.na(x)))]

# write
write.csv(X.all, paste0(expression_dir, "/x-all.csv"))
write.csv(X.ccle, paste0(expression_dir, "/x-ccle.csv"))
