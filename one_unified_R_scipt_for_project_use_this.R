#### Installing Required Packages ####



# BiocManager::install("DESeq2")
# BiocManager::install("EnhancedVolcano")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("clusterProfiler")
# install.packages("stringr")
# BiocManager::install("BiocParallel")
# install.packages("dplyr")
# BiocManager::install("fgsea")
# install.packages("ggplot2")
# install.packages("cowplot")
# install.packages("colorspace")
# BiocManager::install("ComplexHeatmap")


#### Activating Required Packages ####

library(DESeq2)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
library(stats)
library(BiocParallel)
library(dplyr)
library(fgsea)
library(ggplot2)
library(cowplot)
library(colorspace)
library(ComplexHeatmap)

####1st script (For loading data, drawing volcano plot, and conducting GO enrichment analysis) ####

count.Data = read.delim("path/to/final_count_matrix.tsv", row.names=1)

## START - Editing a Count Matrix format ##

# count.Data is a data frame
# Update the row names by removing the part after the dot
rownames(count.Data) <- sub("\\..*", "", rownames(count.Data))
# Extract the gene_symbol column
gene_symbols <- count.Data[, 1]
# Create a new data frame with gene_id and gene_symbol concatenated
new_gene_id <- paste0(rownames(count.Data), "=", gene_symbols)
# Update the row names of count.Data to include the concatenated gene_id and gene_symbol
rownames(count.Data) <- new_gene_id
# Remove the first column
count.Data <- count.Data[, -1]
# Print the updated data frame
print(count.Data)

## End - Editing a Count Matrix format ##



## START - Dividing Data Type ##

count.Data  <- count.Data[rowSums(count.Data)>0,]

colData <- c(rep("EOLP",13),
             rep("NEOLP",27)) %>% as.factor() %>% as.data.frame()
rownames(colData) <- colnames(count.Data)
names(colData) <- "condition"
colData
EOLP.names <- rownames(colData)[1:13]
NEOLP.names <- rownames(colData)[14:40]
countData <- count.Data

## END - Dividing Data Type ##

## START - Running DESeq2 ##

dds = DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design=~condition)
dds = DESeq(dds)
plotDispEsts(dds)
norm.count <- counts(dds, normalized =T)

## END - Running DESeq2 ##

##----EOLP----##

res.EOLP <- results(dds, contrast=c("condition","EOLP","NEOLP"))
res.EOLP <- res.EOLP[order(res.EOLP$padj) ,]
res.EOLP <- res.EOLP[!is.na(res.EOLP$padj) ,]
summary(results(dds,alpha = 0.05))
label_genes_up <- rownames(res.EOLP)[res.EOLP$log2FoldChange > 2 &
                                       res.EOLP$padj <0.05]
label_genes_down <- rownames(res.EOLP)[res.EOLP$log2FoldChange < -2 &
                                         res.EOLP$padj <0.05]

# Get the row names from the data frame
row_names <- rownames(res.EOLP)

# Check each row name for the number of "=" characters
more_than_one_equal <- sapply(row_names, function(name) {
  return(sum(strsplit(name, NULL)[[1]] == "=") > 1)
})

# Extract row names with more than one "="
problematic_names <- row_names[more_than_one_equal]

# Print the results
if (length(problematic_names) == 0) {
  print("No row names have more than one '='.")
} else {
  print("Row names with more than one '=':")
  print(problematic_names)
}

## END - Running DESeq2 ##

## START - Volcano Plot for EOLP vs. NEOLP; Labelled for padj <0.05 & LFC > abs.(2)##
gene_labels <- sapply(row_names, function(name) {
  parts <- strsplit(name, "=")[[1]]
  return(parts[length(parts)])  # Use the last part after "="
})

# Extract gene names from row names
gene_labels <- sapply(rownames(res.EOLP), function(name) {
  parts <- strsplit(name, "=")[[1]]
  return(parts[length(parts)])
})

# Filter for significant genes
significant_genes <- res.EOLP$padj < 0.05 & abs(res.EOLP$log2FoldChange) > 5.0

# Extract labels only for significant genes
significant_labels <- gene_labels[significant_genes]

# Plot using EnhancedVolcano with labels for significant genes only
p <- EnhancedVolcano(
  res.EOLP,
  lab = gene_labels,  # Use all gene names, but only significant ones will be labeled
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.05,
  FCcutoff = 1.0,
  selectLab = significant_labels,  # Only significant genes are labeled
  legendLabels = c('Not sig.', 'Log2 FC', 'padj', 'padj & Log2 FC'),
  title = "EOLP vs. NEOLP\n\nLabeled Genes with Adjusted p-value < 0.05 and |Log2FoldChange| > 5",
  subtitle = NULL,
  cutoffLineType = "blank",
  pointSize = 0.8,
  labSize = 7.0,  # Set an appropriate size for readability
  drawConnectors = TRUE,
  typeConnectors = "open",
  maxoverlapsConnectors = Inf,  # Allow infinite overlaps
  widthConnectors = 0.5,
  colConnectors = 'blue',
  boxedLabels = FALSE,
  colAlpha = 4/5,
  arrowheads = TRUE,
  col = c('grey50', 'blue', 'purple', 'red')
)





# Center the title using theme()
p + theme(
  plot.title = element_text(hjust = 0.5)  # hjust = 0.5 centers the title
)

## END - Volcano Plot for EOLP vs. NEOLP; Labelled for padj <0.05 & LFC > abs.(2) ##


##START - GO Enrichment Analysis##

all_genes <- as.character(rownames(res.EOLP))
all_genes <- sub("=.+$", "", all_genes)

# Extract significant results
signif_res_up <- res.EOLP[res.EOLP$padj < 0.05 & !is.na(res.EOLP$padj) & res.EOLP$log2FoldChange>0, ]
signif_res_down <- res.EOLP[res.EOLP$padj < 0.05 & !is.na(res.EOLP$padj) & res.EOLP$log2FoldChange< -0, ]

signif_genes_up <- as.character(rownames(signif_res_up))
signif_genes_down <- as.character(rownames(signif_res_down))

# Extract the part before the '='
signif_genes_up <- sub("=.+$", "", signif_genes_up)

signif_genes_down <- sub("=.+$", "", signif_genes_down)

# Run GO enrichment analysis

ego_up <- enrichGO(gene = signif_genes_up,
                   universe = all_genes,
                   keyType = "ENSEMBL",
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05
                   ,readable = TRUE
)

ego_down <- enrichGO(gene = signif_genes_down,
                     universe = all_genes,
                     keyType = "ENSEMBL",
                     OrgDb = org.Hs.eg.db,
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.05
                     ,readable = TRUE)

cluster_high <- as.data.frame(ego_up)
cluster_low <- as.data.frame(ego_down)


## START - To sort cluster_high and cluster_low data frames by both p.adjust and GeneRatio in descending order of GeneRatio (and ensuring only significant results with p.adjust < 0.05) ##

# Ensure GeneRatio is treated correctly (convert it from "x/y" format to numeric)

convert_ratio <- function(ratio) {
  sapply(strsplit(ratio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}

# Add a numeric GeneRatio column to cluster_high and cluster_low
cluster_high$GeneRatio_numeric <- convert_ratio(cluster_high$GeneRatio)
cluster_low$GeneRatio_numeric <- convert_ratio(cluster_low$GeneRatio)

# Filter and sort cluster_high
cluster_high_sorted <- cluster_high[cluster_high$p.adjust < 0.05, ] %>%
  dplyr::arrange(desc(GeneRatio_numeric), p.adjust)

# Filter and sort cluster_low
cluster_low_sorted <- cluster_low[cluster_low$p.adjust < 0.05, ] %>%
  dplyr::arrange(desc(GeneRatio_numeric), p.adjust)

# View the sorted data frames
print(head(cluster_high_sorted))
print(head(cluster_low_sorted))

## END - To sort cluster_high and cluster_low data frames by both p.adjust and GeneRatio in descending order of GeneRatio (and ensuring only significant results with p.adjust < 0.05) ##

cluster_high_sorted$Condition <- rep("EOLP", nrow(cluster_high_sorted))
cluster_low_sorted$Condition <- rep("NEOLP", nrow(cluster_low_sorted))


cluster_high_sorted[1:20,]$Description
cluster_low_sorted[1:20,]$Description

head(cluster_high_sorted)
head(cluster_low_sorted)

cluster_high_sorted_bp <- cluster_high_sorted[cluster_high_sorted$ONTOLOGY == "BP", ]
cluster_low_sorted_bp <- cluster_low_sorted[cluster_low_sorted$ONTOLOGY == "BP", ]



cluster <- rbind(cluster_high_sorted_bp[1:7,], cluster_low_sorted_bp[1:7,])

cluster$GeneRatio <- sapply(cluster$GeneRatio, function(x) eval(parse(text=x)))

cluster <- cluster %>% 
  arrange(desc(Condition), GeneRatio, desc(p.adjust))

cluster$idx <- seq(dim(cluster)[1])
cluster.table <- cluster %>% dplyr::select(ONTOLOGY, ID, Description, GeneRatio, pvalue, p.adjust, qvalue, geneID)

## START - Insert line breaks into long pathway descriptions and add new  ##

cluster$Description_wrapped <- sapply(cluster$Description, function(x) {
  if (nchar(x) > 50) {
    paste(strwrap(x, width = 50), collapse = "\n")
  } else {
    x
  }
})

## END - Insert line breaks into long pathway descriptions and add new ##


ggplot(cluster)+
  geom_point(aes(x = Condition,
                 y = reorder(Description_wrapped,idx),
                 color = p.adjust, size = GeneRatio))+
  ggtitle("EOLP vs. NEOLP, GO Enrichment Analysis") +
  scale_color_gradient(low = "darkgreen", high="lightgray") +
  guides(color = guide_colorbar(order=1), size = guide_legend(order=2)) +
  ylab("Pathways") + xlab("") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 50, vjust = 1, hjust=1,
                               size=20),
    axis.text.y = element_text(size=21, family = "serif", margin = margin(t = 10, b = 10)),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5, margin = margin(b = 10)))





#### 2nd Part (For Running gene set enrichment analysis (GSEA) and visualize the result using lollipop plot) ####


## START - Drawing a lollipop plot based on fgsea result ##

res <- res.EOLP
res
summary(res)

res.df <- as.data.frame(res)


res.df <- res.df %>% dplyr::mutate(ranks = (log2FoldChange*(-log10(padj))))

# res.df is a dataframe

# Extract the gene symbols from the row names
res.df$gene_symbol <- sub(".*=", "", rownames(res.df))

# Print the dataframe to verify the changes
print(res.df)

# Find row names with more than one '=' sign
rows_with_multiple_equals <- grep(".*=.*=.*", rownames(res.df), value = TRUE)

# Print the result
if (length(rows_with_multiple_equals) > 0) {
  cat("Row names with more than one '=' sign:\n")
  print(rows_with_multiple_equals)
} else {
  cat("No row names with more than one '=' sign found.\n")
}

lfc_vector <- res.df$ranks

lfc_vector[is.na(lfc_vector)]

lfc_vector <- as.numeric(lfc_vector)

names(lfc_vector) <- res.df$gene_symbol

head(lfc_vector)

any(is.na(lfc_vector))

lfc_vector <- sort(lfc_vector, decreasing = TRUE)

head(lfc_vector)

hist(lfc_vector)

length(lfc_vector)

f1 = read.csv("path/to/inflammation_reference/for_lollipop/Inflammation_Genes_List.csv")

colnames(f1)[1] <- "MitoPathway"
colnames(f1)[2] <- "MitoPathway.Hierarchy"
colnames(f1)[3] <- "Genes"
gg <- as.data.frame(f1[,3])
sym <- as.data.frame(f1[,1])
head(f1)
nrow(f1)
g1 <- NULL
n1 <- NULL
new_list <- NULL

for (i in c(1:37)){
  g1 <- as.vector(gg[i,])
  g1 <- str_split(g1,", ")
  new_list <- append(new_list,g1)
}

head(new_list)

a <- NULL
for (i in c(1:37)){
  a <- c(a,sym[i,])
}

class(a)
length(a)

names(new_list) <- a
length(new_list)
head(new_list)
head(lfc_vector)



set.seed(42)

## START - finding elements with the same gene symbol #

# Find duplicated elements
duplicates <- duplicated(lfc_vector)

# Print duplicates vector
print(duplicates)

# Find all duplicated elements
duplicated_elements <- unique(lfc_vector[duplicated(lfc_vector)])

# Print duplicated elements
if (length(duplicated_elements) == 0) {
  print("No duplicates found.")
} else {
  print("Duplicated elements:")
  print(duplicated_elements)
}

##END - finding elements with the same gene symbol ##

##START - finding elements with the same rank value ##

# Filter rows where 'ranks' is 0
rows_with_zero_rank <- res.df %>%
  filter(ranks == 0) %>%
  rownames()

# Print the row names
if (length(rows_with_zero_rank) == 0) {
  print("No rows with rank value of 0 found.")
} else {
  print("Row names with rank value of 0:")
  print(rows_with_zero_rank)
}



# Assuming 'res.df' has row names set, ensure they're accessible
res.df <- res.df %>%
  mutate(row_names = rownames(res.df))  # Add row names as a new column

# Group by the 'ranks' column, summarize the row names, and filter for duplicates
duplicate_ranks <- res.df %>%
  group_by(ranks) %>%
  filter(n() > 1) %>%  # Keep only groups with more than one member
  summarize(rows = paste(row_names, collapse = " = ")) %>%
  ungroup()

# Print the results
if (nrow(duplicate_ranks) == 0) {
  print("No rows have the same rank value.")
} else {
  duplicate_ranks$rows <- paste0(duplicate_ranks$rows, " = ", duplicate_ranks$ranks)
  print(duplicate_ranks$rows)
}

##END - finding elements with the same rank value##

##START - finding duplicated gene symbols (the one with different ensemble ID but the same gene symbol ##

# Find duplicated names
duplicated_names <- names(lfc_vector)[duplicated(names(lfc_vector))]

# Get unique duplicated names
unique_duplicated_names <- unique(duplicated_names)

# Print the duplicated names
if (length(unique_duplicated_names) == 0) {
  print("No duplicated names found.")
} else {
  print("Duplicated names:")
  print(unique_duplicated_names)
}

# Step 1: Unlist all genes from new_list
all_genes_in_new_list <- unlist(new_list)

# Step 2: Check for overlap
overlap_genes <- unique_duplicated_names[unique_duplicated_names %in% all_genes_in_new_list]

overlap_genes_test <- unique_duplicated_names[unique_duplicated_names %in% all_genes_in_new_list]


# Step 3: Print the results
if (length(overlap_genes) > 0) {
  print("Overlap between unique_duplicated_names and genes in new_list:")
  print(overlap_genes)
} else {
  print("No overlap found between unique_duplicated_names and genes in new_list.")
}


# Assuming unique_duplicated_names exists and adding KLK15 and NOX1
unique_duplicated_names_test <- c(unique_duplicated_names, "NOX1", "KLK15")

# Unlist all genes from new_list
all_genes_in_new_list <- unlist(new_list)

# Check for overlap
overlap_genes <- unique_duplicated_names_test[unique_duplicated_names_test %in% all_genes_in_new_list]

# Print the results
if (length(overlap_genes) > 0) {
  print("Overlap between unique_duplicated_names_test and genes in new_list:")
  print(overlap_genes)
} else {
  print("No overlap found between unique_duplicated_names_test and genes in new_list.")
}


## END - finding duplicated gene symbols (the one with different ensemble ID but the same gene symbol ##

## checked it works, and indeed there is no duplicated gene names between unique_duplicated_names and all_genes_in_new_list ##

# Extract the names from lfc_vector
gene_names <- names(lfc_vector)

# Initialize a new vector to store updated names with suffixes
updated_gene_names <- gene_names

# Use a named numeric vector to count occurrences
name_counts <- setNames(rep(0, length(gene_names)), gene_names)

# Loop through each gene name to assign suffixes to duplicates
for (i in seq_along(gene_names)) {
  current_name <- gene_names[i]
  name_counts[current_name] <- name_counts[current_name] + 1
  if (name_counts[current_name] > 1) {
    updated_gene_names[i] <- paste0(current_name, ".", name_counts[current_name])
  }
}

# Create a new vector with original values and updated names
lfc_vector_duplicated_solved <- setNames(lfc_vector, updated_gene_names)

# Display the modified vector
print(lfc_vector_duplicated_solved)

# Filtering to find genes starting with "MRPS31P5"
MRPS31P5_genes <- lfc_vector_duplicated_solved[grep("^MRPS31P5", names(lfc_vector_duplicated_solved))]
print("Genes starting with MRPS31P5:")
print(MRPS31P5_genes)

# Filtering to find genes starting with "RAET1E-AS1"
RAET1E_AS1_genes <- lfc_vector_duplicated_solved[grep("^RAET1E-AS1", names(lfc_vector_duplicated_solved))]
print("Genes starting with RAET1E-AS1:")
print(RAET1E_AS1_genes)

set.seed(42)
fgseaRes <- fgsea(pathways = new_list,
                  stats = lfc_vector_duplicated_solved,
                  minSize = 5,
                  gseaParam = 1
                  ,maxSize = 500,
                  nproc = 1,
                  nPermSimple = 10000
) 

fgseaRes <- fgseaRes[order(padj),]
nrow(fgseaRes)
fg_df <- as.data.frame(fgseaRes)
head(fg_df)
tail(fg_df)

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(new_list[topPathways], lfc_vector, fgseaRes,
              gseaParam = 1)

## adjusted p-value threshold setting ##
test_path <- as.data.frame(f1[,1:2])
names(test_path)[1] <- "pathway"
head(test_path)
test_merge <- merge(fg_df,test_path,by="pathway")
tail(test_merge)
any(duplicated(test_merge$pathway))

test_pval <- test_merge[test_merge$padj< 0.5,]
nrow(test_pval)

##pathways##


test_pval$MitoPathway.Hierarchy %>% table
# Innate Immune

a<- test_pval[test_pval$MitoPathway.Hierarchy =="Innate Immune",]
c <- a
nrow(c)

fg_central <- c %>%
  arrange(NES) %>%
  mutate(flag = ifelse(NES > 0, "EOLP", "NEOLP"),
         name = factor(pathway, levels =. $pathway))

if (all(c$NES < 0)) {
  color <- "Steel Blue"
} else if(all(c$NES > 0)){
  color <- "Light Coral"
}else{color <- c("Light Coral", "Steel Blue")}



p1 <- ggplot(fg_central, aes(x=NES, y=pathway, color=flag)) +
  geom_segment(aes(x=0, y=pathway, xend=NES, yend=pathway)) +
  geom_point(size=2) + # for default color (red/ blue)
  scale_colour_manual(values = color)+
  coord_equal(ratio = 0.5, xlim = c(-2.5, 2.5), ylim = c(1, nrow(c))) +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size= 15),
        axis.title = element_blank(),
        axis.text.y = element_text(size=5, lineheight = 0.5),
        axis.title.y = element_blank())+
  labs(title = "Innate Immune")
p1



# Downregulated
a<- test_pval[test_pval$MitoPathway.Hierarchy =="Downregulated",]
c <- a
nrow(c)

fg_pish <- c %>%
  arrange(NES) %>%
  mutate(flag = ifelse(NES > 0, "EOLP", "NEOLP
"),
         name = factor(pathway, levels =. $pathway))
if (all(c$NES < 0)) {
  color <- "Steel Blue"
} else if(all(c$NES > 0)){
  color <- "Light Coral"
}else{color <- c("Light Coral", "Steel Blue")}


p2<- ggplot(fg_pish, aes(x=NES, y=pathway, color=flag)) +
  geom_segment(aes(x=0, y=pathway, xend=NES, yend=pathway)) +
  scale_colour_manual(values = color)+
  geom_point(size=2) + # for default color (red/ blue) 
  labs(title = "Downregulated") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size= 15),
        axis.title = element_blank(),
        axis.text.y = element_text(size=8, lineheight = 0.5),
        axis.title.y = element_blank())+
  coord_equal(ratio = 0.5, xlim = c(-2.5, 2.5), ylim = c(1, nrow(c))) 

p2

# Mitochondrial Innate Immune


a<- test_pval[test_pval$MitoPathway.Hierarchy =="Mitochondrial Innate Immune",]
c <- a
head(c)
nrow(c)

fg_oxphos <- c %>%
  arrange(NES) %>%
  mutate(flag = ifelse(NES > 0, "EOLP", "NEOLP"),
         name = factor(pathway, levels =. $pathway))
if (all(c$NES < 0)) {
  color <- "Steel Blue"
} else if(all(c$NES > 0)){
  color <- "Light Coral"
}else{color <- c("Light Coral", "Steel Blue")}


p3 <- ggplot(fg_oxphos, aes(x=NES, y=pathway, color=flag)) +
  geom_segment(aes(x=0, y=pathway, xend=NES, yend=pathway)) +
  geom_point(size=2) + # for default color (red/ blue) 
  scale_colour_manual(values = color)+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size= 15),
        axis.title = element_blank(),
        axis.text.y = element_text(size=8, lineheight = 0.5),
        axis.title.y = element_blank())+
  coord_equal(ratio = 0.5, xlim = c(-2.5, 2.5), ylim = c(1, nrow(c))) + 
  labs(title = "Mitochondrial Innate Immune") 
p3


# Integrated Stress Response (ISR)


a<- test_pval[test_pval$MitoPathway.Hierarchy =="Integrated Stress Response (ISR)",]
c <- a
nrow(c)

fg_smt <- c %>%
  arrange(NES) %>%
  mutate(flag = ifelse(NES > 0, "EOLP", "NEOLP"),
         name = factor(pathway, levels =. $pathway))
if (all(c$NES < 0)) {
  color <- "Steel Blue"
} else if(all(c$NES > 0)){
  color <- "Light Coral"
}else{color <- c("Light Coral", "Steel Blue")}


p4<-ggplot(fg_smt, aes(x=NES, y=pathway, color=flag)) +
  geom_segment(aes(x=0, y=pathway, xend=NES, yend=pathway)) +
  geom_point(size=2) + # for default color (red/ blue) 
  scale_colour_manual(values = color)+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size= 15),
        axis.title = element_blank(),
        axis.text.y = element_text(size=8, lineheight = 0.5),
        axis.title.y = element_blank())+
  coord_equal(ratio = 0.5, xlim = c(-2.5, 2.5), ylim = c(1, nrow(c))) + 
  labs(title = "Integrated Stress Response (ISR)") 
p4


# Unfolded Protein Response (UPR)
a<- test_pval[test_pval$MitoPathway.Hierarchy =="Unfolded Protein Response (UPR)",]
c <- a
nrow(c)

fg_signaling <- c %>%
  arrange(NES) %>%
  mutate(flag = ifelse(NES > 0, "EOLP", "NEOLP"),
         name = factor(pathway, levels =. $pathway))
if (all(c$NES < 0)) {
  color <- "Steel Blue"
} else if(all(c$NES > 0)){
  color <- "Light Coral"
}else{color <- c("Light Coral", "Steel Blue")}


p5 <- ggplot(fg_signaling, aes(x=NES, y=pathway, color=flag)) +
  geom_segment(aes(x=0, y=pathway, xend=NES, yend=pathway)) +
  geom_point(size=2) + # for default color (red/ blue) 
  scale_colour_manual(values = color)+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size= 15),
        axis.title = element_blank(),
        axis.text.y = element_text(size=8, lineheight = 0.5),
        axis.title.y = element_blank())+
  coord_equal(ratio = 0.5, xlim = c(-2.5, 2.5), ylim = c(1, nrow(c))) + 
  labs(title = "Unfolded Protein Response (UPR)") 
p5

# Adaptive/Extracellular Mediated Immunity
a<- test_pval[test_pval$MitoPathway.Hierarchy =="Adaptive/Extracellular Mediated Immunity",]
c <- a
nrow(c)

fg_mds <- c %>%
  arrange(NES) %>%
  mutate(flag = ifelse(NES > 0, "EOLP", "NEOLP"),
         name = factor(pathway, levels =. $pathway))
if (all(c$NES < 0)) {
  color <- "Steel Blue"
} else if(all(c$NES > 0)){
  color <- "Light Coral"
}else{color <- c("Light Coral", "Steel Blue")}


p6 <- ggplot(fg_mds, aes(x=NES, y=pathway, color=flag)) +
  geom_segment(aes(x=0, y=pathway, xend=NES, yend=pathway)) +
  geom_point(size=2) + # for default color (red/ blue) 
  scale_colour_manual(values = color)+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size= 15),
        axis.title = element_blank(),
        axis.text.y = element_text(size=8, lineheight = 0.5),
        axis.title.y = element_blank())+
  coord_equal(ratio = 0.5, xlim = c(-2.5, 2.5), ylim = c(1, nrow(c))) + 
  labs(title = "Adaptive/Extracellular Mediated Immunity")
p6
# Renin-Angiotensin-Aldosterone System (RAAS)

a<- test_pval[test_pval$MitoPathway.Hierarchy =="Renin-Angiotensin-Aldosterone System (RAAS)",]
c <- a
head(c)
nrow(c)



fg_metabolism <- c %>%
  arrange(NES) %>%
  mutate(flag = ifelse(NES > 0, "EOLP", "NEOLP"),
         name = factor(pathway, levels =. $pathway))

if (all(c$NES < 0)) {
  color <- "Steel Blue"
} else if(all(c$NES > 0)){
  color <- "Light Coral"
}else{color <- c("Light Coral", "Steel Blue")}

p7 <- ggplot(fg_metabolism, aes(x=NES, y=pathway, color=flag)) +
  geom_segment(aes(x=0, y=pathway, xend=NES, yend=pathway)) +
  geom_point(size=2) + # for default color (red/ blue) 
  scale_colour_manual(values = color)+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size= 15),
        axis.title = element_blank(),
        axis.text.y = element_text(size=8, lineheight = 0.5),
        axis.title.y = element_blank())+
  coord_equal(ratio = 0.5, xlim = c(-2.5, 2.5), ylim = c(1, nrow(c))) + 
  labs(title = "Renin-Angiotensin-Aldosterone System (RAAS)")
p7

plot_grid(p1, p2, p3, p4, p5, p6, p7, rel_heights = c(1, 1, 1))



##### 3rd Script (For drawing a heat map plot based on z-score (only for genes in "Adaptive_Extracellular_Mediated_Immunity_Genes.csv" data set)) #####



## getting the reference for the heat map ##

ref <- read.csv(file = "path/to/Adaptive_Extracellular_Mediated_Immunity_Genes.csv")

### You can use other datasets to explore gene expression in the selected dataset. Please read the README.md carefully for detailed instructions—especially regarding RAAS_Genes.csv, which includes NOX5 (a gene with two Ensembl IDs) ###

# ref <- read.csv(file = "path/to/inflammation_reference/for_heatmap/Innate_Immune_Genes.csv")
# 
# ref <- read.csv(file = "path/to/inflammation_reference/for_heatmap/ISR_Genes.csv")
# 
# ref <- read.csv(file = "path/to/inflammation_reference/for_heatmap/Mitochondrial_Innate_Immune_Genes.csv")
# 
# ref <- read.csv(file = "path/to/inflammation_reference/for_heatmap/RAAS_Genes.csv")
# 
# ref <- read.csv(file = "path/to/inflammation_reference/for_heatmap/UPR_Genes.csv")




list_a <- norm.count %>% as.data.frame()
list_a %>% head
ttt <- list_a %>% t
anno <- colData %>% subset(select= condition)
tttt <- merge(ttt, anno, by =0)


tttt %>% group_by(condition)

rownames(tttt) <- tttt$Row.names
tttt <- tttt[,-1]


## START - choosing genes of interest ##

# Check if the last column is 'condition' and remove it if true

if (colnames(tttt)[ncol(tttt)] == "condition") {
  ttttt <- tttt[, -ncol(tttt)]
}

# Add row names to a new column called 'Group.1'
ttttt$Group.1 <- rownames(tttt)

# Move 'Group.1' to be the first column
ttttt <- ttttt[, c(ncol(tttt), 1:(ncol(tttt) - 1))]

rownames(ttttt) <- seq_len(nrow(ttttt))  # This explicitly sets the row names to 1, 2, 3, ...


## END - choosing genes of interest ##

tail(ttttt)

cc <- ttttt %>% t

cc %>% head
colnames(cc) <- cc[1,]
cc %>% head
cc<- cc[-1,]

cc.db <- cc  %>% as.data.frame() %>% sapply("as.double") %>% as.data.frame()
rownames(cc.db) <- rownames(cc)


## START - Examine Duplicated Gene Symbol (The one with the same gene symbol but different Ensembl ID) ##

# Data frame is named cc.db and ref

# Step 1: Extract gene names from cc.db row names
current_row_names <- rownames(cc.db)
gene_names <- sub(".*=", "", current_row_names)

# Step 2: Identify duplicated gene names
duplicated_gene_names <- gene_names[duplicated(gene_names) | duplicated(gene_names, fromLast = TRUE)]
unique_duplicated_gene_names <- unique(duplicated_gene_names)

cat("Duplicated gene names in cc.db:\n")
print(unique_duplicated_gene_names)

# Step 3: Check for common genes between duplicated genes and ref's Approved.symbol
common_genes <- intersect(unique_duplicated_gene_names, ref$SYM)

common_genes_test <- intersect("CD74", ref$SYM)

cat("\nCommon genes between duplicated genes and ref's Approved.symbol:\n")
print(common_genes)
print(common_genes_test)

## END - Examine Duplicated Gene Symbol (The one with the same gene symbol but different Ensemble ID)##

cc.db$sym <- rownames(cc.db)

## START - replacing column names ##

# Replace "EOLP" and "NEOLP" with actual column names of cc.db
cc.db.tmp <- cc.db[, !colnames(cc.db) %in% "sym"]

# Check if the number of columns in cc.db.tmp is 40
if (ncol(cc.db.tmp) == 40) {
  print("The number of columns in cc.db.tmp is 40.")
} else {
  print(paste("The number of columns in cc.db.tmp is", ncol(cc.db), "instead of 40."))
}


new_ref <- c("SYM", "pathway")

ref <- ref[,new_ref]

ref <- ref %>% rename(sym = SYM)
ref <- ref[!is.na(ref$sym) ,]

## END -replacing column names ##


## START - removing duplicated entries in the heatmap ##

# Data frame is named 'ref'
ref_cleaned <- ref %>%
  distinct(sym, pathway, .keep_all = TRUE)  # Step 1: Remove duplicates based on 'sym' and 'pathway'

# Save the result back into ref
ref <- ref_cleaned

rm(ref_cleaned)


## END - removing duplicated entries in the heatmap ##



# Dataframes are named ref and cc.db

# Step 1: Extract gene names from the sym column in cc.db
cc.db$gene_name <- sub(".*=", "", cc.db$sym)



# Step 2: Merge ref and cc.db on the extracted gene name
list_ref <- merge(ref, cc.db, by.x = "sym", by.y = "gene_name")

# Print the merged dataframe to verify the merge
print(list_ref)

list_ref %>% dim

list_ref <- list_ref[!duplicated(list_ref),]

# pathway list

#remove all the duplicated pathways
pathways <- unique(list_ref$pathway)
list_ref %>% head

## START - for selecting 13 columns in NEOLP (there are 13 EOLP samples, and to compare z-score between the two groups, chose the same number (13) of samples from NEOLP group ##

# Set seed for reproducibility
set.seed(42)

# Select the columns
selected_columns <- c(
  colnames(list_ref)[3:15],  # Select columns from 3rd to 15th
  sample(colnames(list_ref)[16:42], 13)  # Randomly select 13 columns from 16th to 42nd
)

# Save the selected column names into 'd'
d <- selected_columns


## END - for selecting 13 columns in NEOLP (there are 13 EOLP samples, and to compare z-score between the two groups, chose the same number (13) of samples from NEOLP group ##

pathways %>% length

Surface_Marker_Receptor_Signaling <- list_ref[list_ref$pathway ==
                                                "Surface Marker/Receptor Signaling",]

Interleukins <- list_ref[list_ref$pathway == 
                           "Interleukins",] 

Cytokines <- list_ref[list_ref$pathway ==
                        "Cytokines",]

Antigen_Presentation <- list_ref[list_ref$pathway ==
                                   "Antigen Presentation",]

rename.v <- c("Surface Marker/Receptor Signaling" = "Surface Marker\n/Receptor Signaling",
              "Interleukins" = "Interleukins",
              "Cytokines" = "Cytokines",
              "Antigen Presentation" = "Antigen\nPresentation"
)

Antigen_Presentation$sym


nrow(Surface_Marker_Receptor_Signaling)
nrow(Interleukins)
nrow(Cytokines)
nrow(Antigen_Presentation)

all_list <- list(Surface_Marker_Receptor_Signaling,
                 Interleukins,
                 Cytokines,
                 Antigen_Presentation
)



# repeated code

# Subset the data frame to include only the selected columns
dt_file <- list_ref[, selected_columns]

sym <- as.vector(list_ref$sym)

dt_file <- as.matrix(dt_file)
rownames(dt_file) <- sym
dt_file <- t(dt_file)

all <- colnames(countData)
all_rev <- colData$condition

aka5 = data.frame(Samples = all_rev)

rownames(aka5) <- all
rownames(dt_file)

## START - Examining row names  ##

dt_file_joo_test <- dt_file[rownames(dt_file),]

## END - Examining row names ##


dt_file <- dt_file[rownames(dt_file),]

rm(dt_file_joo_test)


hcl_palettes(plot = TRUE)

#color
path.n <- c()
group_sym <- c()
for (i in all_list){
  tmp <- rep(i$pathway[1], nrow(i))
  group_sym <- c(group_sym, i$sym)
  path.n <- c(path.n, tmp)
}

tmp
head(all_list)
group_sym
path.n
path.n <- path.n %>% as.factor()
sym.n <- levels(path.n)
ncolors <- rainbow(length(sym.n),  s = 0.7, v = 1, alpha = 0.8)
pie(rep(1,length(sym.n)), col = ncolors)
names(ncolors) <- sym.n
aka4 = data.frame(Pathway=path.n)
aka4
aka3 = list(Pathway=ncolors)
aka3
column_split <- path.n

list_1 <- list(Surface_Marker_Receptor_Signaling,
               Interleukins,
               Cytokines,
               Antigen_Presentation
               
)

# repeated code
dt_file <- list_ref[, selected_columns]
sym <- as.vector(list_ref$sym)


dt_file <- as.matrix(dt_file)
rownames(dt_file) <- sym
dt_file <- t(dt_file)

all <- colnames(countData)
all_rev <- colData$condition

# aka5 is not used 
# aka5 = data.frame(Samples = rownames(dt_file))

rownames(dt_file)
dt_file <- dt_file[rownames(dt_file),]

hcl_palettes(plot = TRUE)

# group1 

path.n <- c()
group_sym <- c()
for (i in list_1){
  tmp <- rep(i$pathway[1], nrow(i))
  group_sym <- c(group_sym, i$sym)
  path.n <- c(path.n, tmp)
}
path.n <- path.n %>% as.factor()

sym.n <- levels(path.n)

ncolors <- rainbow(length(sym.n),  s = 0.7, v = 1, alpha = 0.8)


names(ncolors) <- sym.n

aka4 = data.frame(Pathway=path.n)
aka4 <- as.matrix(aka4)
rownames(aka4) <- group_sym


# aka3 is not used 

column_split <- path.n

dt_file <- dt_file[, group_sym] %>% as.matrix()

dt_file <- dt_file[rownames(dt_file), ]


#replace the name on column_split with the values in remane.v
column_split_t <- plyr::revalue(column_split, rename.v)

str(dt_file)

# Extract gene symbols from the row names of res.EOLP
res.EOLP$gene_symbol <- str_extract(rownames(res.EOLP), "(?<=\\=).*")

# Filter the DESeqResults to keep only the significant genes with the desired log2FoldChange
significant_genes <- res.EOLP[res.EOLP$padj < 0.05 & 
                                (res.EOLP$log2FoldChange > 1 | res.EOLP$log2FoldChange < -1), 
                              "gene_symbol"]

# Subset dt_file to keep only the columns with matching gene symbols
dt_file_filtered <- dt_file[, colnames(dt_file) %in% significant_genes, drop = FALSE]


#Step 1: Get the gene symbols from dt_file_filtered
filtered_genes <- colnames(dt_file_filtered)

# Step 2: Filter list_1 to keep only entries with genes present in dt_file_filtered
filtered_list_1 <- lapply(list_1, function(df) {
  df[df$sym %in% filtered_genes, ]
})

# Step 3: Remove any empty elements from filtered_list_1 (optional)
filtered_list_1 <- Filter(function(x) nrow(x) > 0, filtered_list_1)

# Step 4: Verify that the filtering worked
cat("Number of filtered genes:", length(filtered_genes), "\n")
cat("Number of entries in filtered_list_1:", length(filtered_list_1), "\n")

path.n <- c()
group_sym <- c()
for (i in filtered_list_1){
  tmp <- rep(i$pathway[1], nrow(i))
  group_sym <- c(group_sym, i$sym)
  path.n <- c(path.n, tmp)
}
path.n <- path.n %>% as.factor()

sym.n <- levels(path.n)

ncolors <- rainbow(length(sym.n),  s = 0.7, v = 1, alpha = 0.8)


names(ncolors) <- sym.n

aka4 = data.frame(Pathway=path.n)
aka4 <- as.matrix(aka4)
rownames(aka4) <- group_sym

column_split <- path.n

#replace the name on column_split with the values in remane.v
column_split_t <- plyr::revalue(column_split, rename.v)


# Step 5: Redraw the heatmap with the adjusted 'column_split'

ComplexHeatmap::pheatmap(
  dt_file_filtered,
  scale = "column",
  treeheight_row = 0,
  treeheight_col = 20,
  cellwidth = 50,
  cellheight = 18,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  annotation_legend = TRUE,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  border_color = "#BDBDBD",
  main = "",
  name = "z-score for\nAdaptive/ \nExtracellular \nMediated \nImmunity \nGenes", 
  legend = TRUE,
  fontsize = 11,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  column_split = column_split_t,  # Use the filtered column split
  column_gap = unit(1, "mm"),
)

# This step is added to visualize multiple heatmaps if necessary by saving a generated heatmap as a variable
p2 <- ComplexHeatmap::pheatmap(
  dt_file_filtered,
  scale = "column",
  treeheight_row = 0,
  treeheight_col = 20,
  cellwidth = 50,
  cellheight = 18,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  annotation_legend = TRUE,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  border_color = "#BDBDBD",
  main = "",
  name = "z-score for\nAdaptive/ \nExtracellular \nMediated \nImmunity \nGenes", 
  legend = TRUE,
  fontsize = 11,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  column_split = column_split_t,  # Use the filtered column split
  column_gap = unit(1, "mm"),
)


# finish editing for github release (2025.04.03) #