##### 3rd Script (For drawing a heat map plot based on z-score) (only for genes in "RAAS_Genes.csv" data set) #####

ref <- read.csv(file = "path/to/inflammation_reference/for_heatmap/RAAS_Genes.csv")


list_a <- norm.count %>% as.data.frame()
list_a %>% head
ttt <- list_a %>% t
anno <- colData %>% subset(select= condition)
tttt <- merge(ttt, anno, by =0)


tttt %>% group_by(condition)






rownames(tttt) <- tttt$Row.names
tttt <- tttt[,-1]


# Check if the last column is 'condition' and remove it if true
if (colnames(tttt)[ncol(tttt)] == "condition") {
  ttttt <- tttt[, -ncol(tttt)]
}

# Add row names to a new column called 'Group.1'
ttttt$Group.1 <- rownames(tttt)

# Move 'Group.1' to be the first column
ttttt <- ttttt[, c(ncol(tttt), 1:(ncol(tttt) - 1))]

rownames(ttttt) <- seq_len(nrow(ttttt))  # This explicitly sets the row names to 1, 2, 3, ...

tail(ttttt)

cc <- ttttt %>% t

cc %>% head
colnames(cc) <- cc[1,]
cc %>% head
cc<- cc[-1,]



cc.db <- cc  %>% as.data.frame() %>% sapply("as.double") %>% as.data.frame()
rownames(cc.db) <- rownames(cc)


# cc.db and ref are dataframes

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


## END - replacing column names ##




# list_ref


new_ref <- c("SYM", "pathway")

ref <- ref[,new_ref]

ref <- ref %>% rename(sym = SYM)
ref <- ref[!is.na(ref$sym) ,]


# ref is data frame

ref_cleaned <- ref %>%
  distinct(sym, pathway, .keep_all = TRUE)  # Step 1: Remove duplicates based on 'sym' and 'pathway'

# Save the result back into ref
ref <- ref_cleaned

rm(ref_cleaned)


## START - solving issues regarding (NOX5 = two Ensemble IDs) ##

# Create a modified version of the ref data frame
ref_modified <- ref %>%
  # Find rows that match NOX5 using dplyr's filter function explicitly
  dplyr::filter(sym == "NOX5") %>%
  # Replicate the row to create two new rows
  slice(rep(1, 2)) %>%
  # Assign new sym values to each replicated row with the required format
  mutate(sym = c("NOX5_ENSG00000255346", "NOX5_ENSG00000290203"))

# Combine the modified rows back into the original data frame
ref <- ref %>%
  # Remove the original NOX5 row
  dplyr::filter(sym != "NOX5") %>%
  # Add the new rows
  bind_rows(ref_modified)

# Display the updated ref data frame
print(ref)


## END - solving issues regarding (NOX5 = two Ensemble IDs) ##




# ref and cc.db are dataframes

# Step 1: Extract gene names from the sym column in cc.db
cc.db$gene_name <- sub(".*=", "", cc.db$sym)


# Step 2: Merge ref and cc.db on the extracted gene name

## START - change cc.db NOX5 gene name according to ensembl name ##

# Update the gene_name column based on conditions in the sym column
cc.db <- cc.db %>%
  mutate(gene_name = case_when(
    gene_name == "NOX5" & sym == "ENSG00000255346=NOX5" ~ "NOX5_ENSG00000255346",
    gene_name == "NOX5" & sym == "ENSG00000290203=NOX5" ~ "NOX5_ENSG00000290203",
    TRUE ~ gene_name  # Keep the original gene name for other cases
  ))

# Display the modified cc.db data frame
print(cc.db)

## END - change cc.db NOX5 gene name according to ensembl name ##

list_ref <- merge(ref, cc.db, by.x = "sym", by.y = "gene_name")

# Print the merged dataframe to verify the merge
print(list_ref)

list_ref %>% dim

list_ref <- list_ref[!duplicated(list_ref),]




# pathway list 

#remove all the duplicated pathways
pathways <- unique(list_ref$pathway)
list_ref %>% head

# START - for selecting 13 columns in NEOLP ##

# Set seed for reproducibility
set.seed(42)

# Select the columns
selected_columns <- c(
  colnames(list_ref)[3:15],  # Select columns from 3rd to 15th
  sample(colnames(list_ref)[16:42], 13)  # Randomly select 13 columns from 16th to 42nd
)

# Save the selected column names into 'd'
d <- selected_columns


## END - for selecting 13 columns in NEOLP ##

pathways %>% length

AGT_Regulator_Axis <- list_ref[list_ref$pathway ==
                            "AGT Regulator Axis",]

NADPH_Oxidase <- list_ref[list_ref$pathway == 
                         "NADPH Oxidase",] 

PANoptosis <- list_ref[list_ref$pathway ==
                            "PANoptosis",]

Complement_activation_Fibrin_deposition <- list_ref[list_ref$pathway ==
                            "Complement activation / Fibrin deposition",]


Syndecans <- list_ref[list_ref$pathway ==
                         "Syndecans",]

Hyaluronan_Accumulation <- list_ref[list_ref$pathway ==
                         "Hyaluronan Accumulation",]

Bradykinin_Production <- list_ref[list_ref$pathway ==
                         "Bradykinin Production",]




rename.v <- c("AGT Regulator Axis" = "AGT\nRegulator\nAxis",
              "NADPH Oxidase" = "NADPH\nOxidase",
              "PANoptosis" = "PANoptosis",
              "Complement activation / Fibrin deposition" = "Complement\nactivation / \nFibrin deposition",
              "Syndecans" = "Syndecans",
              "Hyaluronan Accumulation" = "Hyaluronan\nAccumulation",
              "Bradykinin Production" = "Bradykinin\nProduction"
)



AGT_Regulator_Axis$sym


nrow(AGT_Regulator_Axis)
nrow(NADPH_Oxidase)
nrow(PANoptosis)
nrow(Complement_activation_Fibrin_deposition)
nrow(Syndecans)
nrow(Hyaluronan_Accumulation)
nrow(Bradykinin_Production)




all_list <- list(AGT_Regulator_Axis,
                 NADPH_Oxidase,
                 PANoptosis,
                 Complement_activation_Fibrin_deposition,
                 Syndecans,
                 Hyaluronan_Accumulation,
                 Bradykinin_Production
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




## START - row names applied test ##

dt_file_joo_test <- dt_file[rownames(dt_file),]

## END - row names applied test ##


dt_file <- dt_file[rownames(dt_file),]

rm(dt_file_joo_test)


hcl_palettes(plot = TRUE)


# color
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


list_1 <- list(AGT_Regulator_Axis,
               NADPH_Oxidase,
               PANoptosis,
               Complement_activation_Fibrin_deposition,
               Syndecans,
               Hyaluronan_Accumulation,
               Bradykinin_Production
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

## START - try to change NOX5 into ensembl ID ##

# Initialize empty vectors
path.n <- c()
group_sym <- c()


for (i in list_1){
  tmp <- rep(i$pathway[1], nrow(i))
  group_sym <- c(group_sym, i$sym)
  path.n <- c(path.n, tmp)
}
path.n <- path.n %>% as.factor()

sym.n <- levels(path.n)



# Convert path.n to a factor
path.n <- as.factor(path.n)

# Get the unique levels of path.n
sym.n <- levels(path.n)

# Display the results
print(group_sym)
print(path.n)
print(sym.n)


### END - try to change NOX5 into ensembl ID ##




ncolors <- rainbow(length(sym.n),  s = 0.7, v = 1, alpha = 0.8)


names(ncolors) <- sym.n

aka4 = data.frame(Pathway=path.n)
aka4 <- as.matrix(aka4)
rownames(aka4) <- group_sym


column_split <- path.n

dt_file <- dt_file[, group_sym] %>% as.matrix()

dt_file <- dt_file[rownames(dt_file), ]


# replace the name on column_split with the values in remane.v
column_split_t <- plyr::revalue(column_split, rename.v)


str(dt_file)


ComplexHeatmap::pheatmap(dt_file,
                         scale="column",
                         treeheight_row = 0,
                         treeheight_col = 20, cellwidth = 10,
                         cellheight = 18, 
                         cluster_rows = F,
                         cluster_cols = F,
                         show_colnames = T,
                         annotation_legend = T,
                         annotation_names_row = F,
                         annotation_names_col = F,
                         border_color = "#BDBDBD",
                         # legend_labels = c(-1,0,1),
                         main = "",
                         # legend_breaks = c(-1,0,1),
                         name ="z-score",
                         legend = T,
                         fontsize = 9, 
                         color = colorRampPalette(c("navy","white","firebrick3"))(50),
                         column_split = column_split_t,
                         column_gap = unit(c(1), "mm",)
                         
)



## START - to check normalized value for NOX5 (duplicated genes) ##

# norm.count is a dataframe
# Extract row names
row_names <- rownames(norm.count)

# Create a new dataframe with rows that have NOX5 after the "=" sign
nox5_df <- norm.count[grepl("=NOX5$", row_names), ]

# Display the new dataframe
print(nox5_df)

# END - to check normalized value for NOX5 (duplicated genes) ##

# Final version for GitHub release (2025.04.03) — NOX5 dual Ensembl ID issue resolved