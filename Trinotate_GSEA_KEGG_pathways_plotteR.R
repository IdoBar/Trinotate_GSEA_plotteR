# Functions from file
source("I:/Academic/Bioinformatics/tools/R_scripts/GSEA/Trinotate_GSEA_plotteR/Trinotate_GSEA_plotteR.R")

#### Perform Goseq analysis ##########

# define parameters
geneset_analysis <- "GO"
max_FDR <- 0.05
min_log2FC <- 2
min_eValue <- 1e-5
minPfamScore <- 20
minBlastScore <- 100


# Fetch ORF data from Trinotate.sqlite
Trinotate <- src_sqlite("../../../de_novo_Trinotate_db/Lentils_Trinotate.sqlite")

# Get specific DE analysis table (either from trinotate or upload from a tab-delimited edgeR or DESeq2 output file with a "contrast" column sepcifying the DE contrast(s))
DE_analysis <- "ORF_DE_kallisto"
DE_table <- tbl(Trinotate, DE_analysis) %>% collect()
# DE_table <- read.delim("DE_Analysis_filename")

TrinityId_type <- ifelse(grepl("m\\.\\d+", DE_table[1,1], perl = TRUE),"orf", "transcript")
colnames(DE_table)[1] <- "Trinity_Id"
if (grepl("edgeR", DE_analysis, ignore.case = TRUE)) {
  DE_table <- DE_table %>% dplyr::rename(log2FoldChange=logFC, padj=FDR)
}

contrasts <- grep("_vs_",levels(factor(DE_table$contrast)), value = TRUE)


############## GO analysis  #####################
# details of GO analysis
geneset_analysis <- "GO"
GSEA_annotation_source <- "PFAM"
GSEA_description <- sprintf("%s annotation from %s, with FullDomainScore>20>%s, based on %s" , geneset_analysis, GSEA_annotation_source, minPfamScore, DE_analysis)

# Fetch and process GO data
geneset_data <- prepare_geneset_data(Trinotate, TrinityId_type, geneset = geneset_analysis)
geneset_results <- bind_rows(lapply(contrasts, function(x) GSEA(DE_table, geneset_data, contras = x, description = GSEA_description, annotation_source = GSEA_annotation_source, geneset = geneset_analysis))) %>% mutate_each_(funs(factor), c("contrast", "Over_represented_in"))

# deposit the results back to a table in Trinotate
# Load data to Trinotate sqlite db
analysis_name <- paste(geneset_analysis, DE_analysis, sep="_")
# GSEA2Trinotate(Trinotate, geneset_results, analysis_name)

# Manually set bar colours (Use Set1, but remove the purple which is too similar to the red):
brewerSet<-brewer.pal(9, "Set1")[-c(4:5,8)]
# Create a named vector with color for each contrast level
bar_color_names <- setNames(brewerSet, levels(geneset_results$Over_represented_in))

# Plot just 1 contrast:
#plotGSEA(geneset_results, cont = contrasts[5], GSEA_filter = "FDR<=0.1", groupOntology = 1, bar_cols = bar_color_names, bar_width = 0.45, savePlot = TRUE, saveFormat = "pdf", rotateSavedPlot = TRUE, prettyTermOpts = "charNum=35, comma=TRUE", plot_width = 15, plot_height = 15, facet=FALSE)

# plot all contrasts (vertical, no facets).
sapply(levels(geneset_results$contrast), function(x) plotGSEA(geneset_results, cont = x, bar_cols = bar_color_names, bar_width = 0.45, GSEA_filter = "FDR<=0.1", savePlot = TRUE, prettyTermOpts = "charNum=35, comma=TRUE", plot_width = 15, plot_height = 15))

# Plot GO ontology legend (as separate plot)
plot_Ont_legend(geneset_analysis = "GO")


###########  COG analysis #################
# details of COG analysis
geneset_analysis <- "COG"
GSEA_annotation_source <- "UniProt_eggNOG"
GSEA_description <- sprintf("%s annotation from %s, with BitScore>%s, based on %s" , geneset_analysis, GSEA_annotation_source, minBlastScore, DE_analysis)

# Fetch and process COG data
geneset_data <- prepare_geneset_data(Trinotate, TrinityId_type, geneset = geneset_analysis)
geneset_results <- bind_rows(lapply(contrasts, function(x) GSEA(DE_table, geneset_data, contras = x, description = GSEA_description, annotation_source = GSEA_annotation_source, geneset = geneset_analysis))) %>% mutate_each_(funs(factor), c("contrast", "Over_represented_in"))

# deposit the results back to a table in Trinotate
# Load data to Trinotate sqlite db
analysis_name <- paste(geneset_analysis, DE_analysis, sep="_")
# GSEA2Trinotate(Trinotate, geneset_results, analysis_name)

# Manually set bar colours (Use Set1, but remove the purple which is too similar to the red):
#brewerSet<-brewer.pal(length(levels(geneset_results$Over_represented_in))+1, "Set1")[-4]
# Create a named vector with color for each contrast level
#bar_color_names <- setNames(brewerSet, levels(geneset_results$Over_represented_in))
# Plot just 1 contrast:
#plotGSEA(geneset_results, cont = contrasts[4], GSEA_filter = "FDR<=0.1", groupOntology = 0, ont_cols = NULL, bar_cols = bar_color_names, savePlot = FALSE, saveFormat = "pdf", rotateSavedPlot = FALSE, prettyTermOpts = "charNum=35, comma=TRUE", plot_width = 15, plot_height = 15)

# plot all contrasts (vertical, no facets).
sapply(levels(geneset_results$contrast), function(x) plotGSEA(geneset_results, cont = x, GSEA_filter = "FDR<=0.1", groupOntology = 0, ont_cols = NULL, bar_cols = bar_color_names, savePlot = TRUE, saveFormat = "pdf", rotateSavedPlot = FALSE, prettyTermOpts = "charNum=35, comma=TRUE", plot_width = 15, plot_height = 15))# GSEA_filter="pvalue<=0.01"

# Plot GO ontology legend (as separate plot)
# plot_Ont_legend(geneset_analysis = "COG")

#####################  KO analysis ######################
# details of KO analysis
ko_table <- "KOBAS_KO"
geneset_analysis <- "KO"
GSEA_annotation_source <- "KOBAS"
GSEA_description <- sprintf("%s annotation from %s, with eValue<%s, based on %s" , geneset_analysis, GSEA_annotation_source, min_eValue, DE_analysis)

# Fetch and process KO data
geneset_data <- prepare_geneset_data(Trinotate, TrinityId_type, geneset = geneset_analysis, ko_table = ko_table)
geneset_results <- bind_rows(lapply(contrasts, function(x) GSEA(DE_table, geneset_data, contras = x, description = GSEA_description, annotation_source = GSEA_annotation_source, geneset = geneset_analysis))) %>% mutate_each_(funs(factor), c("contrast", "Over_represented_in"))


#"ontology"
# deposit the results back to a table in Trinotate
# Load data to Trinotate sqlite db
analysis_name <- paste(geneset_analysis, DE_analysis, sep="_")

# Manually set bar colours (Use Set1, but remove the purple and orange which are too similar to the red):
brewerSet<-brewer.pal(9, "Set1")[-c(4:5,8)]
# Create a named vector with color for each contrast level
bar_color_names <- setNames(brewerSet, levels(geneset_results$Over_represented_in))

# set ontology colours
ko_ont_cols <- list(ko_cols=setNames(brewer.pal(length(unique(geneset_results$ontology)), "Dark2"), unique(geneset_results$ontology)))
# Plot just 1 contrast:
#plotGSEA(geneset_results, cont = contrasts[5], groupOntology = -1, ont_cols = ko_ont_cols, GSEA_filter = "FDR<=0.1", bar_cols = bar_color_names, bar_width = 0.45, prettyTermOpts = "comma=FALSE, charNum=35", savePlot = TRUE, saveFormat = "pdf", rotateSavedPlot = TRUE, plot_width = 15, plot_height = 15)

# plot all contrasts (vertical, no facets).
sapply(levels(geneset_results$contrast), function(x) plotGSEA(geneset_results, cont = x, groupOntology = 1, ont_cols = ko_ont_cols, GSEA_filter = "pvalue<=0.05", bar_cols = bar_color_names, bar_width = 0.45, prettyTermOpts = "comma=FALSE, charNum=35", savePlot = TRUE, saveFormat = "pdf", rotateSavedPlot = FALSE, plot_width = 15, plot_height = 15))# GSEA_filter="pvalue<=0.01"

# Plot ontology legend (as separate plot)
plot_Ont_legend(geneset_analysis = geneset_analysis, ont_cols = ko_ont_cols)

