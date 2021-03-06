# Load GO and DE data to R

## @knitr prep_data

# Install missing CRAN packages if needed

# git.packages <- c("docopt/docopt.R")
# github.packages <- c("argparser")
install.deps <- function(p, repo="cran"){
  call_install <- switch(repo,cran=c("install.packages(package, repos=\"http://cloud.r-project.org/\""),
                         bioc=c("biocLite(package, suppressUpdates=TRUE"),
                         git=c("install.github(package"))
  if (repo=="bioc") eval(parse(text = getURL("http://bioconductor.org/biocLite.R", ssl.verifypeer=FALSE)))
  for (package in p) {
    if (!package %in% installed.packages()) {
      cat(sprintf("Please wait, installing and loading missing package %s and its dependencies from %s.\n",
                  package, repo), file=stderr())
      suppressWarnings(eval(parse(text = sprintf("%s, quiet = TRUE)",call_install))))
      if (!package %in% installed.packages()) {
        ifelse(!interactive(),q("no", 2, TRUE),
               stop(sprintf("Unable to install missing package %s from %s.\n",
                            package, repo), call. = FALSE))
      }
    }
    require(package, character.only=TRUE, quietly=TRUE, warn.conflicts=TRUE)
  }
}


# Install and load cran packages
CRAN_packages <- c("RCurl", "grid", "RColorBrewer", "dplyr", "tidyr"  ,"ggplot2", "gridExtra")
install.deps(CRAN_packages)

# Install and load bioconductor packages
bioc_packages <- c("qvalue","goseq", "GO.db", "gage", "pathview")
install.deps(bioc_packages, repo="bioc")

# Reload dplyr (must be last package loaded)
detach("package:dplyr", unload=TRUE)
library(dplyr)

# GOseq pathway enrichment analysis (based on trinity Goseq script and customized to fetch DE and annotation data from Trinotate sqlite db)

###### Functions ############

# Get and prepare relevant geneset and DE tables function
prepare_geneset_data <- function(Trinotatedb,Id_type, geneset="GO", minPfamScore=20, minBlastScore=100,  ko_table=NULL) {
  # Change SQL query based on DE level (transcript or orf) from DE analysis or specify manualy
  sql_select <- gsub(paste0("as ", Id_type, "_id,"), "as Trinity_Id,",  sub(paste0("as ", Id_type, "_length,"), "as length,", "SELECT O.orf_id as orf_id, O.transcript_id as transcript_id, O.length as orf_length, length(T.sequence) as transcript_length,"), perl=TRUE)

  Id_type_select <- switch(Id_type, "orf"=">", "transcript"="=")

  # GO annotation from PFAM (with Score>20)
  GO_sql <- sprintf(" substr(H.pfam_id,1,7) pfam_acc, G.GO_terms AS term, H.FullDomainScore FROM (SELECT * FROM HMMERDbase WHERE FullSeqScore>%s AND FullDomainScore>%s GROUP BY QueryProtID HAVING MAX(FullSeqScore) AND MAX(FullDomainScore)) H JOIN (SELECT pfam_acc, group_concat(go_id) AS GO_terms FROM pfam2go GROUP BY pfam_acc) G ON substr(H.pfam_id,1,7)=G.pfam_acc JOIN ORF O ON (H.QueryProtID=O.orf_id) JOIN Transcript T on O.transcript_id=T.transcript_id GROUP BY Trinity_Id HAVING MAX(H.FullSeqScore)",minPfamScore, minPfamScore )

  # COG Annotation from UniProt blast (with BitScore>100)
  COG_sql <- sprintf(' S.BitScore, S.Evalue,  E.eggNOGIndexTerm AS term, E.eggNOGDescriptionValue AS desc from BlastDbase S JOIN UniprotIndex U ON S.UniprotSearchString=U.Accession JOIN eggNOGIndex E ON U.LinkId=E.eggNOGIndexTerm JOIN ORF O on S.TrinityID=O.orf_id JOIN Transcript T on O.transcript_id=T.transcript_id WHERE U.AttributeType="E" AND instr(S.TrinityID, "m.")>0 AND S.BitScore>%s GROUP BY Trinity_Id HAVING MAX(S.BitScore)', minBlastScore)

  # KEGG Annotation from UniProt blastp (grouped by best BitScore ORF)
  KEGG_sql <- sprintf(' S.BitScore, S.Evalue,  substr(U.LinkId, 1, U.pos-1) AS ontology, substr(U.LinkId, U.pos+1, U.pos2-1) AS Kegg_species, substr(U.LinkId,U.pos2+U.pos+1) AS term from BlastDbase S JOIN (select *, instr(LinkId,":") AS pos, instr(substr(LinkId, instr(LinkId,":")+1),":") AS pos2 from UniprotIndex) U ON S.UniprotSearchString=U.Accession JOIN ORF O on S.TrinityID=O.orf_id JOIN Transcript T on O.transcript_id=T.transcript_id WHERE U.AttributeType="K" AND instr(S.TrinityID, "m.")>0 AND S.BitScore>%s GROUP BY Trinity_Id HAVING MAX(S.BitScore)', minBlastScore)

  # KEGG KO Annotation from KOBAS
  KO_sql <- sprintf(' K.KO_ID AS term, K.KO_name AS desc from (SELECT * FROM %s WHERE instr(TrinityID, "m.")%s0) K JOIN ORF O on K.TrinityID=O.%s JOIN Transcript T on O.transcript_id=T.transcript_id GROUP BY K.TrinityID HAVING min(O.rowid)',ko_table, Id_type_select, paste0(Id_type, "_id"))

  # construct the whole SQL query
  sql_query <- switch (tolower(geneset),
                       "go" = paste0(sql_select, GO_sql),
                       "cog" = paste0(sql_select, COG_sql),
                       "kegg" = paste0(sql_select, KEGG_sql),
                       "ko" = paste0(sql_select, KO_sql)
  )
 # geneset_tables <- list(geneset_data=NULL, geneset_DE_data=NULL)
  geneset_table <- tbl(Trinotatedb, sql(sql_query)) %>% collect() %>% filter(!is.na(.$term), term!="", term!="NA", term!="None")
  if (nrow(geneset_table)==0) stop(sprintf("No records were found in '%s' database.\nUsing the following query:\n %s\n", Trinotatedb, sql_query), call. = FALSE)
  if (tolower(geneset)=="ko") {
    # Setup KEGG_KO description table and terms
    cat(sprintf("Please wait, preparing KEGG orthologies information...\n"), file=stderr())
    kegg_ko <- kegg.gsets("ko")
    transcript_kegg_info_list <- sapply(geneset_table$term,
                          function(y) names(kegg_ko[[1]])[sapply(kegg_ko[[1]], FUN=function(x) y %in% x)])
    geneset_table <- geneset_table %>% mutate(term=sapply(.$term,
        function(l) paste(gsub('(ko[0-9]*).*','\\1',transcript_kegg_info_list[[l]]), collapse=","),
        USE.NAMES = FALSE)) %>% filter(term!="")
  }
  return(geneset_table)
}

# Calculate GSEA values for DE genes for a specific contrast and geneset analysis
GSEA <- function(de_table, geneset_data, contras, max_FDR=0.05, min_log2FC=2, description, annotation_source, analysis_date=format(Sys.Date(), "%d/%m/%Y"), geneset="GO") {
  contrast_levels <- unlist(strsplit(contras, "_vs_"))
  de_data <- de_table %>% filter(padj<=max_FDR, abs(log2FoldChange)>=min_log2FC,contrast==contras) %>% mutate(term=geneset_data$term[match(.$Trinity_Id, geneset_data$Trinity_Id)]) %>% filter(!is.na(.$term), term!="", term!="NA", term!="None")

  GOseq_result_table <- NULL
  for (i in 1:length(contrast_levels)) {
    # create a binary named list which marks with an 1 an ORF that is DE, and 0 if it's not (from all ORFs with GO)
    multip <- ifelse(i==1, 1, -1)
    de_ids <- de_data %>% filter(multip*log2FoldChange>=min_log2FC) %>%
      dplyr::select(Trinity_Id)
    cat_geneset_data_vec <- as.integer(geneset_data$Trinity_Id %in% de_ids$Trinity_Id)

    geneset_data_pwf <- nullp(cat_geneset_data_vec,bias.data=geneset_data$length,plot.fit = FALSE)
    rownames(geneset_data_pwf) = geneset_data$Trinity_Id
    #ORF_GO_info <- data.frame(orf_go$orf_id, orf_go$term)
    geneset_data_info <- data.frame(geneset_data$term, row.names = geneset_data$Trinity_Id)
    geneset_data_info_listed = apply(geneset_data_info, 1, function(x) unique(unlist(strsplit(x,','))))
    if (tolower(geneset)=="cog") geneset_data_info_listed <- data.frame(geneset_data$Trinity_Id, geneset_data$term)
    GOseq_res = goseq(geneset_data_pwf,gene2cat=geneset_data_info_listed)

    ## over-represented categories:
    pvals = GOseq_res$over_represented_pvalue
    pvals[pvals > 1 -1e-10] = 1-1e-10
    q = qvalue(pvals)
    GOseq_result_table <- GOseq_res %>% mutate(over_represented_FDR = q$qvalues, contrast=contras,
                                               Over_represented_in=contrast_levels[i], DE_padj=max_FDR,
                                               DE_log2FC=min_log2FC,
                                               analysis_date=analysis_date, annotation_source=annotation_source,
                                               DE_analysis=DE_analysis, description=description) %>%
      arrange(over_represented_FDR) %>% bind_rows(GOseq_result_table, .)
  }
    if (tolower(geneset)!="go") {
      # Prepare description dictionaries for COG and KO (Think about a solution for KEGG)
      #       term_dict <- unique(geneset_data[c("term","desc")])
      term_dict <- switch(tolower(geneset),
                          cog = unique(geneset_data[c("term","desc")]),
                          ko = as.data.frame(khier) %>% rename(ontology=category) %>% separate(pathway, c("term", "desc"), extra="merge") %>% mutate(term=paste0("ko", term))
      )
      # Set the appropriate ontology
      GOseq_result_table$ontology <- switch (tolower(geneset),
                  cog = ifelse(grepl("^ENOG.+", GOseq_result_table$category, perl=TRUE), "eggNOG","COG"),
                  kegg = sapply(de_data$Kegg_species[match(GOseq_result_table$category, de_data$term)],
                                  function(s) sub("(^.)", "\\U\\1", s, perl=TRUE),USE.NAMES = FALSE),
                  ko = term_dict$ontology[match(GOseq_result_table$category, term_dict$term)]
      )



      GOseq_result_table <- GOseq_result_table %>% mutate(term=term_dict$desc[match(.$category, term_dict$term)]) %>% filter(!is.na(.$term), term!="", term!="NA", term!="None")
       }

  return(GOseq_result_table)
}

# Load geneset analysis data back to the Trinotate sqlite db
GSEA2Trinotate <- function(trinotateDb, GSEA_results, tableName) {
  copy_to(trinotateDb, GSEA_results, name = tableName, temporary = FALSE,
          indexes = list( c("category", "ontology", "over_represented_FDR", "contrast", "DE_analysis")),
          types = c("TEXT", "REAL", "REAL","INTEGER","INTEGER","TEXT", "TEXT", "REAL","TEXT", "TEXT", "REAL", "REAL", "TEXT" , "TEXT", "TEXT", "TEXT"))
}

# Make term names shorter and prettier - better to do that vectorized and then figure out duplicates as well
pretty_term <- function(term, comma=TRUE, dot=TRUE, paren=TRUE, cap=TRUE, charNum=45) {
  term <- sapply(term, function (s) {
    s_ <- s
    if (cap) s_ <- sub("(^.)", "\\U\\1", s_, perl=TRUE)
    if (paren) s_ <- sub(" \\(.+\\)", "",  s_)
    if (dot) s_ <- unlist(strsplit(s_, ".", fixed = TRUE))[1]
    if (comma) s_ <- unlist(strsplit(s_, ",", fixed = TRUE))[1]
    # Trim to the end of the word if character number limit is reached and remove any trailing commas or dots
    if (charNum>0) s_ <- sub("[,|\\.]$", "", sub(sprintf("(^.{%s}[\\S]*)[ ]*.*", charNum), "\\1", s_, perl=TRUE), perl = TRUE)
    # Add a ^ at the end if the term was trimmed anyhow
    if (nchar(s_)<nchar(s) && !grepl("\\*$", s_)) return(paste0(s_, "^"))
    else return(s_)
  }, USE.NAMES = FALSE)
  return(term)
}

# vertical plot function
GSEA_vert_plot <- function(GSEA_set,geneset_analysis, cols, bar_cols, bar_width=0.4, facet=FALSE) {

  vert_plot <- ggplot(GSEA_set, aes(y=numDEInCat, x=reorder(cont_term, order),
                                    fill=Over_represented_in)) + labs(y="Number of DE genes", x="Term") +
    geom_bar(stat="identity",colour="black", width=bar_width)  +
    scale_y_continuous(expand = c(0,ceiling(max(GSEA_set$numDEInCat)*1.2)/70),                                                                                    limits=c(0,ceiling(max(GSEA_set$numDEInCat)*1.2))) +
    scale_x_discrete(labels = reorder(sub("([^_]+)_.+", "\\1", GSEA_set$cont_term), GSEA_set$order))

  GSEA_vert_theme <- theme_bw(base_size=28) +
    theme(axis.title.y=element_text( face="bold", vjust=2.5, size=rel(0.8)),
          axis.title.x=element_text(angle=90,face="bold", size=rel(0.8), vjust=1, hjust=10),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          legend.title=element_blank(),
          legend.text.align=1,
          #legend.position=c(length(GSEA_set$numDEInCat[GSEA_set$Over_represented_in==GSEA_set$Over_represented_in[1]])/ nrow(GSEA_set)-0.025,0.95), legend.justification=c(1,1),
          legend.position=c(which.min(GSEA_set$numDEInCat[GSEA_set$Over_represented_in==GSEA_set$Over_represented_in[1]])/ nrow(GSEA_set)+0.01,0.95), legend.justification=c(1,1),
          axis.text.y=element_text(angle=90, hjust=0.5), # 0.5
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.2, colour = as.vector(cols[GSEA_set$ontology])),
          panel.grid.major.x = element_blank())

  GSEA_vert_facet_theme <- theme_bw(base_size=28) +
    theme(axis.title.y=element_text( face="bold", vjust=2.5, size=rel(0.8)),
          axis.title.x=element_text(angle=90,face="bold", size=rel(0.8), vjust=1, hjust=10),
          strip.background=element_rect(fill="lightblue"),
          strip.text=element_text(size=rel(0.75)),
          legend.position="none",
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          axis.text.y=element_text(angle=90, hjust=0.5), # 0.5
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.2),
          panel.grid.major.x = element_blank())

  p <- vert_plot  +
       scale_fill_manual(values = bar_cols,
                      guide = guide_legend(direction = "horizontal",
                      label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,                                                   label.theme = element_text(angle = 90, size=20))) + GSEA_vert_theme
  #title.theme = element_text(angle = 90, face="bold", size=20)

  facet_plot <- vert_plot +
    scale_fill_manual(values = bar_cols[GSEA_set$Over_represented_in]) +
    facet_grid(. ~ Over_represented_in, scales="free_x", space="free_x") + GSEA_vert_facet_theme

  if (facet) facet_plot
  else p
}
# p

## @knitr horizontalPlot
# Horizontal plot function
GSEA_horiz_plot <- function(GSEA_set, cols,  bar_cols, bar_width=0.4, facet=FALSE) {

  # Define the theme
  GSEA_horizontal_theme <- theme_bw(base_size=18) +
    theme(axis.title.y=element_text(face="bold", vjust = 1.5, size=rel(0.8)),
          axis.title.x=element_text(face="bold", vjust = 0.1, size=rel(0.8)),
          legend.title=element_text(size=rel(0.8), hjust=0.5),
          legend.text=element_text(size = rel(0.7),lineheight = 1.5),
          panel.grid.major.x = element_blank(),
          legend.position=c(0.975,0.95), legend.justification=c(1,1),
          axis.text.x=element_text(angle=45, hjust=1, colour = as.vector(cols[GSEA_set$ontology])))


  horiz_plot <- ggplot(GSEA_set,
                       aes(y=numDEInCat, x=reorder(cont_cat, order), fill=Over_represented_in))
  hz <- horiz_plot +
    geom_bar(colour="black", width=bar_width, stat="identity", position="identity") +
    coord_fixed() + scale_fill_manual(values = bar_cols) + #[levels(GSEA_set$Over_represented_in)]
    geom_hline(yintercept=0) +
    scale_y_continuous(expand = c(0,0), limits=c(0,ceiling(max(GSEA_set$numDEInCat)*1.2))) +
    scale_x_discrete(labels = reorder(sub("([^_]+)_.+", "\\1", GSEA_set$cont_cat), GSEA_set$order)) +
    #  scale_fill_brewer(palette = "Set2", limits=c("Testis", "Ovary")) +
    labs(x="Term", y="Number of DE genes", fill="Over-represented\nin group\n") + GSEA_horizontal_theme

  # Define horizontal facet theme
  GSEA_horizontal_facet_theme <- theme_bw(base_size=18) +
    theme(axis.title.y=element_text(face="bold", vjust = 1.5, size=rel(0.8)),
          axis.title.x=element_text(face="bold", vjust = 0.1, size=rel(0.8)),
          strip.background=element_rect(fill="lightblue"),
          strip.text=element_text(size=rel(0.75)),
          legend.position="none",
          legend.title=element_text(size=rel(0.8), hjust=0.5),
          legend.text=element_text(size = rel(0.7),lineheight = 1.5),
          panel.grid.major.x = element_blank(),
          legend.position=c(0.975,0.95), legend.justification=c(1,1),
          axis.text.x=element_text(angle=45, hjust=1))

  hf <- horiz_plot +
    geom_bar(colour="black", width=bar_width, stat="identity", position="identity") +
    facet_grid(. ~ Over_represented_in, scales="free_x", space="free_x") +
    coord_fixed() + scale_fill_manual(values = bar_cols[GSEA_set$Over_represented_in]) +
    geom_hline(yintercept=0) +
    scale_y_continuous(expand = c(0,0), limits=c(0,ceiling(max(GSEA_set$numDEInCat)*1.2))) +
    scale_x_discrete(labels = reorder(sub("([^_]+)_.+", "\\1", GSEA_set$cont_cat), GSEA_set$order)) +
    #  scale_fill_brewer(palette = "Set2", limits=c("Testis", "Ovary")) +
    labs(x="Term", y="Number of DE genes", fill="Over-represented\nin group\n") + GSEA_horizontal_facet_theme
  if (facet) hf
  else hz
}

# Plot ontology legend (as separate plot)
plot_Ont_legend <- function(geneset_analysis=geneset_analysis, savePlot=TRUE, plotFormat="pdf", width=15, height=10, ont_cols=list(go_cols=c("Biological Process"="darkred", "Cell Cycle"="darkgreen", "Molecular Function"="darkblue"), cog_cols=c(COG="black", eggNOG="purple4"))){
  library(gridExtra)
  # sort out ontology colors for each level
  if (length(ont_cols)>0) {
    col_list <- paste0(tolower(geneset_analysis), "_cols")
    ontology_cols <- ont_cols[[col_list]]
  } else stop("No colours were defined for Ontology.", call. = FALSE)
  ont_legend <- data.frame(labels=names(ontology_cols), cols=as.vector(ontology_cols), x=0.05, stringsAsFactors = FALSE)
  ont_legend$y <- 0.8 - as.numeric(row.names(ont_legend))/10
  if (savePlot) {
    outDir=paste0(geneset_analysis, "_output_plots")
    dir.create(outDir, showWarnings = FALSE)
    plotFileName <- paste(paste(geneset_analysis, "ontology_legend", format(Sys.Date(), "%d_%m_%y"), sep="_"), plotFormat, sep=".")
    plotString <- paste0(  plotFormat, "(", "'",file.path(outDir,plotFileName), "'" ,",width=",width, ",height=",height, ")")
    eval(parse(text=plotString))
  }
  grid.newpage()
  vp1 <- viewport(width=0.45, height = 0.45)
  pushViewport(vp1)
  grid.rect()
  grid.text(paste(geneset_analysis,"Ontologies"), x = 0.05, y=0.825, just="left", gp=gpar(fontsize=23, fontface="bold"))
  grid.text(ont_legend$labels,ont_legend$x, ont_legend$y, just="left", gp=gpar(fontsize=20, col=ont_legend$cols))
  if (savePlot) dev.off()

}

# Produce pathview plots
produce_pathview <- function(deTable=DE_table,contras, fdr=max_FDR, koTable,pathway, outDir=paste0(geneset_analysis, "_pathway_output")){
  # create the data set for the de file

  deTable <- deTable %>% filter(contrast==contras, padj<=fdr)
  de <-  setNames(deTable$log2FoldChange, deTable$Trinity_Id)
  #     filter(KO_id!="None")
  ko_matrix <- mol.sum(de, as.data.frame(koTable[koTable$term!="None",c("Trinity_Id", "term")]), sum.method = "mean")
  colnames(ko_matrix) <- contras

  # Create an output folder for each DE comparison
  pathDir <- file.path(getwd(),outDir)
  try(dir.create(pathDir), silent = TRUE)
  currentDir <- getwd()
  setwd(ifelse(file.exists(pathDir), pathDir, getwd()))
  # use this to determine pathview limit: gene = max(abs(geneList))
  #genelist <- ko_matrix # change to sex matrix to analyse pathways in the other DE set
  # cat(sprintf("Producing top %s significant KEGG pathways from file %s in folder %s\n ",
  #            opts$pathview, deFile, pathDir), file=stderr())
  #topPathways <- koResults[1:as.numeric(opts$pathview),1]
  # Run Pathview to visualize expression data on KEGG pathways
  pv.out <- pathview(gene.data = ko_matrix, pathway.id = pathway, species = tolower(geneset_analysis),
             out.suffix = "DE_KO_pathview", keys.align = "y",
             kegg.native = T, match.data = F, multi.state = F, same.layer = T,
             limit = list(gene = max(abs(ko_matrix[,1])), cpd = 1))
  pathwayFiles <- paste(pathway, c("png", "xml"), sep=".")
  try(file.remove(pathwayFiles), silent = TRUE)
  setwd(currentDir)
}

# Save plot function
save_GSEA_Plot <- function(GSEA_set, orientation, rotate=FALSE, plotFormat="pdf", geneset=geneset_analysis, sig_value, width=10, height=15){
  outDir=paste0(geneset, "_output_plots")
  dir.create(outDir, showWarnings = FALSE)
  sig <- sub("(\\w+)[<=]+(\\d*.+)", "\\1_\\2", sig_value)
  plotFileName <- paste(paste(geneset, unique(GSEA_set$contrast), sig, orientation, format(Sys.Date(), "%d_%m_%y"), sep="_"), plotFormat, sep=".")
  ggsave(file.path(outDir,plotFileName), width=width, height=height)
  # flip the image (needed for markdown output and requires powershell in windows and ImageMagick in Linux)
  imageFlip <- function(angle=90, flipDir=outDir, pngName=plotFileName){
    fsep <- ifelse(.Platform$OS.type=="windows", "\\", "/")
    infile <- file.path(outDir,pngName, fsep=fsep)
    outfile <- file.path(outDir,paste0("Rotated_", pngName), fsep=fsep)
    if (.Platform$OS.type=="windows") {
      switch (plotFormat,
        png = system2("powershell", args=c("-ExecutionPolicy ByPass -command lib\\FlipImage.ps1",
                                           infile, outfile), wait = TRUE, invisible = TRUE),
        pdf = system2("powershell", args=c("-ExecutionPolicy ByPass -command lib\\Rotate_PDF_win.ps1",
                                         infile, outfile), wait = TRUE, invisible = TRUE)
      )

      error_msg <- sprintf("Could not rotate image or save it into <%s>, please verify that Powershell is installed and lib folder is available in the current working directory", outfile)
      } else {
        system2("convert",args=c(infile, "-rotate",angle, outfile),
                   wait = TRUE, invisible = TRUE)
        error_msg <- sprintf("Could not rotate image or save it into <%s>, please verify that ImageMagick is available in the PATH\n or download it from http://www.imagemagick.org/script/binary-releases.php#unix", outfile)
      }
    if (!file.exists(outfile)) message(error_msg)
  }
  if (rotate) imageFlip()

}


# Produce plots function
plotGSEA <- function(geneset_results, GSEA_filter="FDR<=0.1", cont, ont_cols=list(go_cols=c(BP="darkred", CC="darkgreen", MF="darkblue"), cog_cols=c(COG="black", eggNOG="purple4")), bar_cols="Set1", bar_width=0.4, groupOntology=1, orientation="vert", prettyTermOpts="comma=TRUE, charNum=35", facet=FALSE, savePlot=TRUE, saveFormat="pdf", rotateSavedPlot=FALSE, plot_width=10, plot_height=15) {
  cont_levels <- unlist(strsplit(cont, split = "_vs_", fixed = TRUE))
  # Prepare and filter the geneset by FDR or pvalue and then by contrast. Add needed columns
  GSEA_comparison <- geneset_results %>% filter_(paste0("over_represented_", GSEA_filter)) %>% filter(contrast==cont, !is.na(term), term!="", term!="NA")
  GSEA_comparison <- GSEA_comparison %>% mutate(short_term=eval(parse(text=sprintf("pretty_term(GSEA_comparison$term, %s)", prettyTermOpts))), cont_ont=paste(Over_represented_in, ontology, sep="_"), cont_term=paste(short_term, Over_represented_in,  sep="_"), cont_cat=paste(category, Over_represented_in,  sep="_"), ontology=as.character(ontology))
  # fix duplicate short terms
  for (u in unique(GSEA_comparison$cont_term[duplicated(GSEA_comparison$cont_term)])) {
    count=0
    for (j in 1:length(GSEA_comparison$cont_term)) {
      if (GSEA_comparison$cont_term[j]==u) {
        count=count+1
        GSEA_comparison$cont_term[j] <- sub("(^[^_]+)", paste0("\\1", paste0(rep("*", count), collapse = "")),sub("\\^+$", "", u),  perl=TRUE)
      }
    }
  }

  # Sort table and add plotting order
 # GSEA_comparison <- GSEA_comparison %>% arrange(Over_represented_in, desc(numDEInCat)) %>% mutate(order=1:nrow(.))
  # mutate(cont_term=factor(cont_term)) %>%


  # Assign ontology and bar colors for each contrast level
  geneset_analysis <- switch(sub("(^.).+", "\\U\\1", GSEA_comparison$category[1], perl=TRUE),
                             G="GO", C="COG",E="COG",K="KO")

  defaultBrewerPal <- "Set1"
  #brewerPal <- switch(geneset_analysis,"GO" = "Set1", "COG"="Set2")
  if(length(bar_cols)==length(levels(GSEA_comparison$Over_represented_in)) && all(levels(GSEA_comparison$Over_represented_in) %in% names(bar_cols))) {
    named_bar_cols <- bar_cols
    } else if (is.character(bar_cols) && length(bar_cols)==1) {
    tryCatch(brewerSet<-brewer.pal(length(levels(GSEA_comparison$Over_represented_in)), bar_cols), error=function(e) {message(sprintf("Warning: %sUsing %s instead", e$message, defaultBrewerPal)); brewerSet<-brewer.pal(length(levels(GSEA_comparison$Over_represented_in)), defaultBrewerPal)} )
    #print(brewerSet)
    named_bar_cols <- setNames(brewerSet, levels(GSEA_comparison$Over_represented_in))
    } else {
      message(sprintf("Warning: Bar colours for each contrast level were not defined, using RColorBrewer \"%s\" instead", defaultBrewerPal))
      named_bar_cols <- setNames(brewer.pal(length(levels(GSEA_comparison$Over_represented_in)), defaultBrewerPal), levels(GSEA_comparison$Over_represented_in))
    }
  # sort out ontology colors for each level
  # ontology_cols <- NULL
  if (length(ont_cols)>0) {
    #for (i in unique(GSEA_comparison$Over_represented_in)) {
      col_list <- paste0(tolower(geneset_analysis), "_cols")
      ontology_cols <-ont_cols[[col_list]]
      #ontology_cols <- c(ontology_cols, setNames(as.vector(ont_cols[[col_list]]), outer(i, names(ont_cols[[col_list]]),paste, sep="_")))
    #}
  }
  # Remove unused and reorder levels to plot by correct contrast order
  GSEA_comparison <- GSEA_comparison %>% mutate(Over_represented_in=factor(Over_represented_in, levels=rev(cont_levels)))
  # Sort by functional groups
  sort_funGroup <- switch (as.character(groupOntology),
                           "0" = NULL,
                           "1" = ifelse(orientation=="vert", "cont_ont", "cont_cat"),
                           "-1" = sprintf("desc(%s)", ifelse(orientation=="vert", "cont_ont", "cont_cat"))
  )
  sort_order <- c("Over_represented_in", sort_funGroup, "desc(numDEInCat)")
  GSEA_comparison <- GSEA_comparison %>% arrange_(.dots = as.list(sort_order)) %>% mutate(order=1:nrow(.))

  print(eval(parse(text=paste("GSEA", orientation, "plot(GSEA_comparison, cols=ontology_cols, bar_cols=named_bar_cols, bar_width=bar_width, facet=facet)", sep="_"))))
  if (savePlot) save_GSEA_Plot(GSEA_comparison, orientation, rotateSavedPlot, saveFormat, width=plot_width, height=plot_height, geneset = geneset_analysis, sig_value=GSEA_filter)
}

# Copy lib folder to the current working directory
