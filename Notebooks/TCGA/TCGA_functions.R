prepare_datasets <- function(project_name, filter = list(sample_type = c("Primary Tumor"))) {
  library(data.table)
  
  # from https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/493#issuecomment-1083311661
  # Defines the query to the GDC
  query <- GDCquery(project = project_name,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    experimental.strategy = "RNA-Seq",
                    workflow.type = "STAR - Counts")
  
  # Get metadata matrix
  metadata <- query[[1]][[1]]
  
  # Download data using api
  GDCdownload(query, method = "api")
  
  # Get main directory where data is stored
  main_dir <- file.path("GDCdata", project_name)
  # Get file list of downloaded files
  file_list <- file.path("GDCdata", project_name,list.files(main_dir,recursive = TRUE)) 
  
  # Read first downloaded to get gene names
  test_tab <- read.table(file = file_list[1], sep = '\t', header = TRUE)
  # Delete header lines that don't contain usefull information
  test_tab <- test_tab[-c(1:4),]
  # STAR counts and tpm datasets
  tpm_data_frame <- data.frame(test_tab[,2])
  count_data_frame <- data.frame(test_tab[,2])
  
  # Append cycle to get the complete matrix
  for (i in c(1:length(file_list))) {
    # Read table
    test_tab <- fread(file = file_list[i],sep = "\t",header = T,data.table = F)
    # Delete not useful lines
    test_tab <- test_tab[-c(1:4),]
    # Column bind of tpm and counts data
    tpm_data_frame <- cbind(tpm_data_frame, test_tab[,9])
    sample_id = strsplit(file_list[i],split = "/")[[1]][6]
    sample_name = metadata[metadata$id == sample_id ,"cases.submitter_id"]
    colnames(tpm_data_frame)[i+1] = sample_name # col 1 is gene
    # count_data_frame <- cbind(count_data_frame, test_tab[,4])
    # Print progres from 0 to 1
    print(i/length(file_list))
  }
  # update col names
  colnames(tpm_data_frame)[1] = "Gene"
  genes = tpm_data_frame$Gene
  tpm_data_frame$Gene <- NULL
  
  
  
  # remove all 0 genes
  not_all_zero <-apply(tpm_data_frame, 1, function(row) any(row != 0))
  tpm_data_frame = tpm_data_frame[not_all_zero,]
  genes = genes[not_all_zero]
  
  
  # remove all rows that are not gene symbols
  ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  results <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol","transcript_biotype"), mart=ensembl)
  not_symbol = !genes %in% results$hgnc_symbol
  genes = genes[!not_symbol]
  tpm_data_frame = tpm_data_frame[!not_symbol,]
  rownames(tpm_data_frame) = make.unique(genes) 
  
  for (i in 1:length(filter)) {
    col = names(filter)[i]
    values = filter[[i]]
    primary_samples = metadata[metadata[[col]] %in% values,"cases.submitter_id"]
    tpm_data_frame = tpm_data_frame[,primary_samples]
  }
  return(tpm_data_frame)
}

set_clinical_data <- function(clin_data, genes, tpm_data_frame,stratify = ("Q"),signature_name = NULL) {
  
  if (stratify == "Q") {
    bins = seq(0,1,0.25) #split to quantiles
    low = 2; high=4 #  Q1 low expression, Q3 high
  } 
  else if (stratify == "M") {
    bins = seq(0,1,0.25) #split to quantiles
    low = 3
    high = 3 #  stratify by median
  } 
  else if(stratify == "T"){
    bins = seq(0,1,1/3) #split to terciles
    low = 2
    high = 3 #  stratify by median
  }
  
  
  if (length(genes) == 1) {
    name = genes
  }else if (! is_null(signature_name)) {
    name = signature_name
  }else{
    name = "score"
  }
  
  #filter only for samples in clin_data
  gene_data = tpm_data_frame[rownames(tpm_data_frame) %in% genes, colnames(tpm_data_frame) %in% clin_data$submitter_id]
  genes_average = as.data.frame(rowMeans(as.data.frame(t(gene_data))))

  colnames(genes_average) = "mean_score"
  
  
  
  
  gene_status = genes_average %>% mutate (
    gene_status = case_when(
      mean_score > quantile(mean_score,bins)[high] ~ paste("High", name),
      mean_score < quantile(mean_score,bins)[low] ~  paste("Low", name),
      TRUE ~ "same"
    )
  )
  rownames(gene_status) = gsub(x = rownames(gene_status),
                               pattern = "\\.",
                               replacement = "-")
  gene_status = gene_status[gene_status$gene_status != "same",]
  
  colnames(gene_status)[colnames(gene_status) == "mean_score"] <- name
  return(gene_status)
}
