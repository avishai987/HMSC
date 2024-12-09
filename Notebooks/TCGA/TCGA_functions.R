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
