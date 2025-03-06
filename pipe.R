pipeline = list()


# pipeline[["test"]] = list(
#   input = list(
#     script= "/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/HMSC/test/test.Rmd"
#   ),
#   output = list(
#     report = "./test/my_testt.html"
#   )
# )
# 
# 
# make_with_recipe(
#   recipe = my_render(notebook_path =  "/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/HMSC/test/test.Rmd"),
#   targets = unlist(get_output( "/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/HMSC/test/test.Rmd")),
#   dependencies = unlist(get_input( "/sci/labs/yotamd/lab_share/avishai.wizel/R_projects/HMSC/test/test.Rmd")),
#   label = "create_data",build = F
# )

####################################### HMSC-ACC preprocess ####################################################


pipeline[["create_data"]] = list(
  input = list(
    script= "./Notebooks/HMSC/01_preprocess/01_create_data.Rmd",
    counts= "Input_data/all.txt"
  ),
  output = list(
    report = "./Reports/HMSC/01_preprocess/01_create_data/01_create_data.html",
    all_cells_seurat = "./Reports/HMSC/01_preprocess/01_create_data/acc_tpm_nCount_mito_no146_15k_with_ACC1.RDS"
  )
)


pipeline[["Cancer_filtering"]] = list(
  input = list(
    script= "./Notebooks/HMSC/01_preprocess/02_Cancer_filtering.Rmd",
    all_cells_seurat= pipeline$create_data$output$all_cells_seurat,
    kaye_acc_genes = "./Input_data/preprocess/oncotarget-05-12528-s001_acchigh.txt"
  ),
  output = list(
    report = "./Reports/HMSC/01_preprocess/02_Cancer_filtering/02_Cancer_filtering.html",
    cell_types = "./Reports/HMSC/01_preprocess/02_Cancer_filtering/cell_types.RDS",
    acc_cancer_cells = "./Reports/HMSC/01_preprocess/02_Cancer_filtering/acc_cancer_cells.RDS"
  )
)

pipeline[["Cancer_filtering"]] = list(
  input = list(
    script= "./Notebooks/HMSC/01_preprocess/02_Cancer_filtering.Rmd",
    all_cells_seurat= pipeline$create_data$output$all_cells_seurat,
    kaye_acc_genes = "./Input_data/preprocess/oncotarget-05-12528-s001_acchigh.txt"
  ),
  output = list(
    report = "./Reports/HMSC/01_preprocess/02_Cancer_filtering/02_Cancer_filtering.html",
    cell_types = "./Reports/HMSC/01_preprocess/02_Cancer_filtering/cell_types.RDS",
    acc_cancer_cells = "./Reports/HMSC/01_preprocess/02_Cancer_filtering/acc_cancer_cells.RDS"
  )
)

pipeline[["HMSC_cells_preprocess"]] = list(
  input = list(
    script = "./Notebooks/HMSC/01_preprocess/03_HMSC_cells_preprocess.Rmd",
    acc_cancer_cells = pipeline$Cancer_filtering$output$acc_cancer_cells
  ),
  output = list(
    report = "./Reports/HMSC/01_preprocess/03_HMSC_cells_preprocess/03_HMSC_cells_preprocess.html",
    hmsc_integrated = "./Reports/HMSC/01_preprocess/03_HMSC_cells_preprocess/acc1_cancer_cells_2500features_integrated_V5.RDS"
  )
)

pipeline[["build_datasets"]] = list(
  input = list(
    script = "./Notebooks/HMSC/01_preprocess/04_build_datasets.Rmd",
    hmsc_integrated = pipeline$HMSC_cells_preprocess$output$hmsc_integrated,
    acc_cancer_cells = pipeline$Cancer_filtering$output$acc_cancer_cells
  ),
  output = list(
    report = "./Reports/HMSC/01_preprocess/04_build_datasets/04_build_datasets.html",
    acc_hmsc_integrated = "./Reports/HMSC/01_preprocess/04_build_datasets/hmsc_acc_pri_cancer_processed.RDS"
  )
)

####################################### HMSC analysis ####################################################

pipeline[["Lum-Myo_HMSC"]] = list(
  input = list(
    script = "./Notebooks/HMSC/03_Lum-Myo_HMSC.Rmd",
    hmsc_cancer_cells = pipeline$HMSC_cells_preprocess$output$hmsc_integrated,
    hmsc_acc_pri_cancer = pipeline$build_datasets$output$acc_hmsc_integrated,
    genesets_h = "./Input_data/h.all.v7.0.symbols.pluscc.gmt"
  ),
  output = list(
    report = "./Reports/HMSC/03_Lum-Myo_HMSC/03_Lum-Myo_HMSC.html"
  )
)

pipeline[["HMSC_vs_ACC"]] = list(
  input = list(
    script = "./Notebooks/HMSC/04_HMSC_vs_ACC.Rmd",
    hmsc_cancer_cells = pipeline$HMSC_cells_preprocess$output$hmsc_integrated,
    hmsc_acc_pri_cancer = pipeline$build_datasets$output$acc_hmsc_integrated,
    genesets_h = "./Input_data/h.all.v7.0.symbols.pluscc.gmt",
    ensembl =  "./Input_data/hsapiens_gene_ensembl_version_113",
    fc_params  = "./Input_data/FC_params.R"
  ),
  output = list(report = "./Reports/HMSC/04_HMSC_vs_ACC/04_HMSC_vs_ACC.html")
)



####################################### HMSC cNMF ####################################################
pipeline[["run_cnmf"]] = list(
  input = list(
    script = "./Notebooks/HMSC/06_cNMF/cNMF_HMSC_01_run.Rmd",
    hmsc_cancer_cells = pipeline$HMSC_cells_preprocess$output$hmsc_integrated
  ),
  output = list(report = "./Reports/HMSC/06_cNMF/cNMF_HMSC_01_run.html",
                cnmf_object = "./Reports/HMSC/06_cNMF/HMSC_BatchCorrected_NMF/cnmf_obj.pkl")
)

pipeline[["cnmf_get_results"]] = list(
  input = list(
    script = "./Notebooks/HMSC/06_cNMF/cNMF_HMSC_02_get_results.Rmd",
    cnmf_object = pipeline$run_cnmf$output$cnmf_object
  ),
  output = list(report = "./Reports/HMSC/06_cNMF/cNMF_HMSC_02_get_results.html",
                cnmf_object = "./Reports/HMSC/06_cNMF/HMSC_BatchCorrected_NMF/cnmf_obj.pkl")
)

#############################################OPSCC###############################
pipeline[["OPSCC_preprocess"]] = list(
  input = list(
    script = "./Notebooks/HPV_OPSCC_Analysis_PMC10191634/01_preprocess.Rmd",
    matrix = "./Input_data/HPV_OPSCC_Analysis_PMC10191634/GSE182227_OPSCC.mtx.gz",
    barcodes = "./Input_data/HPV_OPSCC_Analysis_PMC10191634/GSE182227_barcodes.txt.gz",
    genes = "./Input_data/HPV_OPSCC_Analysis_PMC10191634/GSE182227_genes.txt.gz"
  ),
  output = list(report = "./Reports/HPV_OPSCC_Analysis_PMC10191634/01_preprocess/01_preprocess.html",
                opscc = "./Reports/HPV_OPSCC_Analysis_PMC10191634/01_preprocess/opscc_suerat.RDS")
)

pipeline[["OPSCC_deg"]] = list(
  input = list(
    script = "./Notebooks/HPV_OPSCC_Analysis_PMC10191634/02_findDEG.Rmd",
    opscc = pipeline$OPSCC_preprocess$output$opscc,
    opscc_metadata = "./Input_data/HPV_OPSCC_Analysis_PMC10191634/NIHMS1892753-supplement-Supp_Tables.xlsx"
  ),
  output = list(
    report = "./Reports/HPV_OPSCC_Analysis_PMC10191634/02_findDEG/02_findDEG.html",
    opscc_deg = "./Reports/HPV_OPSCC_Analysis_PMC10191634/02_findDEG/OPSCC_HPV_DEG.csv")
)



################## HMSC HPV analysis ########################################

pipeline[["HPV_analysis"]] = list(
  input = list(
    script = "./Notebooks/HMSC/05_HPV_analysis.Rmd",
    hmsc_cancer_cells = pipeline$HMSC_cells_preprocess$output$hmsc_integrated,
    hpv_conuts_data_p3 = "./Input_data/HMSC_HPV_data/HPV33_P3.txt",
    hpv_conuts_data_p2 = "./Input_data/HMSC_HPV_data/HPV33_P2.txt",
    opscc_deg = pipeline$OPSCC_deg$output$opscc_deg,
    genesets_h = "./Input_data/h.all.v7.0.symbols.pluscc.gmt",
    genesets_cp_pid = "./Input_data/c2.cp.pid.v2024.1.Hs.symbols.gmt",
    ensembl =  "./Input_data/hsapiens_gene_ensembl_version_113",
    myb_targets = "./Input_data/MYB_tergets_PMC4767593_sup_table_4.xlsx",
    fc_params  = "./Input_data/FC_params.R"
  ),
  output = list(
    report = "./Reports/HMSC/05_HPV_analysis/05_HPV_analysis.html",
    hmsc_deg = "./Reports/HMSC/05_HPV_analysis/HMSC_hpv_deg_df.csv"
  )
)


##################################### OPSCC analysis continue ###################################
pipeline[["OPSCC_volcano"]] = list(
  input = list(
    script = "./Notebooks/HPV_OPSCC_Analysis_PMC10191634/03_volcano_plot.Rmd",
    hmsc_deg = pipeline$HPV_analysis$output$hmsc_deg,
    opscc_deg = pipeline$OPSCC_deg$output$opscc_deg,
    fc_params  = "./Input_data/FC_params.R"
  ),
  output = list(
    report = "./Reports/HPV_OPSCC_Analysis_PMC10191634/03_volcano_plot/03_volcano_plot.html"  )
)

pipeline[["OPSCC_MYB_fisher"]] = list(
  input = list(
    script = "./Notebooks/HPV_OPSCC_Analysis_PMC10191634/05_MYB_HPV_fisher.Rmd",
    opscc = pipeline$OPSCC_preprocess$output$opscc,
    myb_targets = "./Input_data/MYB_tergets_PMC4767593_sup_table_4.xlsx"
  ),
  output = list(
    report = "./Reports/HPV_OPSCC_Analysis_PMC10191634/05_MYB_HPV_fisher/05_MYB_HPV_fisher.html"  )
)

pipeline[["OPSCC_HMSC_HPV_signature"]] = list(
  input = list(
    script = "./Notebooks/HPV_OPSCC_Analysis_PMC10191634/06_HMSC_HPV_signature.Rmd",
    opscc = pipeline$OPSCC_preprocess$output$opscc,
    hmsc_deg = pipeline$HPV_analysis$output$hmsc_deg,
    fc_params  = "./Input_data/FC_params.R"
  ),
  output = list(
    report = "./Reports/HPV_OPSCC_Analysis_PMC10191634/06_HMSC_HPV_signature/06_HMSC_HPV_signature.html"  )
)

############################## HPV signature compare ####################

pipeline[["HPV_signature_compare"]] = list(
  input = list(
    script = "./Notebooks/integrative_analysis/HPV_signature_compare.Rmd",
    opscc_deg = pipeline$OPSCC_deg$output$opscc_deg,
    hmsc_deg = pipeline$HPV_analysis$output$hmsc_deg,
    genesets_h = "./Input_data/h.all.v7.0.symbols.pluscc.gmt",
    fc_params  = "./Input_data/FC_params.R"
  ),
  output = list(
    report = "./Reports/integrative_analysis/HPV_signature_compare.html"  )
)

####################################### TCGA #####################################################
pipeline[["TCGA_HNSC_analysis"]] = list(
  input = list(
    script =  "./Notebooks/TCGA/OPSCC_TCGA.Rmd",
    opscc = pipeline$OPSCC_preprocess$output$opscc,
    hmsc_deg = pipeline$HPV_analysis$output$hmsc_deg,
    functions = "Notebooks/TCGA/TCGA_functions.R",
    opscc_deg =  pipeline$OPSCC_deg$output$opscc_deg,
    myb_targets = "./Input_data/MYB_tergets_PMC4767593_sup_table_4.xlsx",
    tp53_data = "./Input_data/TCGA/cbioportal_hnsc_tcga_gdc_TP53_data.tsv",
    genesets_h = "./Input_data/h.all.v7.0.symbols.pluscc.gmt",
    HNSC_tpm = "./Reports/TCGA/build_TCGA_datasets/TCGA_HNSC_TPM.RDS",
    cbp_data = "./Input_data/TCGA/cbioportal_hnsc_tcga_pan_can_atlas_2018_clinical_data.tsv",
    clinical_follow_up = "./Input_data/TCGA/clinical.project-tcga-hnsc.2025-03-04/follow_up.tsv",
    fc_params  = "./Input_data/FC_params.R"
    
  ),
  output = list(
    report = "./Reports/TCGA/OPSCC_TCGA/OPSCC_TCGA.html"  )
)


######################################## functions ###############################################
#get input/output from pipeline with script name
get_input <- function(script) {
  for (i in seq_along(pipeline)) {
    item <- pipeline[[i]]
    if (is.list(item) && "input" %in% names(item) && "script" %in% names(item$input)) {
      if (item$input$script == script) {
        return(pipeline[[i]]$input)
      }
    }
  }
}

get_output <- function(script) {
  for (i in seq_along(pipeline)) {
    item <- pipeline[[i]]
    if (is.list(item) && "input" %in% names(item) && "script" %in% names(item$input)) {
      if (item$input$script == script) {
        return(pipeline[[i]]$output)
      }
    }
  }
}

######################################## set makepipe ######################################

library(makepipe)
makepipe::reset_pipeline()
pipe <- get_pipeline()


make_with_recipe(
  recipe = my_render(notebook_path = "./Notebooks/HMSC/01_preprocess/01_create_data.Rmd"),
  targets = unlist(get_output("./Notebooks/HMSC/01_preprocess/01_create_data.Rmd")),
  dependencies = unlist(get_input("./Notebooks/HMSC/01_preprocess/01_create_data.Rmd")),
  label = "create_data",build = F
)


# 02_Cancer_filtering.Rmd
script = "./Notebooks/HMSC/01_preprocess/02_Cancer_filtering.Rmd"
make_with_recipe(
  recipe = my_render(notebook_path ="./Notebooks/HMSC/01_preprocess/02_Cancer_filtering.Rmd"),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = "Cancer_filtering",build = F
)


#03_HMSC_cells_preprocess.Rmd
script = "./Notebooks/HMSC/01_preprocess/03_HMSC_cells_preprocess.Rmd"
make_with_recipe(
  recipe = my_render(notebook_path = "./Notebooks/HMSC/01_preprocess/03_HMSC_cells_preprocess.Rmd"),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = "HMSC_cells_preprocess",build = F
)

#04_build_datasets.Rmd
script = "./Notebooks/HMSC/01_preprocess/04_build_datasets.Rmd"
make_with_recipe(
  recipe = my_render("./Notebooks/HMSC/01_preprocess/04_build_datasets.Rmd"),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = "build_datasets",build = F
)


# #02_CNV.Rmd
# notebook_path = "./Notebooks/HMSC/02_CNV.Rmd"
# targets_and_depens = get_targets_and_depens(notebook_path)
# 
# 
# make_with_recipe(
#   recipe = my_render("./Notebooks/HMSC/02_CNV.Rmd"),
#     targets = c(targets_and_depens$targets),
#     dependencies = c(targets_and_depens$dependencies,notebook_path))


#03_Lum-Myo_HMSC.Rmd
script = "./Notebooks/HMSC/03_Lum-Myo_HMSC.Rmd"
make_with_recipe(
  recipe = my_render("./Notebooks/HMSC/03_Lum-Myo_HMSC.Rmd"),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = "Lum-Myo_HMSC",build = F
)


#04_HMSC_vs_ACC.Rmd
script = "./Notebooks/HMSC/04_HMSC_vs_ACC.Rmd"

make_with_recipe(
  recipe = my_render( "./Notebooks/HMSC/04_HMSC_vs_ACC.Rmd"),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = "HMSC_vs_ACC",build = F
)

#05_HPV_analysis.Rmd
script = "./Notebooks/HMSC/05_HPV_analysis.Rmd"
make_with_recipe(
  recipe = my_render( "./Notebooks/HMSC/05_HPV_analysis.Rmd"),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = "HPV_analysis",build = F
)

#OPSCC preprocess
script = "./Notebooks/HPV_OPSCC_Analysis_PMC10191634/01_preprocess.Rmd"
make_with_recipe(
  recipe = my_render("./Notebooks/HPV_OPSCC_Analysis_PMC10191634/01_preprocess.Rmd"),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = "OPSCC_preprocess",build = F
)


#OPSCC DEG
script = "./Notebooks/HPV_OPSCC_Analysis_PMC10191634/02_findDEG.Rmd"
make_with_recipe(
  recipe = my_render("./Notebooks/HPV_OPSCC_Analysis_PMC10191634/02_findDEG.Rmd"),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = "OPSCC_DEG",build = F
)

#OPSCC volcano

script = "./Notebooks/HPV_OPSCC_Analysis_PMC10191634/03_volcano_plot.Rmd"
make_with_recipe(
  recipe = my_render("./Notebooks/HPV_OPSCC_Analysis_PMC10191634/03_volcano_plot.Rmd"),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = "OPSCC_volcano",build = F
)


#OPSCC MYB fisher

script = "./Notebooks/HPV_OPSCC_Analysis_PMC10191634/05_MYB_HPV_fisher.Rmd"
make_with_recipe(
  recipe = my_render("./Notebooks/HPV_OPSCC_Analysis_PMC10191634/05_MYB_HPV_fisher.Rmd"),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = "OPSCC_MYB_fisher",build = F
)

#OPSCC HMSC_HPV_signature

script = "./Notebooks/HPV_OPSCC_Analysis_PMC10191634/06_HMSC_HPV_signature.Rmd"

make_with_recipe(
  recipe = my_render("./Notebooks/HPV_OPSCC_Analysis_PMC10191634/06_HMSC_HPV_signature.Rmd"),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = "OPSCC_MYB_fisher",build = F
)


# #TCGA built datasets
# notebook_path = "./Notebooks/TCGA/build_TCGA_datasets.Rmd"
# targets_and_depens = get_targets_and_depens(notebook_path)
# 
# make_with_recipe(
#   recipe = my_render("./Notebooks/TCGA/build_TCGA_datasets.Rmd"),
#     targets = c(targets_and_depens$targets),
#     dependencies = c(targets_and_depens$dependencies,notebook_path),
#   label = "TCGA built datasets")

#HPV signature compare
script = "./Notebooks/integrative_analysis/HPV_signature_compare.Rmd"
make_with_recipe(
  recipe = my_render("./Notebooks/integrative_analysis/HPV_signature_compare.Rmd"),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = "HPV_signature_compare",build = F
)


#TCGA HNSC
script = "./Notebooks/TCGA/OPSCC_TCGA.Rmd"
make_with_recipe(
  recipe = my_render("./Notebooks/TCGA/OPSCC_TCGA.Rmd"),
  targets = unlist(get_output(script)),
  dependencies = unlist(get_input(script)),
  label = "TCGA_HNSC_analysis",build = F
)


