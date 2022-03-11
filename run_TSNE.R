# Name: Kiku Koyano
# Date: 05/2021
# Description: run tsne and pca on TPM matrix (or whichever matrix, values will be log2 transformed). 
# also correlate samples by TPM values to check for sample reproducibility

helper_function_file="/Users/kiku/mount/scripts/R_scripts/helper_functions.R"
source(helper_function_file)

libraries = c("data.table", "dplyr", "optparse", "ggpubr", "tidyr", "ggplot2", 
              "ggthemes", "purrr", "forcats", "RColorBrewer", "grid", "gridExtra", "wesanderson",
              "psych", "corrplot", "Rtsne", "ggfortify", "GGally") 
load_libs = lapply(libraries, require, quietly = T, character.only = TRUE) 


option_list <- list(make_option(c("-i", "--input_file"), action = "store", type = "character", default = NULL,
                                help = "TPM matrix (rows are genes, cols are sample)"),
                    make_option(c("-o", "--output_dir"), action = "store", type = "character",
                                help = "location of output directory for files to go to"),
                    make_option(c("-m", "--metadata_fn"), action = "store", type = "character",
                                help = "metadata file, needs to have at least 2 columns 1) sample name and 2) biofluid")
                    
                    
                    
)

opt <- parse_args(OptionParser(option_list=option_list))
input_file = opt$input_file
output_dir = opt$output_dir
metadata_fn = opt$metadata_fn

# testing variables 
# input_file = "/Users/kiku/mount/projects/long_csf_2020_03/analysis/correlate_replicates/read_count-20/TPM_matrix.txt"
# output_dir = "/Users/kiku/mount/projects/long_csf_2020_03/analysis/correlate_replicates/read_count-20/"
# metadata_fn = "/Users/kiku/mount/projects/long_csf_2020_03/data/metadata/metadata_mapping_info.txt"

# tpm_df = fread("/Users/kiku/mount/projects/long_csf_2020_03/analysis/correlate_replicates/read_count-20/TPM_matrix.txt", header = T)


# run tsne for these biotypes separately
biotype_list = c("protein_coding", "lincRNA")

tpm_df = fread(input_file, header = T )
x = lapply(biotype_list , function(biotype){
  
  tpm_protein_coding = tpm_df %>% filter(gene_type == !! biotype)
  tpm_protein_coding$gene_type = NULL 
  row.names(tpm_protein_coding) = tpm_protein_coding$gene_name
  tpm_protein_coding$gene_name= NULL
  tpm_protein_coding = tpm_protein_coding + 1 #pseudocount
  tpm_protein_coding = log2(tpm_protein_coding) #log2 TPMs 
  
  # make correlation_df 
  correlation_df = lapply(seq(colnames(tpm_protein_coding)), function(i){
    sample1_name = colnames(tpm_protein_coding)[i] 
    message(sample1_name)
    sample_df = data.frame()
    for (j in seq(ncol(tpm_protein_coding)) ){
      sample2_name = names(tpm_protein_coding)[j]
      # message(sample2_name)
      tpm_protein_coding.mat = as.matrix(tpm_protein_coding)
      corr_df = data.frame("V1" = tpm_protein_coding.mat[,i], "V2" = tpm_protein_coding.mat[,j] ) 
      corr_df.filtered = corr_df[apply(corr_df, 1, function(x) !all(x==0)),]
      
      corr_value = cor(corr_df.filtered$V1, corr_df.filtered$V2)
      result = data.frame("sample1_name" = corr_value, sample_name = sample2_name)
      sample_df = bind_rows(sample_df, result)
      
    }
    colnames(sample_df)[1] <- sample1_name
    sample_df
  })
  correlation_matrix = Reduce(function(x,y, ...) merge(x, y, by = "sample_name"), 
                              correlation_df )
  
  # for each correlation of the log2(TPM) between each sample, remove values from genes that are both 0 between both samples 
  row.names(correlation_matrix) = correlation_matrix$sample_name
  
  correlation_matrix$sample_name = NULL
  correlation_matrix = correlation_matrix[row.names(correlation_matrix)]

  
  # TSNE of samples based on the log2(TPM)
  tsne_matrix <- normalize_input(as.matrix(t(tpm_protein_coding))) #for tsne, the rows need to be samples 
  tsne_out = Rtsne(tsne_matrix, check_duplicates = T)
  tsne_plot <- data.frame(x = tsne_out$Y[,1],
                          y = tsne_out$Y[,2])
  my_gene_row = metadata_df %>% filter(sample_name %in% row.names(tsne_matrix)) 
  tsne_df = cbind(tsne_plot, my_gene_row)
  plot_title = paste(biotype, "TSNE on log2(TPM)")
  
  tsne_plot = ggplot(tsne_df) + 
    geom_point(aes(x=x, y=y, color=biofluid), size = 3) + 
    labs(title =plot_title) + theme_classic2(base_size = 10)
  
  # PCA on the log2(TPM)
  pca_obj = prcomp(x = t(tpm_protein_coding) ) 
  pca_plot_title = paste(biotype, "PCA on log2(TPM)")
  pca.plot <- autoplot(pca_obj, data = my_gene_row, colour = 'biofluid', size = 3 ) + 
    labs(title =pca_plot_title) + theme_classic2(base_size = 10)
  
  # correlating samples
  filename = paste("corr_samples_matrix.", biotype, ".pdf", sep ="")
  corrplot1 = corrplot(as.matrix(correlation_matrix), method="number", type="upper")
  corrplot2 = ggpairs(tpm_protein_coding, lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1)))
  output_plots(output_dir, filename, plot_list = list(corrplot1, corrplot2) , width = 12, height = 10)
  
  # output plots 
  filename = paste( "tsne_pca.log2TPM.", biotype, ".pdf", sep ="")
  output_plots(output_dir, filename, plot_list = list(tsne_plot, pca.plot) , width = 4, height = 3)

  
  
})

message("Job Completed")

