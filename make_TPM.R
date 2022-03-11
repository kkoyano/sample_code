# Name: Kiku Koyano
# Date: 05/2021
# Description: covert counts to TPM for each RNA-seq sample 
# input: 1) list of all samples and their respecive gene count files to inlcude in TPM calculation and normalization 
#        2) read count minimum 
#        3) file containing gene lengths 
#        4) metadata file with the sample name and biofluid/tissue type

# output: 
# TPM_matrix.txt 
  # gene_name	gene_type	Ascites_1	Amniotic_fluid_1 .... N_samples
  # A1BG	protein_coding	19.4599881098698	0 ....
  # A1CF	protein_coding	0	9.45522454848814  ....
  # A2M	protein_coding	55.3669595749306 ....

# TPM_matrix.txt 
# gene_name	gene_type	sample_name TPM 
# A1BG	protein_coding	Ascites_1 19.4599881098698	0
# A1CF	protein_coding	Ascites_1 0	9.45522454848814
# A2M	protein_coding	Ascites_1 55.3669595749306

libraries = c("data.table", "dplyr", "optparse", "ggpubr", "tidyr", "ggplot2", 
              "ggthemes", "purrr", "forcats", "RColorBrewer", "grid", "gridExtra", "wesanderson",
              "psych", "corrplot", "Rtsne", "ggfortify", "GGally") 
load_libs = lapply(libraries, require, quietly = T, character.only = TRUE) 


option_list <- list(make_option(c("-i", "--input_file"), action = "store", type = "character", default = NULL,
                                help = "input file with 2 columns 1) sample_name 2) filenames of gene counts *.GeneFeatureCount file."),
                    make_option(c("-o", "--output_dir"), action = "store", type = "character",
                                help = "location of output directory for files to go to"),
                    make_option(c("-m", "--metadata_fn"), action = "store", type = "character",
                                help = "metadata file, needs to have at least 2 columns 1) sample name and 2) biofluid"),
                    make_option(c("-r", "--read_count"), action = "store", type = "numeric", default = 20, 
                                help = "total read count threshold each gene"),
                    make_option(c("-g", "--gene_lengths_file"), action = "store", type = "character",
                                help = "contains the length of genes to calculate TPM:
                                ~/mount/annotations/hg19/genes/ensembl/gene_lengths_summary.txt")
                  
)

opt <- parse_args(OptionParser(option_list=option_list))
input_file = opt$input_file
output_dir = opt$output_dir
metadata_fn = opt$metadata_fn
READ_COUNT = opt$read_count
gene_lengths_fn = opt$gene_lengths_file

# Developer variables 
# input_file = "/Users/kiku/mount/projects/long_csf_2020_03/analysis/correlate_replicates/correlate_replicates.input.txt"
# output_dir = "/Users/kiku/mount/projects/long_csf_2020_03/analysis/correlate_replicates/read_count-20/"
# READ_COUNT=20
# metadata_fn = "/Users/kiku/mount/projects/long_csf_2020_03/data/metadata/metadata_mapping_info.txt"
# gene_lengths_fn = "~/mount/annotations/hg19/genes/ensembl/gene_lengths_summary.txt"


# get TPM 
gene_lengths_df = fread(gene_lengths_fn, header = T ) %>% dplyr::select(gene_name, gene_biotype, tx_len_mean)
colnames(gene_lengths_df)[2] = "gene_type"

input_filenames = fread(input_file, header = F, sep = "\t")
# input_filenames$V2 = gsub("/u/home/k/kkoyano/gxxiao3", "~/mount", input_filenames$V2)

metadata_df = fread(metadata_fn, header = T ) %>% dplyr::select(sample_name, biofluid)
biol_groups = unique(metadata_df$biofluid)

# Get the TPM and RPKM of each gene
df_list = lapply(input_filenames$V1, function(filename){
  sample_name = filename
  gene_count_file = input_filenames %>% filter(V1 == !! sample_name)
  gene_count_file = gene_count_file$V2
  gene_count_df = fread(gene_count_file, header = T ) %>% dplyr::select(gene_name, gene_type, gene_count) 
  gene_count_df.lengths = inner_join(gene_count_df, gene_lengths_df, by = c("gene_name", "gene_type") ) %>% 
    filter(gene_count >= READ_COUNT)
  
  # TPM step1: divide all read counts by length of gene 
  gene_count_df.lengths$RKB =  gene_count_df.lengths$gene_count/(gene_count_df.lengths$tx_len_mean/1000)
  # TPM step2: sum all RKB and divide by 1,000,000
  tpm_scale_factor =  sum(gene_count_df.lengths$RKB)/1000000
  gene_count_df.lengths$TPM = gene_count_df.lengths$RKB/tpm_scale_factor
  
  # FPKM step1 = count total reads in sample and divide by 1000000
  fpkm_scale = sum(gene_count_df.lengths$gene_count)/1000000
  gene_count_df.lengths$fpm = gene_count_df.lengths$gene_count/fpkm_scale
  gene_count_df.lengths$FPKM = gene_count_df.lengths$fpm/(gene_count_df.lengths$tx_len_mean/1000)
  gene_count_df.lengths$RKB = NULL
  gene_count_df.lengths$fpm = NULL
  gene_count_df.lengths$sample_name = sample_name
  gene_count_df.lengths
})


# TPM data frame 
tpm_list = lapply(df_list, function(df){
  sample_name = df$sample_name %>% unique()
  df.select = df %>% dplyr::select(gene_name, gene_type, TPM)
  colnames(df.select)[3] = sample_name
  # df.select$gene_name_type = paste(df.select$gene_name, df.select$gene_type, sep = '.')
  df.select
})
tpm_df = Reduce(function(...) merge(..., all = TRUE), tpm_list)  


tpm_df[is.na(tpm_df)] <- 0 

output_fn = paste(output_dir, "TPM_matrix.txt", sep = "")
fwrite(file = output_fn, x = tpm_df, sep = "\t", quote = F)

tpm_df.long = bind_rows(tpm_list)
tpm_df.long[is.na(tpm_df.long)] <- 0 
output_fn = paste(output_dir, "TPM_matrix.long.txt", sep = "")
fwrite(file = output_fn, x = tpm_df.long, sep = "\t", quote = F)
rm(tpm_df.long)




message("Job Completed")
