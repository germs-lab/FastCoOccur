

# read phyloseq object
phyloseq_object <- readRDS("sample_phyloseq_object.rds")

trt <- as.character(sample_data(phyloseq_object)$Set)

f_table <- cbind(trt, t(otu_table(phyloseq_object)))

# write csv file
write.csv(f_table, "sample_input_file.csv")