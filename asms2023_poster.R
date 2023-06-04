library(promor)
library(VennDiagram)
library(ggplot2)

################################################################################

# 1. Benchmarking

################################################################################
################################################################################

# Perseus data parsing

################################################################################

# Upload Perseus results
perseus_results <- read.csv("https://raw.githubusercontent.com/caranathunge/promor_bioRxiv_preprint/main/Perseus_MinProb_noNormalization_DEresults.txt",
                            sep = "\t"
)


# Reduce the data frame to significant hits and limit the data frame to
#only those columns we need
de_perseus <- perseus_results[perseus_results$H_vs_L_Significant == "+", c(
  "Majority.protein.IDs",
  "H_vs_L_P.Value",
  "H_vs_L_logFC"
)]


# Add a Protein.IDs column to promor results by extracting the first protein
#from majority_protein_ids
de_perseus$Protein.IDs <- sapply(strsplit(as.character
                                          (de_perseus$Majority.protein.IDs),
                                          ";"), "[", 1)

# remove maj prot id column
de_perseus <- subset(de_perseus, select = -Majority.protein.IDs)

# Add a new column with the name of the method used
de_perseus$method <- "Perseus"

# Let's give both data frames similar column names
colnames(de_perseus) <- c("p_val", "log_fc", "protein", "method")

# Make a list object to build a venn diagram
de_perseus_prot <- de_perseus$protein


################################################################################

# promor analysis

################################################################################

# Create a raw_df object with proteinGroups.txt and exp_design file
raw_df <- create_df(
  prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
  exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt",
  uniq_pep = 2
)


# Filter out proteins with higher than 40% missing data in either one group
#(in other words - requires 60% valid data in each group to retain the protein)
raw_df_filt <- filterbygroup_na(raw_df, set_na = 0.40, filter_condition = "either")


# Impute missing data
imp_df <- impute_na(raw_df_filt, seed = 327, method = "minProb")

# Find DE proteins
fit_df <- find_dep(imp_df)

# Save results from all DE proteins
fit_df <- find_dep(imp_df,
                   n_top = 1294,
                   save_tophits = TRUE,
                   save_output = TRUE,
                   file_path = "./")

# Upload the TopHits
de_promor <- read.csv("./TopHits.txt", sep = "\t")

# Add a Protein.IDs column to promor results by extracting the first protein
#from majority_protein_ids
de_promor$Protein.IDs <- sapply(strsplit
                                (as.character(de_promor$majority_protein_id),
                                  ";"), "[", 1)

# Extract only those columns we need from de_promor
de_promor <- de_promor[, c("Protein.IDs", "P.Value", "logFC")]

# Add a new column with the method information
de_promor$method <- "promor"

# Let's give both data frames similar column names
colnames(de_promor) <- c("protein", "p_val", "log_fc", "method")

# Make a list object to build a venn diagram
de_promor_prot <- de_promor$protein


################################################################################

# Figure 1A. Comparsion - Number of DE proteins

################################################################################


venn.diagram(list("promor" = de_promor_prot, "Perseus" = de_perseus_prot),
             fill = c("#17456B", "#ACF0F2"),
             alpha = c(0.5, 0.5),
             resolution = 400,
             lwd = 5,
             filename = "./venn_diagram.tiff",
             scaled = TRUE,
             ext.pos = 0,
             ext.percent = 0.5,
             fontface = "bold",
             ext.line.lwd = 3
)


################################################################################

# Figure 1B. Comparsion - log FC

################################################################################

# combine data from both dataframes into one
df_all <- merge(de_promor, de_perseus, by = "protein")

# Convert non-numeric values to numeric
df_all$log_fc.y <- as.numeric(df_all$log_fc.y)
df_all$p_val.y <- as.numeric(df_all$p_val.y)
attach(df_all)

#calculate and annotate pearson correlation
grob1 <- grobTree(textGrob(paste("Pearson Correlation : ",
                                 round(cor(log_fc.x, log_fc.y), 4)),
                           x = 0.5, y = 0.97, hjust = 0,
                           gp = gpar(
                             col = "black",
                             fontsize = 15,
                             fontface = "bold"
                           )
))

#Make the plot
ggplot(df_all, aes(x = log_fc.x, y = log_fc.y)) +
  geom_point(size = 10, shape = 20, alpha = .4, col = "#17456B") +
  # ggtitle(bquote('promor vs Perseus - log' [2]~ 'fold change')) +
  geom_smooth(method = lm, se = FALSE, lwd = 1) +
  xlab(bquote("promor log"[2] ~ "fold-change")) +
  ylab(bquote("Perseus log"[2] ~ "fold-change")) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 2),
    axis.line.x = element_line(color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  annotation_custom(grob1)


################################################################################

# Figure 1C. Comparsion - p-value

################################################################################

#calculate and annotate pearson correlation
grob2 <- grobTree(textGrob(paste("Pearson Correlation : ",
                                 round(cor(log(p_val.x), log(p_val.y)), 4)),
                           x = 0.5, y = 0.97, hjust = 0,
                           gp = gpar(
                             col = "black",
                             fontsize = 15,
                             fontface = "bold"
                           )
))

#Make the plot
ggplot(df_all, aes(x = log(p_val.x), y = log(p_val.y))) +
  geom_point(size = 10, shape = 20, alpha = .4, col = "#17456B") +
  # ggtitle(bquote('promor vs Perseus - log' [10]~ 'p-value')) +
  geom_smooth(method = lm, se = FALSE, lwd = .5) +
  xlab(bquote("promor log"[10] ~ "p-value")) +
  ylab(bquote("Perseus log"[10] ~ "p-value")) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 2),
    axis.line.x = element_line(color = "black"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  annotation_custom(grob2)

################################################################################

# 2. Quality control & Visualization

################################################################################

################################################################################

# Figure 2A. Missing Data Heatmap

################################################################################
#Upload the data
raw_df1 <- create_df(prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg2.txt",
                     exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed2.txt",
                     uniq_pep = 1,
                     tech_reps = TRUE)

#Calculate average across tech reps
rawdf1_avg <- aver_techreps(raw_df1)

#Filter proteins by group level missing data
rawdf1_filt <- filterbygroup_na(rawdf1_avg, set_na = 0.50, filter_condition = "each")

#Make missing data heatmap
heatmap_na(rawdf1_filt, text_size = 15, reorder_y = TRUE,  save = TRUE, file_path = ".", file_type = "png", dpi = 300, palette = "mako")

################################################################################

# Figure 2B. Density plots to visualize the impact of imputation

################################################################################
#Impute missing data using minDet method
imp_df1 <- impute_na(rawdf1_filt, method = "minDet", seed = 327)

#Visualize missing data imputation
impute_plot(original = rawdf1_filt, imputed = imp_df1,
            global = FALSE,  n_col = 2, n_row = 3,
            dpi = 300, save = TRUE, file_path = ".",
            file_type = "png", palette = "mako", text_size = 20)

################################################################################

# Figure 2C. Density plots to visualize the impact of normalization

################################################################################

#Normalize the data set
norm_df1 <- normalize_data(imp_df1, method = "quantile")

#Make density plots to compare before and after
norm_plot(original = imp_df1, normalized = norm_df1,
          type = "density", save = TRUE, file_path = ".", dpi = 300,
          file_type = "png", palette = "mako", text_size = 30)


################################################################################

# Figure 2D. Scatter plots to visualize correlationbetween pairs of tech.replicates

################################################################################
#Make correlation plots
corr_plot(raw_df1, rep_1 = 1, rep_2 = 2,
          file_type = "png", save = TRUE, file_path = ".",
          dpi = 300, text_size = 20,
          n_row = 3, n_col = 2,
          palette = "mako")

################################################################################

# 3. Differential expression analysis

################################################################################

################################################################################

# Figure 3A. Volcano plot

################################################################################

#Upload data set 2
raw_df2 <- create_df(prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
                     exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt",
                     uniq_pep = 2)

#Filter by group level missing data
rawdf2_filt <- filterbygroup_na(raw_df2, set_na = 0.34, filter_condition = "each")

#Impute missing data
imp_df2 <- impute_na(rawdf2_filt, method = "kNN", seed = 327)

#Normalize data
norm_df2 <- normalize_data(imp_df2)

#Find DE proteins
fit_df2 <- find_dep(norm_df2)

#Make volcano plot
volcano_plot(fit_df2, save = TRUE, file_path = ".",
             file_name = "volcano_plot_ecoli", dpi = 300,
             file_type = "png", palette = "mako", text_size = 15)


################################################################################

# Figure 3B. Heatmap of DE proteins

################################################################################
#Upload data set 3
raw_df3 <- create_df(prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg3.txt",
                     exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed3.txt",
                     uniq_pep = 2)

#Filter by group level missing data
rawdf3_filt <- filterbygroup_na(raw_df3, set_na = 0.34, filter_condition = "each")

#Impute missing data
imp_df3 <- impute_na(rawdf3_filt, method = "kNN", seed = 327)

#Normalize data
norm_df3 <- normalize_data(imp_df3)

#Find DE proteins
fit_df3 <- find_dep(norm_df3, cutoff = 0.1)

#Make a heatmap of top 20 DE proteins
heatmap_de(fit_df = fit_df3, df = norm_df3,
           n_top = 20, cutoff = 0.1,
           save = TRUE, file_path = ".", dpi = 300, file_name = "heatmap_covid",
           file_type = "png", palette = "mako")

################################################################################

# 4. Feature selection

################################################################################

################################################################################

# Figure 4A. Visualize feature variation - density plots

################################################################################

#create a model_df object
model_df3 <- pre_process(fit_df = fit_df3, norm_df = norm_df3, sig_cutoff = 0.06)

#Make feature plots
feature_plot(model_df3,  save = TRUE,
             type = "density", dpi = 300, file_name = "feature_covid_density",
             file_path = ".", file_type = "png",  n_col = 2, n_row = 3,
             plot_width = 3, plot_height = 15,
             palette = "mako", text_size = 20)


################################################################################

# Figure 4B. Variable importance plots

################################################################################

#split the model_df object into training and test data
split_df3 <- split_data(model_df3, train_size = 0.5, seed = 8314)

#train models on training data
model_list <- train_models(split_df3, resample_method = "repeatedcv", seed = 351)

#Make variable importance plots
varimp_plot(model_list, save = TRUE,
            plot_width = 28, plot_height = 20,
            n_col = 2, n_row = 2 ,
            text_size = 10, dpi = 300, file_path = ".", file_type = "png",
            file_name = "varimp_covid", palette = "mako")

################################################################################

# 5. Model Building & Evaluation

################################################################################

################################################################################

# Figure 5A. Performance plots

################################################################################

performance_plot(model_list, type = "dot",
                 dpi = 300, file_name = "covid_performance",
                 file_type = "png", save = TRUE, file_path = ".", palette = "mako",
                 text_size = 20, plot_height = 5)



################################################################################

# Figure 5B. ROC plots

################################################################################
#Test the models on the test data
prob_list <- test_models(model_list = model_list, split_df = split_df3)

#Make ROC curves
roc_plot(probability_list = prob_list, split_df = split_df3,
         save = TRUE, file_path = ".", file_name = "covid_roc", file_type = "png",
         dpi = 300, palette = "mako",
         plot_height = 14, plot_width = 14, text_size = 10)
