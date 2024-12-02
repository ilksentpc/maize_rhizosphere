##alpha
setwd("~/......")

library(ggplot2)
library(phyloseq)
library(vegan)

otu_leafbac <- read.csv("otu_table_soil.csv", row.names = 1)

tax_leafbac <- read.csv("taxa_soil.csv", row.names = 1)

otumat <- as.matrix(otu_leafbac)

taxmat <- as.matrix(tax_leafbac)

metadata <- read.csv("metadata_rhizo_only.csv", row.names = 1)
sampledata = sample_data(data.frame(metadata))


OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU

physeq = phyloseq(OTU, TAX, sampledata)
physeq

ASV_TREE = read_tree("bacterial_soil.tre")

physeq = phyloseq(OTU, TAX, sampledata, ASV_TREE)
physeq

physeq = phyloseq(OTU, TAX, sampledata)
physeq

library(Biostrings)
seqtab<-readDNAStringSet("bacterialsoil_tree_file.fasta", format="fasta")
physeq = phyloseq(OTU, TAX, sampledata, ASV_TREE, seqtab)
physeq

library(microeco)
library(file2meco)

# from phyloseq to microtable object
meco_dataset <- phyloseq2meco(physeq)
meco_dataset
library(tidyverse)
library(magrittr)
library(ggplot2)
meco_dataset$sample_table$Genotype  %<>% factor(levels = c("Perennial teosinte", "Balsas teosinte", "Mexican landrace", "Mexican inbred", "US landrace", "US inbred"))

meco_dataset$sample_table$Plant.Group  %<>% factor(levels = c("Teosinte", "Mexican maize", "US maize"))

# Cumulative Sum Scaling (CSS)
dt <- trans_norm$new(meco_dataset)
dt2 <- dt$norm(method = "CSS")
View(dt2$otu_table)

t1 <- trans_alpha$new(dataset = dt2, group = "Plant.Group")
# return t1$data_stat
head(t1$data_stat)
t1$data_alpha

t1$data_alpha$Plant.Group  %<>% factor(levels = c("Teosinte", "Mexican maize", "US maize"))

t1$data_alpha$Genotype <- factor(
  t1$data_alpha$Genotype,
  levels = c("Teosinte", "Mexican maize", "US maize")
)
##
shannon_plot <- t1$plot_alpha(
  measure = "Shannon",       # Specify Shannon diversity
  group = "Plant.Group",     # Group by plant group
  use_boxplot = TRUE,        # Add boxplot
  point_alpha = 0.5,         # Add jitter points with some transparency
  add_sig_text_size = 6,     # Size of significance annotations
  xtext_angle = 0           # Rotate x-axis text for better visibility
)

png("rhizospehere_plantgroup.png", width=13, height=8, units="in", res=300)
shannon_plot + theme(
  text = element_text(family = "Arial", size = 40),
  axis.title = element_text(size = 24),
  axis.text = element_text(size = 12),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 16)
)

dev.off()



anova_model <- aov(Value ~ Plant.Group, data = t1$data_alpha)
tukey_results <- TukeyHSD(anova_model)
print(tukey_results)
plot(tukey_results)


#plant group comparisons
c1 <- c(1, -1, 0)  # Perennial teosinte vs Balsas teosinte
c2 <- c(0, 1, -1)  # Balsas teosinte vs Mexican maize

contrast_matrix <- cbind(c1, c2)
contrasts(t1$data_alpha$Plant.Group) <- contrast_matrix
contrast_model <- aov(Value ~ Plant.Group, data = t1$data_alpha)
summary(contrast_model)


summary.aov( contrast_model,
split = list(
  Plant.Group = list(
      "Perennial vs Balsas" = 1,
      "Balsas vs Mexican" = 2)))

pairwise.t.test(
  t1$data_alpha$Value,
  t1$data_alpha$Plant.Group,
  p.adjust.method = "fdr"
)


shannon_data <- t1$data_alpha %>%
  filter(Measure == "Shannon") %>%
  dplyr::select(Plant.Group, Value) # Use dplyr::select explicitly
shannon_data$Plant.Group <- factor(shannon_data$Plant.Group)
shannon_aov <- aov(Value ~ Plant.Group, data = shannon_data)
summary(shannon_aov)





#####Genotype within PLant group

t2 <- trans_alpha$new(dataset = dt2, group = "Genotype")
head(t2$data_stat)

t2$data_alpha$Genotype <- factor(
  t2$data_alpha$Genotype,
  levels = c("Perennial teosinte", "Balsas teosinte", "Mexican landrace", "Mexican inbred", "US landrace", "US inbred")
)


shannon_plot2 <- t2$plot_alpha(
  measure = "Shannon",       # Specify Shannon diversity
  group = "Genotype",     # Group by plant group
  use_boxplot = TRUE,        # Add boxplot
  point_alpha = 0.5,         # Add jitter points with some transparency
  add_sig_text_size = 6,     # Size of significance annotations
  xtext_angle = 45           # Rotate x-axis text for better visibility
)


png("rhizospehere_genotype_shannon.png", width=13, height=8, units="in", res=300)
shannon_plot2 + theme(
  text = element_text(family = "Arial", size = 40),
  axis.title = element_text(size = 24),
  axis.text = element_text(size = 12),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 16)
)

dev.off()



t2$data_alpha$Plant.Group <- as.factor(t2$data_alpha$Plant.Group)
t2$data_alpha$Genotype <- as.factor(t2$data_alpha$Genotype)

#nested
nested_aov <- aov(Value ~ Plant.Group / Genotype, data = t2$data_alpha)
summary(nested_aov)


#specific comparisons
c1 <- c(1, -1, 0, 0, 0, 0)  # Perennial teosinte vs Balsas teosinte
c2 <- c(0, 1, -1, 0, 0, 0)  # Balsas teosinte vs Mexican landrace
c3 <- c(0, 0, 1, -1, 0, 0)  # Mexican landrace vs Mexican inbred
c4 <- c(0, 0, 1, 0, -1, 0)  # Mexican landrace vs US landrace
c5 <- c(0, 0, 0, 0, 1, -1)  # US landrace vs US inbred

contrast_matrix <- cbind(c1, c2, c3, c4, c5)
contrasts(t2$data_alpha$Genotype) <- contrast_matrix
contrast_model <- aov(Value ~ Genotype, data = t2$data_alpha)
summary(contrast_model)

summary.aov(
  contrast_model,
  split = list(
    Genotype = list(
      "Perennial vs Balsas" = 1,
      "Balsas vs MX Landrace" = 2,
      "MX Landrace vs MX Inbred" = 3,
      "MX Landrace vs US Landrace" = 4,
      "US Landrace vs US Inbred" = 5
    )
  )
)

pairwise.t.test(
  t2$data_alpha$Value,
  t2$data_alpha$Genotype,
  p.adjust.method = "fdr"
)



####for richness, don't do normalization 

t3 <- trans_alpha$new(dataset = meco_dataset, group = "Plant.Group")
head(t3$data_stat)


t3$data_alpha$Plant.Group  %<>% factor(levels = c("Teosinte", "Mexican maize", "US maize"))

t3$data_alpha$Genotype <- factor(
  t3$data_alpha$Genotype,
  levels = c("Teosinte", "Mexican maize", "US maize")
)


shannon_plot3 <- t3$plot_alpha(
  measure = "Chao1",       
  group = "Plant.Group",     
  use_boxplot = TRUE,        
  point_alpha = 0.5,         
  add_sig_text_size = 6,  
  xtext_angle = 0          
)


png("rhizospehere_plantgroup_chao1.png", width=13, height=8, units="in", res=300)
shannon_plot3 + theme(
  text = element_text(family = "Arial", size = 40),
  axis.title = element_text(size = 24),
  axis.text = element_text(size = 12),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 16)
)

dev.off()



anova_model3 <- aov(Value ~ Plant.Group, data = t3$data_alpha)
tukey_results3 <- TukeyHSD(anova_model3)
print(tukey_results3)
plot(tukey_results3)


#####Genotype within PLant group_chao1

t4 <- trans_alpha$new(dataset = meco_dataset, group = "Genotype")
head(t4$data_stat)

t4$data_alpha$Genotype <- factor(
  t4$data_alpha$Genotype,
  levels = c("Perennial teosinte", "Balsas teosinte", "Mexican landrace", "Mexican inbred", "US landrace", "US inbred")
)

shannon_plot4 <- t4$plot_alpha(
  measure = "Chao1",      
  group = "Genotype",    
  use_boxplot = TRUE,     
  point_alpha = 0.5,     
  add_sig_text_size = 6,  
  xtext_angle = 45   
)


png("rhizospehere_genotype_chao1.png", width=13, height=8, units="in", res=300)
shannon_plot4 + theme(
  text = element_text(family = "Arial", size = 40),
  axis.title = element_text(size = 24),
  axis.text = element_text(size = 12),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 16)
)

dev.off()


t4$data_alpha$Plant.Group <- as.factor(t4$data_alpha$Plant.Group)
t4$data_alpha$Genotype <- as.factor(t4$data_alpha$Genotype)

#nested
nested_aov4 <- aov(Value ~ Plant.Group / Genotype, data = t4$data_alpha)
summary(nested_aov4)

pairwise.t.test(
  t4$data_alpha$Value,
  t4$data_alpha$Genotype,
  p.adjust.method = "fdr"
)


#Beta Diversity
#beta:CSS

setwd("~/.......")

library(ggplot2)
library(phyloseq)
library(vegan)

otu_leafbac <- read.csv("otu_table_soil.csv", row.names = 1)

tax_leafbac <- read.csv("taxa_soil.csv", row.names = 1)

otumat <- as.matrix(otu_leafbac)

taxmat <- as.matrix(tax_leafbac)

metadata <- read.csv("metadata_rhizo_only.csv", row.names = 1)
sampledata = sample_data(data.frame(metadata))


OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU

physeq = phyloseq(OTU, TAX, sampledata)
physeq

ASV_TREE = read_tree("bacterial_soil.tre")

physeq = phyloseq(OTU, TAX, sampledata, ASV_TREE)
physeq

physeq = phyloseq(OTU, TAX, sampledata)
physeq

library(Biostrings)

seqtab<-readDNAStringSet("bacterialsoil_tree_file.fasta", format="fasta")

physeq = phyloseq(OTU, TAX, sampledata, ASV_TREE, seqtab)
physeq

library(microeco)
library(file2meco)

# from phyloseq to microtable object
meco_dataset <- phyloseq2meco(physeq)
meco_dataset

library(tidyverse)
library(magrittr)
library(ggplot2)

meco_dataset$sample_table$Genotype  %<>% factor(levels = c("Perennial teosinte", "Balsas teosinte", "Mexican landrace", "Mexican inbred", "US landrace", "US inbred"))

meco_dataset$sample_table$Plant.Group  %<>% factor(levels = c("Teosinte", "Mexican maize", "US maize"))

meco_dataset$sample_table$geno.pcoa  %<>% factor(levels = c("Perennial teosinte", "Balsas teosinte", "Mexican landrace", "Mexican inbred", "US landrace", "US inbred", "GOS", "LNC", "RYD", "Hp301", "B73", "Mo17", "W438"))

dt <- trans_norm$new(meco_dataset)
dt2 <- dt$norm(method = "CSS")
View(dt2$otu_table)

dt2$cal_betadiv(method = "bray", unifrac = FALSE, binary = FALSE)

t1 <- trans_beta$new(dataset = dt2, group = "geno.pcoa", measure = "bray")


beta_div_matrix <- dt2$beta_diversity
print(beta_div_matrix)

t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)
t1$plot_ordination(plot_color = "geno.pcoa", plot_shape = "geno.pcoa", plot_type = c("point", "ellipse"))


axis_title_size <- 20
custom_text_size <- 16

ordination_plot <- t1$plot_ordination(
  plot_color = "geno.pcoa", 
  plot_shape = "geno.pcoa", 
  plot_type = c("point", "ellipse")
) + 
  
  theme(
    text = element_text(family = "Arial"),          
    legend.text = element_text(size = custom_text_size, family = "Arial"), 
    axis.text = element_text(size = 15, family = "Arial"),
    axis.title = element_text(size = axis_title_size, family = "Arial"), 
    strip.text = element_text(size = 15, family = "Arial") 
  )  +
  
  guides(
    color = guide_legend(override.aes = list(size = 7)), 
    shape = guide_legend(override.aes = list(size = 7))  
  )

png("rhizospehere_genotype_bray_curtis.png", width=13, height=8, units="in", res=300)
print(ordination_plot)
dev.off()

library(vegan)

bray_curtis_matrix <- beta_div_matrix$bray
class(bray_curtis_matrix)  # Should return "matrix"
is.numeric(bray_curtis_matrix)  # Should return TRUE


nested_formula <- as.formula("bray_curtis_matrix ~ Plant.Group / geno.pcoa")
sample_df <- dt2$sample_table

#nested PERMANOVA
nested_adonis_result <- adonis2(
  nested_formula,
  data = sample_df,
  permutations = 999
)

print(nested_adonis_result)

library(pairwiseAdonis)

#Loop 
pairwise_results <- lapply(unique(sample_df$Plant.Group), function(group) {
  # Subset the data for the current Plant.Group
  subset_indices <- sample_df$Plant.Group == group
  subset_dist_matrix <- bray_curtis_matrix[subset_indices, subset_indices]
  subset_sample_df <- sample_df[subset_indices, ]
  
  pairwise_result <- pairwise.adonis(
    x = subset_dist_matrix,
    factors = subset_sample_df$geno.pcoa,
    p.adjust.m = "fdr"  # Use False Discovery Rate adjustment
  )
  
  pairwise_result$Plant.Group <- group
  return(pairwise_result)
})

pairwise_results_combined <- do.call(rbind, pairwise_results)
print(pairwise_results_combined)

specific_comparisons <- c(
  "Perennial teosinte vs Balsas teosinte",
  "Balsas teosinte vs Mexican landrace",
  "Mexican landrace vs Mexican inbred",
  "GOS vs RYD",
  "LNC vs RYD",
  "Hp301 vs B73",
  "Hp301 vs W438",
  "Mo17 vs W438",
  "Mo17 vs B73"
)

filtered_results <- pairwise_results_combined[
  pairwise_results_combined$pairs %in% specific_comparisons, 
]


print(filtered_results)

global_pairwise_results <- pairwise.adonis(
  x = bray_curtis_matrix,
  factors = sample_df$geno.pcoa,
  p.adjust.m = "fdr"  # Use False Discovery Rate adjustment
)


specific_comparisons <- c(
  "Perennial teosinte vs Balsas teosinte",
  "Balsas teosinte vs Mexican landrace",
  "Mexican landrace vs Mexican inbred",
  "GOS vs RYD",
  "LNC vs RYD",
  "Hp301 vs B73",
  "Hp301 vs W438",
  "Mo17 vs W438",
  "Mo17 vs B73"
)


filtered_results <- global_pairwise_results[
  global_pairwise_results$pairs %in% specific_comparisons, 
]

print(filtered_results)
