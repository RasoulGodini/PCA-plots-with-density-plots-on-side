
library(dplyr)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggside)
library(ggpubr)
library(grid)



#Random data generation
set.seed(123)  # For reproducibility

n_participants <- 200
n_genes <- 500

sample_id <- sample(paste0("Sample_", sprintf("%03d", 1:200)), n_participants, replace = FALSE)
sex <- sample(c("Male", "Female"), n_participants, replace = TRUE)

# Initialize matrix for gene expression
gene_data <- matrix(nrow = n_participants, ncol = n_genes)

# Generate gene expression values based on Sex
for (i in 1:n_participants) {
  if (sex[i] == "Male") {
    gene_data[i, ] <- rnorm(n_genes, mean = 5, sd = 3.6)
  } else {
    gene_data[i, ] <- rnorm(n_genes, mean = 6, sd = 4)
  }
}

# Convert to dataframe and add Sex column
colnames(gene_data) <- paste0("Gene_", 1:n_genes)  # Name genes
df <- data.frame(Sampe_ID = sample_id, Sex = sex, gene_data)  # Make the datframe


res.pca <- PCA(df[,3:ncol(df)], graph = FALSE, scale.unit = TRUE, ncp = 5)


pca_data <- as.data.frame(res.pca$ind$coord)
pca_data_annot <- merge(df[,1:2], pca_data, by = "row.names")


# To get the value of sex status
condition <- factor(pca_data_annot$Sex)

#To get median for density plots
pca_data_median <- pca_data_annot %>% group_by(Sex) %>% summarize_at(c("Dim.1", "Dim.2", "Dim.3", "Dim.4", "Dim.5"), median)


# Perform t-test (two-sided) and calculate p-value to compare eigenvectors
eigenve_pca_data <- as.data.frame(res.pca$ind$coord)


pc_list <- c("Dim.1", "Dim.2")
# To store the results
p_value_df <- data.frame(PC = character(), P_Value_Base = numeric(), 
                         P_Value_Exponent = numeric(), stringsAsFactors = FALSE)

for (pc in pc_list) {
  # Perform t-test
  p_value <- t.test(as.formula(paste(pc, "~ Sex")), data = pca_data_annot, var.equal = TRUE, alternative = "two.sided")$p.value
  
  # Convert p-value to scientific notation
  p_value_sci <- formatC(p_value, format = "e", digits = 1)
  
  # Split into base and exponent
  p_value_split <- strsplit(p_value_sci, "e")[[1]]
  p_value_base <- as.numeric(p_value_split[1])
  p_value_exponent <- as.numeric(p_value_split[2])
  
  # Store results in dataframe
  p_value_df <- rbind(p_value_df, data.frame(PC = pc, P_Value_Base = p_value_base, P_Value_Exponent = p_value_exponent))
}



# Create PCA plot with side density plots
p <- ggplot(pca_data_annot, aes(x = Dim.1, y = Dim.2)) +
    geom_point(aes(fill = condition), size = 2, shape = 21, , 
               stroke = 0.1,  alpha = 0.7, show.legend = TRUE) +
    scale_fill_manual(values = c("royalblue", "darkorange1"), 
                      aesthetics = c("color", "fill"), labels = c("Female", "Male")) +
    geom_xsidedensity(aes(colour = condition, fill = condition), 
                      alpha = .3, linewidth = 0.7, show.legend = FALSE) +
    geom_ysidedensity(aes(colour = condition, fill = condition), 
                      alpha = .3, linewidth = 0.7, , show.legend = FALSE) +
    geom_xsidevline(xintercept = as.numeric(pca_data_median[1,"Dim.1"]), 
                    colour = "royalblue", linewidth = 0.7, linetype = "solid") +
    geom_xsidevline(xintercept = as.numeric(pca_data_median[2,"Dim.1"]), 
                    colour = "darkorange1", linewidth = 0.7, , linetype = "solid") +
    geom_ysidehline(yintercept = as.numeric(pca_data_median[1,"Dim.2"]), 
                    colour = "royalblue", linewidth = 0.7, , linetype = "solid") +
    geom_ysidehline(yintercept = as.numeric(pca_data_median[2,"Dim.2"]), 
                    colour = "darkorange1", linewidth = 0.7, , linetype = "solid") +
    xlab(paste0("PC1: ", "(", round(res.pca$eig[1,2], 2), "%)")) +
    ylab(paste0("PC2: ", "(", round(res.pca$eig[2,2], 2), "%)")) +
    theme(
          plot.background = element_rect(fill = NA, color = NA),
          plot.margin = margin(t = 2, r = 20, b = 2, l = 2),
          panel.border = element_rect(linetype = "solid", size = 0.5, colour = "black", fill = NA),
          ggside.panel.border = element_blank(),
          ggside.panel.grid = element_blank(),
          ggside.panel.background = element_blank(),
          ggside.axis.text = element_blank(),
          ggside.axis.ticks = element_blank(),
          ggside.panel.scale = .2,
          axis.line = element_line(size = 0.3, colour = "black"),
          axis.text = element_text(size = 7, colour = "black"),
          axis.title.y = element_text(size = 8, angle = 90, hjust = 0.38),
          axis.title.x = element_text(size = 8, angle = 0, hjust = 0.32),
          plot.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          legend.key.width = unit(0.6, "mm"),
          legend.key.spacing = unit(1, "mm"),
          legend.box.margin=margin(t = 0, r = 0, b = -18, l = 0),
          legend.background = element_rect(fill = NA)
          ) +
    guides(colour = "none") +
    # Set the position of the Pvalues of the density plots
    annotation_custom(grob = textGrob(label = bquote(~italic("P") ~ "value:"~ .(p_value_df$P_Value_Base[1]) ~ "x" ~ 
                                                       10^.(p_value_df$P_Value_Exponent[1])),
                                      hjust = 0, rot = 0, gp = gpar(cex = 0.6, fontface = "italic")),
                      xmin = -9, xmax = -9,
                      ymin = 12, ymax = 12) + # Set the position of the lable manually or via calculation of PC1
    annotation_custom(grob = textGrob(label = bquote(~italic("P") ~ "value:"~ .(p_value_df$P_Value_Base[2]) ~ "x" ~ 
                                                       10^.(p_value_df$P_Value_Exponent[2])), 
                                      hjust = 0, rot = 270, gp = gpar(cex = 0.6, fontface = "italic")),
                      xmin = 10, xmax = 10,
                      ymin = 10, ymax =10) # Set the position of the lable manually or via calculation of PC2
  
# Code to override clipping
gt <- ggplot_gtable(ggplot_build(p))

grob_names <- gt$layout$name

# Loop through each grob name and set clip to "off"
for (name in grob_names) {
  gt$layout$clip[gt$layout$name == name] <- "off"
}
grid::grid.draw(gt)


pdf("PCA_random_data.pdf", width = 2.6, height = 2.6,bg = "transparent")
grid::grid.draw(gt)
dev.off()

