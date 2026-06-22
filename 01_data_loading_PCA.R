# ===================================================
# DNA Methylation Aging Analysis - GSE40279
# Step 1: Data Loading, Metadata Merge, PCA Exploration
# ===================================================

library(data.table)
library(GEOquery)
library(ggplot2)

# ---- Load beta matrix ----
beta <- fread("/Users/nousheenjahanshaik/Downloads/GSE40279_average_beta.txt")

# set CpG IDs as rownames
rownames(beta) <- beta$ID_REF
beta$ID_REF <- NULL

dim(beta)

# ---- Load metadata from GEO ----
#gse <- getGEO("GSE40279", GSEMatrix = TRUE)
gse <- getGEO("GSE40279", GSEMatrix = TRUE, getGPL = FALSE)
metadata <- pData(phenoData(gse[[1]]))

# ---- Build clean pheno table ----
pheno <- data.frame(
  sample_id = metadata$source_name_ch1,
  age = as.numeric(metadata$`age (y):ch1`),
  gender = metadata$`gender:ch1`,
  ethnicity = metadata$`ethnicity:ch1`
)

# ---- Check alignment ----
all(colnames(beta) == pheno$sample_id)  # should be TRUE

# ---- PCA ----
pca <- prcomp(t(beta), scale. = FALSE)
summary(pca)$importance[, 1:5]

# ---- PCA dataframe for plotting ----
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  age = pheno$age,
  gender = pheno$gender,
  ethnicity = pheno$ethnicity
)

# ---- Age groups for visualization ----
pca_df$age_group <- cut(pheno$age,
                        breaks = c(0, 40, 60, 80, 101),
                        labels = c("20-40", "41-60", "61-80", "81-101"))

# ---- PCA plots ----
ggplot(pca_df, aes(x=PC1, y=PC2, color=age_group)) +
  geom_point(size=2) +
  scale_color_manual(values = c("20-40" = "blue",
                                "41-60" = "green",
                                "61-80" = "orange",
                                "81-101" = "red")) +
  labs(title="PCA colored by Age Group") +
  theme_bw()

ggplot(pca_df, aes(x=PC1, y=PC2, color=gender)) +
  geom_point(size=2) +
  labs(title="PCA colored by Gender") +
  theme_bw()

ggplot(pca_df, aes(x=PC1, y=PC2, color=ethnicity)) +
  geom_point(size=2) +
  labs(title="PCA colored by Ethnicity") +
  theme_bw()

# ---- Identify outliers ----
outliers <- rownames(pca_df[pca_df$PC1 > 30, ])
outliers

# ---- Sample size balance ----
table(pheno$gender)
table(pheno$ethnicity)
table(pheno$gender, pheno$ethnicity)

# ---- Confounding check: age vs ethnicity/gender ----
boxplot(pheno$age ~ pheno$ethnicity)
boxplot(pheno$age ~ pheno$gender)
