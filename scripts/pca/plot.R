# Setup

library(ggplot2)
base_path = "/data/biol-silvereye/ball6625/norfolk-island/outputs/pca"

# Read covariance matrix and sample names
cov <- as.matrix(read.table(file.path(base_path,"output.cov")))
sample_names <- readLines("/data/biol-silvereye/ball6625/norfolk-island/samples.txt")

# Convert eigenvectors to df
e<-eigen(cov)
eigen_df <- as.data.frame(e$vectors[, 1:2])
colnames(eigen_df) <- c("PC1", "PC2")
eigen_df$Index <- 1:nrow(eigen_df)

# Add the sample labels
eigen_df$Label <- sample_names

# Plot
plot <- ggplot(eigen_df, aes(x = PC1, y = PC2, label = Label)) +
  geom_point(size=3) +
  geom_text(vjust = -0.2, hjust = -0.2) +  # Adjust label position
  #ylim(-1,1.3)+
  #xlim(-1,1.3)+
  theme_classic() +
  labs(x = "PC1",
       y = "PC2")

plot

ggsave(plot, width=4, height= 4, filename=file.path(base_path,"plot.png"))