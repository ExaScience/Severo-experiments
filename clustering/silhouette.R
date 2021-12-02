library(patchwork)
library(ggplot2)
library(reshape2)
library(cowplot)
library(dplyr)
library(tidyr)
library(rhdf5)

imec.cols <- c("#3F98BD", "#929497", "#90298D", "#36337D", "#1582BE", "#99BDE4", "#C778AD", "#52BDC2", "#3F98BD", "#2D6C85")
cols <- c("#ff7f0e","#922b21","#76448a","#2874a6","#148f77","#1d8348","#d68910","#ba4a00", "#3F98BD")
cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

palette <- c("red4", "darkslategray3", "dodgerblue1", "darkcyan",
               "gray79", "black", "skyblue2", "dodgerblue4",
               "purple4", "maroon", "chocolate1", "bisque3", "bisque",
               "seagreen4", "lightgreen", "skyblue4", "mediumpurple3",
               "palevioletred1", "lightsalmon4", "darkgoldenrod1")

cols <- c("seurat"="#F8766D", "scanpy"="#00BA38", "severo"="#619CFF")

load_data <- function(implementation, dataset, params) {
  n <- paste(implementation, dataset, params, sep="/")
  lbls <- h5read("results.h5", paste0(n, "/lbls"))
  width <- h5read("results.h5", paste0(n, "/width"))
  data.frame(cluster=as.factor(lbls), width=width)
}

plot_silhouette <- function(df) {
  df <- df[order(df$cluster, -df$width), ]
  avg_silh <- mean(df$width)
  ggplot(df, aes(x=1:length(width), y=width, color=cluster, fill=cluster)) +
  	geom_bar(stat = "identity") +
	labs(y = "Silhouette width", x = "", title = "Clusters silhouette plot") +
	annotate("text", label=paste0("Average silhouette width: ", round(avg_silh, 2)),  x=Inf, y = Inf, vjust=1, hjust=1) +
	ggplot2::ylim(c(NA, 1)) +
	geom_hline(yintercept=avg_silh, linetype="dashed", color="red") +
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

df <- load_data("seurat", "Human_M1_10xV3-D1", "10_20_20_0.1")
p <- plot_silhouette(df) + ggtitle("Silhouette width: optimal ASW")
df <- load_data("seurat", "Human_M1_10xV3-D1", "20_10_40_0.05")
q <- plot_silhouette(df) + ggtitle("Silhouette width: optimal ARI")
comb <- (p + q)
