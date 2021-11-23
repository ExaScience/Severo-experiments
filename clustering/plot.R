library(patchwork)
library(ggplot2)
library(reshape2)
library(cowplot)
library(dplyr)
library(tidyr)

imec.cols <- c("#3F98BD", "#929497", "#90298D", "#36337D", "#1582BE", "#99BDE4", "#C778AD", "#52BDC2", "#3F98BD", "#2D6C85")
cols <- c("#ff7f0e","#922b21","#76448a","#2874a6","#148f77","#1d8348","#d68910","#ba4a00", "#3F98BD")
cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

palette <- c("red4", "darkslategray3", "dodgerblue1", "darkcyan",
               "gray79", "black", "skyblue2", "dodgerblue4",
               "purple4", "maroon", "chocolate1", "bisque3", "bisque",
               "seagreen4", "lightgreen", "skyblue4", "mediumpurple3",
               "palevioletred1", "lightsalmon4", "darkgoldenrod1")

cols <- c("seurat"="#F8766D", "scanpy"="#00BA38", "severo"="#619CFF")

X <- read.csv("comparison.csv", stringsAsFactors=T)
X <- X %>% tidyr::pivot_longer(c(ari, purity, ri, classpurity, clusters), names_to="score", values_to="value")

Z <- X %>% dplyr::group_by(dataset, implementation, score) %>% dplyr::summarize(value=max(value)) %>% dplyr::ungroup()

Y <- X %>% dplyr::group_by(dataset, implementation, resolution, dims, k, tables, score) %>% dplyr::summarize(mu=mean(value), sd=sd(value), n=length(value)) %>% dplyr::ungroup()
Y <- Y %>% dplyr::filter(score == "ari")
Y %>% dplyr::group_by(dataset, implementation, score) %>% dplyr::summarize(sd=mean(sd, na.rm=T))

p <- ggplot(Y, aes(x=resolution, y=mu, color=implementation)) +
     geom_pointrange(aes(ymin=mu-sd, ymax=mu+sd), position=position_dodge(width=0.1)) +
#    scale_fill_manual(values=cols) +
    facet_grid(score~dataset, scales="free_y") + ylab("Score") + xlab(NULL) + theme(panel.background=element_rect("#F8f8f8"))

ggsave(file="comparison_clusters.pdf", plot=p, width=20, height=15)
