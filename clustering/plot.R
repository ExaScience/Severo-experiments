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
X <- X %>% tidyr::pivot_longer(c(ari, purity, ri, classpurity, clusters, time, avgsilh), names_to="score", values_to="value") %>%
	dplyr::group_by(dataset, implementation, resolution, dims, k, tables, score) %>%
	dplyr::summarize(med=median(value), mu=mean(value), sd=sd(value), ymin=min(value), ymax=max(value), n=length(value)) %>%
	dplyr::ungroup()

Z <- X %>% tidyr::pivot_wider(names_from=score, values_from=c(med, mu, sd, ymin, ymax, n)) %>%
	dplyr::group_by(dataset, implementation) %>% dplyr::top_n(n=1, wt=med_ari) %>% dplyr::ungroup() %>%
   	tidyr::pivot_longer(contains("_"), names_pattern="(.*)_(.*)", names_to=c(".value", "score"))

Y <- X %>% tidyr::pivot_wider(names_from=score, values_from=c(med, mu, sd, ymin, ymax, n)) %>%
	dplyr::group_by(dataset, implementation) %>% dplyr::top_n(n=1, wt=med_avgsilh) %>% dplyr::ungroup() %>%
   	tidyr::pivot_longer(contains("_"), names_pattern="(.*)_(.*)", names_to=c(".value", "score"))

p <- ggplot(Z, aes(x=dataset, fill=implementation, group=implementation, y=med)) +
   	geom_bar(stat="identity", position="dodge") +
	geom_errorbar(aes(ymin=mu-sd, ymax=mu+sd), width=.2, position=position_dodge(.9)) +
	facet_wrap(~score, scales="free_y")

ggsave(file="comparison_optimal.pdf", plot=p, width=20, height=20)

X %>% dplyr::group_by(dataset, implementation, score) %>% dplyr::summarize(sd=mean(sd, na.rm=T))

p <- ggplot(X %>% dplyr::filter(dataset == "Human_M1_10xV3-D2", score == "ari"), aes(x=factor(resolution), y=mu, color=implementation)) +
     geom_boxplot(aes(middle=med, lower=mu-sd, upper=mu+sd, ymin=ymin, ymax=ymax), stat="identity") +
#    scale_fill_manual(values=cols) +
    facet_grid(dims~k, labeller = label_both) + ylab("ARI") + xlab("resolution") + theme(panel.background=element_rect("#F8f8f8"))

ggsave(file="comparison_clusters.pdf", plot=p, width=20, height=20)

ggplot(X %>% dplyr::filter(dataset == "Human_M1_10xV3-D1", score == "avgsilh" | score == "ari"), aes(x=resolution, y=med, color=implementation)) +
	geom_pointrange(aes(ymin=mu-sd, ymax=mu+sd)) +
   	facet_grid(dims~score,labeller = label_both) +
   	ylab("avg silhouette width")

res <- h5read("umap.h5", "Human_M1_10xV3-D1/10")
lbls_true <- h5read("umap.h5", "Human_M1_10xV3-D1/lbls")
lbls_R <- h5read("results.h5", "seurat/Human_M1_10xV3-D1/10_20_20_0.1/lbls")
lbls_py <- h5read("results.h5", "scanpy/Human_M1_10xV3-D1/10_20_20_0.1/lbls")
lbls_jl <- h5read("results.h5", "severo/Human_M1_10xV3-D1/10_20_20_0.1/lbls")
Q <- data.frame(x=res[,1], y=res[,2], lbls_R=factor(lbls_R), lbls_true=factor(lbls_true), lbls_py=factor(lbls_py), lbls_jl=factor(lbls_jl))
wrap_plots(
	ggplot(Q, aes(x=x, y=y)) + geom_point(aes(color=lbls_R)),
	ggplot(Q, aes(x=x, y=y)) + geom_point(aes(color=lbls_true)),
	ggplot(Q, aes(x=x, y=y)) + geom_point(aes(color=lbls_py)),
	ggplot(Q, aes(x=x, y=y)) + geom_point(aes(color=lbls_jl)))

