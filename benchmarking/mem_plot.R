library(patchwork)
library(ggplot2)
library(reshape2)
library(cowplot)
library(dplyr)
library(tidyr)
library(lmtest)
library(lme4)

cols <- c("seurat"="#F8766D", "scanpy"="#00BA38", "severo"="#619CFF")

round_size <- function(size) {
    sizes <- c(3000, 6000, 12000, 25000, 50000, 75000, 100000, 125000, 250000, 500000, 750000, 1000000)
    i <- findInterval(size, sizes, left.open=T)
    sizes[i+1]
}

X <- read.csv("comparison_metrics.csv") %>%
    dplyr::mutate(implementation = factor(implementation, levels=c("R", "py", "jl"), labels=c("seurat", "scanpy", "severo")),
        dataset = factor(dataset,
                    levels=c("1M_gz", "210129_raw_BALPBMC", "l5_all", "pbmc_68k", "3k"),
                    labels=c("Brain cells (1.3M)", "Covid-19 (500k)", "Mouse Brain Atlas (120k)", "pbmc_68k"="PBMC (68k)", "3k"="PBMC (3k)"),
        ),
    ) %>% tidyr::drop_na(implementation, dataset)

Y <- X %>% dplyr::group_by(dataset, size, implementation) %>% dplyr::summarize(ari=median(ari, na.rm=T), jaccard=median(jaccard, na.rm=T), peakmem=median(peakmem, na.rm=T)) %>% dplyr::ungroup()

X1 <- read.csv("ari_nohvf.csv")
X2 <- read.csv("ari_hvf.csv")
XX <- rbind(cbind(hvf=F,X1), cbind(hvf=T,X2)) %>% dplyr::filter(implementation %in% c("R", "jl", "py")) %>%
    dplyr::mutate(implementation = factor(implementation, levels=c("R", "py", "jl"), labels=c("seurat", "scanpy", "severo")), size = round_size(size))
YY <- XX %>% dplyr::group_by(dataset, size, implementation) %>% dplyr::summarize(ari=median(ari, na.rm=T), jaccard=median(jaccard, na.rm=T), peakmem=median(peakmem, na.rm=T)) %>% dplyr::ungroup()

p <- ggplot(Y, aes(x=dataset, y=peakmem, group=implementation, fill=implementation, color=implementation)) +
    geom_col(position="dodge") +
    scale_color_manual(values=cols) +
    geom_hline(aes(yintercept=256), color="red") + annotate("text", x="PBMC (3k)", y=256, label="maximum system memory", color="red", vjust=-1) +
    ylab("Peak memory (GB)") + xlab(NULL)

Q <- rbind(Y, YY) %>% dplyr::filter(!is.na(peakmem))


model_comparison <- function(x) {
	lm1 <- lmList(peakmem ~ size - 1 | implementation, x)
	lm2 <- lmList(peakmem ~ I(size^2) + size - 1 | implementation, x)

	t <- mapply(lrtest, lm1, lm2, SIMPLIFY=F)
	z <- mapply(function(t, a, b) {
			pval <- t$`Pr(>Chisq)`[2]
			if(pval < 0.01) b else a
		}, t, lm1, lm2, SIMPLIFY=F)
	list(best=z, test=t)
}

predict_lmList <- function(mod, newdata, ...) {
	newdata <- as.data.frame(newdata)
	z <- do.call(rbind, mapply(function(x,i) predict(x, newdata[newdata$implementation == i,], ...), mod, names(mod)))
	merge(newdata, z, by="row.names", all.x=T)
}

o <- Q %>% model_comparison()
D <- predict_lmList(o$best, Q, interval = 'confidence')

q <- ggplot(D, aes(x=size, y=peakmem, group=implementation, fill=implementation, color=implementation)) +
    geom_line(aes(x=size, y=fit)) +
    geom_ribbon(aes(x=size, y=fit, ymin=lwr, ymax=upr), alpha=0.3) +
    geom_point() +
    geom_hline(aes(yintercept=256), color="red") + annotate("text", x=100000, y=256, label="maximum system memory", color="red", vjust=-1) +
    scale_color_manual(values=cols) +
    ylab("Peak memory (GB)") + xlab("Number of cells")

ggsave(file="memory_usage_bars.pdf", plot=p, width=20, height=15)
ggsave(file="memory_usage_scaling.pdf", plot=q, width=20, height=15)

comb <- (p | q) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom')
p_ranges_y <- c(ggplot_build(comb[[1]])$layout$panel_scales_y[[1]]$range$range,
                ggplot_build(comb[[2]])$layout$panel_scales_y[[1]]$range$range)
comb <- comb & ylim(min(p_ranges_y), max(p_ranges_y))

ggsave(file="memory_usage_comb.pdf", plot=comb, width=20, height=15)
