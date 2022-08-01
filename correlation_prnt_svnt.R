library(tidyverse)
library(boot)
library(patchwork)
library(blandr)

df_meta <- read_csv("../data/sample_groups.csv")
df_meta <- df_meta %>% filter(!grepl("Don’t work on these", note))
df_meta <- df_meta %>% filter(!grepl("D700 Rec", group))
df_meta$parent_group2 <- gsub(" +\\+", "\\+", df_meta$parent_group2)
df_meta$parent_group2 <- gsub("\\+ +", "+", df_meta$parent_group2)
df_meta$parent_group2 <- gsub("\\+", " \\+ ", df_meta$parent_group2)
colors_t <- c("#000000", "#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00", "#a65628","#f781bf","#999999") # colorbrewer2 color for groups
names(colors_t) <- unique(df_meta$parent_group2)

df_correlation <- read_csv("../data/correlation_data.csv")

response_all <- df_correlation %>% select(WT_PRNT_50:`BA1_sVNT`) %>% colnames()
RBD_all <- sapply(response_all, function(x){strsplit(x, "_")[[1]][1]})

df_correlation <- left_join(df_correlation, df_meta %>% select(-sample), "group")
df_correlation <- df_correlation %>% filter(!is.na(parent_group2))
df_correlation$color <- sapply(df_correlation$parent_group2, function(x){colors_t[names(colors_t)==x]})

unique(df_correlation$group)
df_correlation$open_circle <- FALSE
df_correlation$open_circle[grepl("Pre", df_correlation$group)] <- TRUE
df_correlation$open_circle[grepl("^D0", df_correlation$group)] <- TRUE
df_correlation$open_circle[grepl("Acute", df_correlation$group)] <- TRUE
df_correlation$open_circle[grepl("Neg", df_correlation$group)] <- TRUE


p1 <- ggplot(df_correlation) +
	geom_boxplot(aes(x=group, y=BA1_beads, color=parent_group2))+
	scale_color_manual(name = "Group", values=colors_t)+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
	guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
	scale_x_discrete(guide = guide_axis(n.dodge = 1), expand=c(0.1,1.5)) +
	xlab("")+
	ylab("BA1 beads")+
	NULL
p2 <- ggplot(df_correlation) +
	geom_boxplot(aes(x=group, y=BA1_PRNT_50, color=parent_group2))+
	scale_color_manual(name = "Group", values=colors_t)+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
	guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
	scale_x_discrete(guide = guide_axis(n.dodge = 1), expand=c(0.1,1.5)) +
	xlab("")+
	ylab("BA1 PRNT 50")+
	NULL
p3 <- ggplot(df_correlation) +
	geom_boxplot(aes(x=group, y=BA1_sVNT, color=parent_group2))+
	scale_color_manual(name = "Group", values=colors_t)+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
	guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
	scale_x_discrete(guide = guide_axis(n.dodge = 1), expand=c(0.1,1.5)) +
	xlab("")+
	ylab("BA1 sVNT")+
	NULL

p_out <- ((p1+ggtitle("A"))/(p2+ggtitle("B"))/(p3+ggtitle("C"))) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("../results/BA1_responses.pdf", width = 10, height = 12, plot=p_out)

ggplot(df_correlation) +
	geom_point(aes(x=BA1_beads, y=WT_beads, color=parent_group2), alpha=0.8)+
	scale_color_manual(name = "Group", values=colors_t)+
	theme_bw()+
	theme(legend.position = "top")+
	guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
	xlab("BA1 beads response")+
	ylab("WT beads response")+
	NULL
ggsave("../results/WT_BA1_beads_response_2d.pdf", width = 8, height = 8)


# Boot
func_boot_cor <- function(data, indices, var1, var2){
	return(cor(data[indices,][[var1]], data[indices,][[var2]], method="spearman"))
}

plots <- list()
plots_meta <- c()
i <- 1
types <- c("PRNT_50", "sVNT", "beads")
names(types) <- c("Viral PRNT", "Plate sVNT", "Bead sVNT")
lapply(unique(RBD_all), function(this_RBD){
	print(this_RBD)
	# this_RBD <- "BA2"
	df_tmp <- df_correlation %>% select(contains(this_RBD)) 
	df_tmp_bind <- df_tmp %>% bind_cols(df_correlation %>% select(group, sample, parent_group2, open_circle))

	experiments_all_i <- gsub(paste0(this_RBD, "_"), "", colnames(df_tmp))
	if(length(experiments_all_i)<2){return(NA)}
	experiments_pairs_i <- combn(experiments_all_i, 2)
	plots_i <- apply(experiments_pairs_i, 2, function(pair){
		# pair <- c("PRNT_50", "beads")
		# pair <- c("sVNT", "beads")
		print(pair)
		
		var1 <- paste0(this_RBD, "_", pair[1])
		var2 <- paste0(this_RBD, "_", pair[2])

		label1 <- names(types)[types==pair[1]]
		label2 <- names(types)[types==pair[2]]
		
		check_row <- df_tmp_bind %>% select(all_of(var1), all_of(var2)) %>% apply(1,function(x){any(is.na(x))})
		df_tmp_bind_filter <- df_tmp_bind[!check_row,]

		set.seed(2022)
		rst_boot <- boot(df_tmp_bind_filter, func_boot_cor, R=10000, var1=var1, var2=var2)
		# print(rst_boot)
		rst_boot_ci <- boot.ci(rst_boot, type="perc")
		text_i <- paste0("r = ", round(mean(rst_boot$t),2), "\n", "95% CI = ", round(rst_boot_ci$percent[4],2), "~", round(rst_boot_ci$percent[5],2))

		## plot for correlation
		colors_t_i <- colors_t[names(colors_t) %in% unique(df_tmp_bind_filter$parent_group2)]
		p_i <- ggplot(df_tmp_bind_filter, aes_string(x=var1, y=var2))+
			geom_point(aes_string(color="parent_group2"), alpha=0.8, shape=1, size=2, data=. %>% filter(open_circle))+
			geom_point(aes_string(color="parent_group2", fill="parent_group2"), alpha=0.8, shape=16, size=2, data=. %>% filter(!open_circle))+
			geom_smooth(method = 'loess', span = 10, color="black")+
			# geom_smooth(method = 'lm', color="black")+
			theme_bw()+
			# scale_color_manual(name = "Group", values=colors_t_i)+
			# scale_fill_manual(name = "Group", values=colors_t_i)+
			scale_color_manual(name = "Group", values=colors_t)+
			scale_fill_manual(name = "Group", values=colors_t)+
			# theme(legend.position = "bottom")+
			xlab(label1)+
			ylab(label2)

		x_range <- layer_scales(p_i)$x$range$range
		x_pos <- x_range[2]-(x_range[2]-x_range[1])*0.2
		x_pos2 <- x_range[1]+(x_range[2]-x_range[1])*0.05
		y_range <- layer_scales(p_i)$y$range$range
		y_pos <- y_range[1]+(y_range[2]-y_range[1])*0.1
		y_pos2 <- y_range[2]-(y_range[2]-y_range[1])*0.05
		
		if(grepl("PRNT", pair[1])){
			p_i <- p_i + scale_x_continuous(breaks = c(20, 40, 80, 160, 320)) +
				annotate("text", x=x_pos, y=y_pos, label=text_i)+
				annotate("text", x=x_pos2, y=y_pos2, label=this_RBD)
		} else {
			p_i <- p_i + 
				annotate("text", x=x_pos, y=y_pos, label=text_i)+
				annotate("text", x=x_pos2, y=y_pos2, label=this_RBD)
		}

		p_i_out <- p_i + theme(legend.position = "none")
		ggsave(paste0("../results/cor_",this_RBD, "_", label1, "_", label2, ".pdf"), plot=p_i_out, width=6, height=6)

		## plot for Bland Altman analyses 
		if(!grepl("PRNT", var1)){
			response1 <- df_tmp_bind_filter[[var1]]
			response2 <- df_tmp_bind_filter[[var2]]
			p_bla <- blandr.draw(response1, response2)
			p_bla <- p_bla+ylab(paste0(label1, " - ", label2))
			ggsave(paste0("../results/Bland_Altman_",this_RBD, "_", label1, "_", label2, ".pdf"), plot=p_bla, width=6, height=6)
			blandr.output.text(response1, response2)
		}
		

		p_i <- p_i + ggtitle(LETTERS[i])
		plots <<- c(plots, list(p_i))
		plots_meta <<- c(plots_meta, paste0(this_RBD, " ", label1, " ", label2))
		i <<- i + 1
	})
	return("")
})

plots_meta
p_out <- (plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + guide_area() + plots[[7]] + plots[[8]] + plot_layout(guides = 'collect'))
# p_out <- (plots[[1]]/plots[[2]]/plots[[3]]/plots[[4]])|(plots[[5]]/plots[[6]]/plots[[7]]/plots[[8]])
ggsave("../results/correlation_experiments.pdf", width = 15, height = 15, plot=p_out)
