# Quantify antigenic diversity on different RBDs for different groups,
# and comparison between groups.

library(readxl)
library(tidyverse)
library(boot)
library(parallel)
library(scico)
library(ggrepel)
library(patchwork)

# read data
df_meta <- read_csv("../data/sample_groups.csv")
df_meta <- df_meta %>% filter(!grepl("Don’t work on these", note))
df_meta <- df_meta %>% filter(!grepl("D700 Rec", group))
df_meta$parent_group2 <- gsub(" +\\+", "\\+", df_meta$parent_group2)
df_meta$parent_group2 <- gsub("\\+ +", "+", df_meta$parent_group2)
df_meta$parent_group2 <- gsub("\\+", " \\+ ", df_meta$parent_group2)

df_inhibition <- read_excel("../data/data_inhibition.xlsx")
df_inhibition <- df_inhibition %>% pivot_longer(-c(group, sample), names_to = "RBD")

df_inhibition <- df_inhibition %>% filter(group %in% df_meta$group)# filter out the groups which were not needed
df_inhibition$Response <- ifelse(df_inhibition$value<20, "negative", "positive")

# double check
RBDs_all <- unique(df_inhibition$RBD)
groups_all <- unique(df_meta$group)
parent_groups_all <- unique(df_meta$parent_group2)
stopifnot(all(df_inhibition$group %in% groups_all))

# first let's test the negative controls:
## whether the distribution of negative controls are similar between different settings. (inhibition)
df_tmp <- df_inhibition %>% filter(group%in%c("Neg", "Pre-BB", "Pre-CC"))
ggplot(df_tmp, aes(x= RBD, y=value)) +
	geom_boxplot()+
	geom_point(aes(color=Response), size=1, alpha=0.6)+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1))+
	scale_x_discrete(guide = guide_axis(n.dodge = 1), expand=c(0,1.5)) +
	facet_wrap(vars(group), ncol=1)+
	ylab("% Inhibition")+
	ggtitle("Neg")+
	NULL
ggsave("../results/negative_samples_3.pdf", width = 8, height=12)
kruskal.test(df_tmp$value, df_tmp$RBD) # p = 1, The differences between the medians are not statistically significant

df_tmp <- df_inhibition %>% filter(group=="Neg")
ggplot(df_tmp, aes(x= RBD, y=value)) +
	geom_boxplot()+
	geom_point(aes(color=Response), size=1, alpha=0.6)+
	theme_minimal()+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1))+
	scale_x_discrete(guide = guide_axis(n.dodge = 1), expand=c(0,1.5)) +
	ylab("% Inhibition")+
	ggtitle("Neg")+
	NULL
ggsave("../results/negative_samples.pdf", width = 8, height=6)

df_tmp$RBD2 <- df_tmp$RBD
df_tmp$RBD2[df_tmp$Response=="negative"] <- "Negatives"
df_tmp$RBD2 <- factor(df_tmp$RBD2, levels=c(RBDs_all, "Negatives"))
mean(df_tmp$value[df_tmp$Response=="positive"])
df_tmp <- df_tmp %>% group_by(RBD2) %>% summarise(N=n())
df_tmp_freq <- df_tmp
N_total <- sum(df_tmp_freq$N)
df_tmp_freq$p <- df_tmp_freq$N/N_total
(hs <- sum(-df_tmp_freq$p*log(df_tmp_freq$p)))
(hpi <- 1-(sum(df_tmp_freq$N*(df_tmp_freq$N-1))/(N_total*(N_total-1))))

ggplot(df_tmp) +
	geom_col(aes(x=RBD2,y=N))+
	geom_text(aes(x=RBD2,y=N+10, label=N))+
	scale_x_discrete(guide = guide_axis(n.dodge = 1), expand=c(0,1.5), drop=F) +
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1))+
	xlab("Responses")+
	ylab("Number of samples")+
	ggtitle("Neg")+
	NULL
ggsave("../results/negative_samples_responses.pdf", width = 8, height=6)


df_tmp <- df_inhibition %>% filter(group=="BB")
ggplot(df_tmp, aes(x= RBD, y=value)) +
	geom_boxplot()+
	geom_point(aes(color=Response), size=1, alpha=0.6)+
	theme_minimal()+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1))+
	scale_x_discrete(guide = guide_axis(n.dodge = 1), expand=c(0,1.5)) +
	ylab("% Inhibition")+
	ggtitle("BB")+
	NULL
ggsave("../results/BB_samples.pdf", width = 8, height=6)

df_tmp$RBD2 <- df_tmp$RBD
df_tmp$RBD2[df_tmp$Response=="negative"] <- "Negatives"
df_tmp$RBD2 <- factor(df_tmp$RBD2, levels=c(RBDs_all, "Negatives"))
df_tmp <- df_tmp %>% group_by(RBD2) %>% summarise(N=n())

ggplot(df_tmp) +
	geom_col(aes(x=RBD2,y=N))+
	geom_text(aes(x=RBD2,y=N+10, label=N))+
	scale_x_discrete(expand=c(0,1.5), drop=F) +
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1))+
	xlab("Responses")+
	ylab("Number of samples")+
	ggtitle("BB")+
	NULL
ggsave("../results/BB_samples_responses.pdf", width = 8, height=6)


df_tmp <- df_inhibition %>% filter(group=="D30-BBB")
ggplot(df_tmp, aes(x= RBD, y=value)) +
	geom_boxplot()+
	geom_point(aes(color=Response), size=1, alpha=0.6)+
	theme_minimal()+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1))+
	scale_x_discrete(guide = guide_axis(n.dodge = 1), expand=c(0,1.5)) +
	ylab("% Inhibition")+
	ggtitle("D30-BBB")+
	NULL
ggsave("../results/D30-BBB_samples.pdf", width = 8, height=6)

df_tmp$RBD2 <- df_tmp$RBD
df_tmp$RBD2[df_tmp$Response=="negative"] <- "Negatives"
df_tmp$RBD2 <- factor(df_tmp$RBD2, levels=c(RBDs_all, "Negatives"))
mean(df_tmp$value[df_tmp$Response=="positive"])
df_tmp <- df_tmp %>% group_by(RBD2) %>% summarise(N=n())
df_tmp_freq <- df_tmp
N_total <- sum(df_tmp_freq$N)
df_tmp_freq$p <- df_tmp_freq$N/N_total
(hs <- sum(-df_tmp_freq$p*log(df_tmp_freq$p)))
(hpi <- 1-(sum(df_tmp_freq$N*(df_tmp_freq$N-1))/(N_total*(N_total-1))))

ggplot(df_tmp) +
	geom_col(aes(x=RBD2,y=N))+
	geom_text(aes(x=RBD2,y=N+3, label=N))+
	scale_x_discrete(expand=c(0,1.5), drop=F) +
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1))+
	xlab("Responses")+
	ylab("Number of samples")+
	ggtitle("D30-BBB")+
	NULL
ggsave("../results/D30-BBB_samples_responses.pdf", width = 8, height=6)


# When analyzing antigenic diversity, we have two questions to be answered:
# 1. Which groups have higher magnitude responses?
# 2. Which groups have boarder antibody responses?

# before starting the analysis, we have to introduce the definition of
# "positive response", which many of the following analyses are based on.
# Only responses (% inhibition) >= 20 are treated as "positive responses",
# otherwise the responses are treated as "negative responses".
# There are 16 different types (on 16 different RBDs) of "positive responses"
# in this study, but all negative responses (even on different RBDs)
# are classified in one "negative responses" group.


## 1. Which groups have higher magnitude responses?
## 2. Which groups have boarder antibody response?
func_boot <- function(data, indices){
	# df_tmp <- df_inhibition_i
	df_tmp <- data[indices,]
	### 1. Which groups have higher magnitude responses?
	#### we estimated the mean and 95CI of %inhibition of the "all responses"
	mean_res <- mean(df_tmp$value)
	#### we estimated the mean and 95CI of %inhibition of the "positive responses"
	mean_pos_res <- mean(df_tmp$value[df_tmp$Response=="positive"])
	
	### 2. Which groups have boarder antibody response?
	### https://academic.oup.com/ve/article/5/1/vey041/5304643
	df_tmp$RBD2 <- df_tmp$RBD
	df_tmp$RBD2[df_tmp$Response=="negative"] <- "Negatives"
	df_tmp_freq <- df_tmp %>% group_by(RBD2) %>% summarise(N=n())
	N_total <- nrow(df_tmp)
	df_tmp_freq$p <- df_tmp_freq$N/N_total

	#### we use shannon entropy to estimate the diversity
	hs <- sum(-df_tmp_freq$p*log(df_tmp_freq$p))

	#### we also use "antigenic diveristy (pi)" which is simialr to nucleotide diversity pi, to estimate the diversity
	hpi <- 1-(sum(df_tmp_freq$N*(df_tmp_freq$N-1))/(N_total*(N_total-1)))

	return(c(mean_res, mean_pos_res, hs, hpi))
}

### full spectrum/panel, including Sarbecovirus
num_repeat <- 10000
num_cpus <- 8
set.seed(2022)
list_boot_rst_full <- mclapply(groups_all, function(group_i){
	# print(group_i)
	# group_i <- "Neg"
	df_inhibition_i <- df_inhibition %>% filter(group == group_i)
	rst_boot <- boot(df_inhibition_i, func_boot, R=num_repeat, parallel = "no")
	ci_mean_res <- boot.ci(rst_boot, type=c("perc"), index = 1)
	ci_mean_pos_res <- boot.ci(rst_boot, type=c("perc"), index = 2)
	ci_hs <- boot.ci(rst_boot, type=c("perc"), index = 3)
	ci_hpi <- boot.ci(rst_boot, type=c("perc"), index = 4)
	# str(rst_boot)
	# str(ci_mean_pos_res)
	rst <- c(group_i, apply(rst_boot$t, 2, mean, na.rm=T), apply(rst_boot$t, 2, sd, na.rm=T), ci_mean_res$percent[4:5], ci_mean_pos_res$percent[4:5], ci_hs$percent[4:5], ci_hpi$percent[4:5])
	print(rst)
	rst
}, mc.cores = num_cpus)

df_boot_rst <- do.call(rbind, list_boot_rst_full)
colnames(df_boot_rst) <- c("group", "e_mean_res", "e_mean_pos_res", "e_hs", "e_hpi", "sd_mean_res", "sd_mean_pos_res", "sd_hs", "sd_hpi", "ci_mean_res_low", "ci_mean_res_high", "ci_mean_pos_res_low", "ci_mean_pos_res_high", "ci_hs_low", "ci_hs_high", "ci_hpi_low", "ci_hpi_high")
df_boot_rst <- as_tibble(df_boot_rst)
write_tsv(df_boot_rst, "../results/df_boot_rst_full.tsv")

### partial spectrum/panel, including only SARS-CoV-2 ancestral/VoC panel
num_repeat <- 10000
num_cpus <- 8
set.seed(2022)
RBDs_all
df_inhibition_partial <- df_inhibition %>% filter(grepl("SARS-CoV-2", RBD))

list_boot_rst_scov2 <- mclapply(groups_all, function(group_i){
	# print(group_i)
	# group_i <- "Neg"
	df_inhibition_i <- df_inhibition_partial %>% filter(group == group_i)
	rst_boot <- boot(df_inhibition_i, func_boot, R=num_repeat, parallel = "no")
	ci_mean_res <- boot.ci(rst_boot, type=c("perc"), index = 1)
	ci_mean_pos_res <- boot.ci(rst_boot, type=c("perc"), index = 2)
	ci_hs <- boot.ci(rst_boot, type=c("perc"), index = 3)
	ci_hpi <- boot.ci(rst_boot, type=c("perc"), index = 4)
	# str(rst_boot)
	# str(ci_mean_pos_res)
	rst <- c(group_i, apply(rst_boot$t, 2, mean, na.rm=T), apply(rst_boot$t, 2, sd, na.rm=T), ci_mean_res$percent[4:5], ci_mean_pos_res$percent[4:5], ci_hs$percent[4:5], ci_hpi$percent[4:5])
	print(rst)
	rst
}, mc.cores = num_cpus)

df_boot_rst <- do.call(rbind, list_boot_rst_scov2)
colnames(df_boot_rst) <- c("group", "e_mean_res", "e_mean_pos_res", "e_hs", "e_hpi", "sd_mean_res", "sd_mean_pos_res", "sd_hs", "sd_hpi", "ci_mean_res_low", "ci_mean_res_high", "ci_mean_pos_res_low", "ci_mean_pos_res_high", "ci_hs_low", "ci_hs_high", "ci_hpi_low", "ci_hpi_high")
df_boot_rst <- as_tibble(df_boot_rst)
write_tsv(df_boot_rst, "../results/df_boot_rst_scov2.tsv")

### partial spectrum/panel, including only Sarbecovirus panel
num_repeat <- 10000
num_cpus <- 8
set.seed(2022)
RBDs_all
df_inhibition_partial <- df_inhibition %>% filter(!grepl("SARS-CoV-2", RBD))

list_boot_rst_sarbecov <- mclapply(groups_all, function(group_i){
	# print(group_i)
	# group_i <- "Neg"
	df_inhibition_i <- df_inhibition_partial %>% filter(group == group_i)
	rst_boot <- boot(df_inhibition_i, func_boot, R=num_repeat, parallel = "no")
	ci_mean_res <- boot.ci(rst_boot, type=c("perc"), index = 1)
	ci_mean_pos_res <- boot.ci(rst_boot, type=c("perc"), index = 2)
	ci_hs <- boot.ci(rst_boot, type=c("perc"), index = 3)
	ci_hpi <- boot.ci(rst_boot, type=c("perc"), index = 4)
	# str(rst_boot)
	# str(ci_mean_pos_res)
	rst <- c(group_i, apply(rst_boot$t, 2, mean, na.rm=T), apply(rst_boot$t, 2, sd, na.rm=T), ci_mean_res$percent[4:5], ci_mean_pos_res$percent[4:5], ci_hs$percent[4:5], ci_hpi$percent[4:5])
	print(rst)
	rst
}, mc.cores = num_cpus)

df_boot_rst <- do.call(rbind, list_boot_rst_sarbecov)
colnames(df_boot_rst) <- c("group", "e_mean_res", "e_mean_pos_res", "e_hs", "e_hpi", "sd_mean_res", "sd_mean_pos_res", "sd_hs", "sd_hpi", "ci_mean_res_low", "ci_mean_res_high", "ci_mean_pos_res_low", "ci_mean_pos_res_high", "ci_hs_low", "ci_hs_high", "ci_hpi_low", "ci_hpi_high")
df_boot_rst <- as_tibble(df_boot_rst)
write_tsv(df_boot_rst, "../results/df_boot_rst_sarbecov.tsv")

# Draw
### full spectrum/panel, including Sarbecovirus
df_boot_rst <- read_tsv("../results/df_boot_rst_full.tsv")
df_boot_rst <- df_boot_rst %>% mutate_at(vars(-group), as.numeric)
df_boot_rst <- left_join(df_boot_rst, df_meta)
df_boot_rst$parent_group2 <- factor(df_boot_rst$parent_group2, levels=parent_groups_all)
df_boot_rst <- df_boot_rst %>% arrange(parent_group2)
df_boot_rst$group <- factor(df_boot_rst$group, levels=unique(df_boot_rst$group))

colors_t <- scico(length(unique(df_boot_rst$parent_group2)), palette = 'batlow')

p0 <- ggplot(df_boot_rst, aes(x=group)) +
	# geom_point(aes())+
	geom_pointrange(aes(y=e_mean_res, ymin=ci_mean_res_low, ymax=ci_mean_res_high, fill=parent_group2), size=0.8, shape=21, stroke=0.6)+
	# geom_point(aes(y=e_mean_pos_res), size=2, shape=21)+
	scale_x_discrete(expand=c(0,1.5), drop=F) +
	ylab("Average %inhibition of responses")+
	# scale_color_manual(name="Group", values=colors_t)+
	scale_fill_manual(name="Group", values=colors_t)+
	# facet_wrap(vars(parent_group2), nrow=3, scales="free_x")+
	theme_bw()+
	xlab("")+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
	guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
	NULL
ggsave("../results/Mean_inhibition_responses_full.pdf", width = 8, height=6)

p1 <- ggplot(df_boot_rst, aes(x=group)) +
	# geom_point(aes())+
	geom_pointrange(aes(y=e_mean_pos_res, ymin=ci_mean_pos_res_low, ymax=ci_mean_pos_res_high, fill=parent_group2), size=0.8, shape=21, stroke=0.6)+
	# geom_point(aes(y=e_mean_pos_res), size=2, shape=21)+
	scale_x_discrete(expand=c(0,1.5), drop=F) +
	ylab("Average %inhibition of positive responses")+
	# scale_color_manual(name="Group", values=colors_t)+
	scale_fill_manual(name="Group", values=colors_t)+
	# facet_wrap(vars(parent_group2), nrow=3, scales="free_x")+
	theme_bw()+
	xlab("")+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
	guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
	NULL
ggsave("../results/Mean_inhibition_positive_responses_full.pdf", width = 8, height=6)

p2 <- ggplot(df_boot_rst, aes(x=group)) +
	# geom_point(aes())+
	geom_pointrange(aes(y=e_hs, ymin=ci_hs_low, ymax=ci_hs_high, fill=parent_group2), size=0.8, shape=21, stroke=0.6)+
	scale_x_discrete(expand=c(0,1.5), drop=F) +
	ylab("Shannon entropy")+
	scale_fill_manual(name="Group", values=colors_t)+
	theme_bw()+
	xlab("")+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
	guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
	NULL
ggsave("../results/hs_full.pdf", width = 8, height=6)

p3 <- ggplot(df_boot_rst, aes(x=group)) +
	# geom_point(aes())+
	geom_pointrange(aes(y=e_hpi, ymin=ci_hpi_low, ymax=ci_hpi_high, fill=parent_group2), size=0.8, shape=21, stroke=0.6)+
	scale_x_discrete(expand=c(0,1.5), drop=F) +
	ylab(expression("RBD cross-reactivity ("~pi~")"))+
	scale_fill_manual(name="Group", values=colors_t)+
	theme_bw()+
	xlab("")+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
	guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
	NULL
ggsave("../results/hpi_full.pdf", width = 8, height=6)

p_out <- (p0+ggtitle("A"))/(p1+ggtitle("B")) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("../results/Antigenic_response_compare_full.pdf", width = 8, height=10)
p_out <- (p0+ggtitle("A"))/(p3+ggtitle("B")) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("../results/Antigenic_response_combined_full.pdf", width = 8, height=10)

p_out <- (p1+ggtitle("A"))/(p3+ggtitle("B"))/(p2+ggtitle("C")) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("../results/Antigenic_response_combined_add_hs_full.pdf", width = 8, height=12)



### partial spectrum/panel, including only SARS-CoV-2 ancestral/VoC panel
### partial spectrum/panel, including only Sarbecovirus panel
df_boot_rst0 <- read_tsv("../results/df_boot_rst_scov2.tsv")
df_boot_rst0$RBDs <- "SARS-CoV-2"
df_boot_rst1 <- read_tsv("../results/df_boot_rst_sarbecov.tsv")
df_boot_rst1$RBDs <- "Sarbecoviruses"

df_boot_rst <- bind_rows(df_boot_rst0, df_boot_rst1)
df_boot_rst <- df_boot_rst %>% mutate_at(vars(-group, -RBDs), as.numeric)
df_boot_rst <- left_join(df_boot_rst, df_meta)
df_boot_rst$parent_group2 <- factor(df_boot_rst$parent_group2, levels=parent_groups_all)
df_boot_rst <- df_boot_rst %>% arrange(parent_group2)
df_boot_rst$group <- factor(df_boot_rst$group, levels=unique(df_boot_rst$group))
df_boot_rst$group_rename <- factor(df_boot_rst$group_rename, levels=unique(df_boot_rst$group_rename))

colors_t <- scico(length(unique(df_boot_rst$parent_group2)), palette = 'batlow')

p0 <- ggplot(df_boot_rst, aes(x=group)) +
	# geom_point(aes())+
	geom_pointrange(aes(y=e_mean_res, ymin=ci_mean_res_low, ymax=ci_mean_res_high, fill=parent_group2), size=0.8, shape=21, stroke=0.6)+
	# geom_point(aes(y=e_mean_pos_res), size=2, shape=21)+
	scale_x_discrete(expand=c(0,1.5), drop=F) +
	ylab("Average %inhibition of responses")+
	# scale_color_manual(name="Group", values=colors_t)+
	scale_fill_manual(name="Group", values=colors_t)+
	facet_wrap(vars(RBDs), ncol=2, scales="free_x")+
	theme_bw()+
	xlab("")+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
	guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
	NULL
ggsave("../results/Mean_inhibition_responses_seperate.pdf", width = 12, height=6)

p1 <- ggplot(df_boot_rst, aes(x=group)) +
	# geom_point(aes())+
	geom_pointrange(aes(y=e_mean_pos_res, ymin=ci_mean_pos_res_low, ymax=ci_mean_pos_res_high, fill=parent_group2), size=0.8, shape=21, stroke=0.6)+
	# geom_point(aes(y=e_mean_pos_res), size=2, shape=21)+
	scale_x_discrete(expand=c(0,1.5), drop=F) +
	ylab("Average %inhibition of positive responses")+
	# scale_color_manual(name="Group", values=colors_t)+
	scale_fill_manual(name="Group", values=colors_t)+
	facet_wrap(vars(RBDs), ncol=2, scales="free_x")+
	theme_bw()+
	xlab("")+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
	guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
	NULL
ggsave("../results/Mean_inhibition_positive_responses_seperate.pdf", width = 12, height=6)

p2 <- ggplot(df_boot_rst, aes(x=group)) +
	# geom_point(aes())+
	geom_pointrange(aes(y=e_hs, ymin=ci_hs_low, ymax=ci_hs_high, fill=parent_group2), size=0.8, shape=21, stroke=0.6)+
	scale_x_discrete(expand=c(0,1.5), drop=F) +
	ylab("Shannon entropy")+
	scale_fill_manual(name="Group", values=colors_t)+
	facet_wrap(vars(RBDs), ncol=2, scales="free_x")+
	theme_bw()+
	xlab("")+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
	guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
	NULL
ggsave("../results/hs_seperate.pdf", width = 12, height=6)

p3 <- ggplot(df_boot_rst, aes(x=group)) +
	# geom_point(aes())+
	geom_pointrange(aes(y=e_hpi, ymin=ci_hpi_low, ymax=ci_hpi_high, fill=parent_group2), size=0.8, shape=21, stroke=0.6)+
	scale_x_discrete(expand=c(0,1.5), drop=F) +
	ylab(expression("RBD cross-reactivity ("~pi~")"))+
	scale_fill_manual(name="Group", values=colors_t)+
	facet_wrap(vars(RBDs), ncol=2, scales="free_x")+
	theme_bw()+
	xlab("")+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
	guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
	NULL
ggsave("../results/hpi_seperate.pdf", width = 12, height=6)

p_out <- (p0+ggtitle("A"))/(p1+ggtitle("B")) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("../results/Antigenic_response_compare_seperate.pdf", width = 14, height=10)
p_out <- (p0+ggtitle("A"))/(p3+ggtitle("B")) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("../results/Antigenic_response_combined_seperate.pdf", width = 14, height=10)
p_out <- (p1+ggtitle("A"))/(p3+ggtitle("B")) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("../results/Antigenic_response_combined_seperate_2.pdf", width = 14, height=10)

p_out <- (p0+ggtitle("A_1"))/(p1+ggtitle("A_2"))/(p3+ggtitle("B_1"))/(p2+ggtitle("B_2")) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("../results/Antigenic_response_combined_add_hs_seperate.pdf", width = 14, height=18)

## simulate data to illustrate differnece between 
## positive response and all response

num_each_group <- 20
n_RBD <- 4
value_high <- 90
df <- tibble(group=rep(c("Group 1", "Group 2"), each=num_each_group*n_RBD), RBD=rep(rep(c("RBD_1", "RBD_2", "RBD_3", "RBD_4"), each=num_each_group), 2), value=NA)
df$value[df$group=="Group 1" & df$RBD=="RBD_1"] <- rnorm(num_each_group, mean=value_high, sd=5)
df$value[df$group=="Group 1" & df$RBD=="RBD_2"] <- rnorm(num_each_group, mean=0, sd=5)
df$value[df$group=="Group 1" & df$RBD=="RBD_3"] <- rnorm(num_each_group, mean=0, sd=5)
df$value[df$group=="Group 1" & df$RBD=="RBD_4"] <- rnorm(num_each_group, mean=0, sd=5)

expected_mean <- value_high/(n_RBD-1)+10
df$value[df$group=="Group 2" & df$RBD=="RBD_1"] <- rnorm(num_each_group, mean=expected_mean, sd=5)
df$value[df$group=="Group 2" & df$RBD=="RBD_2"] <- rnorm(num_each_group, mean=expected_mean, sd=5)
df$value[df$group=="Group 2" & df$RBD=="RBD_3"] <- rnorm(num_each_group, mean=expected_mean, sd=5)
df$value[df$group=="Group 2" & df$RBD=="RBD_4"] <- rnorm(num_each_group, mean=0, sd=5)

df$Response <- ifelse(df$value<20, "negative", "positive")

ggplot(df, aes(x=RBD, y=value)) +
	geom_boxplot()+
	geom_point(aes(color=Response), size=1, alpha=0.6)+
	# theme_minimal()+
	theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1))+
	scale_x_discrete(guide = guide_axis(n.dodge = 1), expand=c(0,1.5)) +
	ylab("% Inhibition")+
	facet_wrap(vars(group), ncol=2)+
	NULL
ggsave("../results/similated_data_responses.pdf", width=8, height=6)

## compare mean response
df %>% group_by(group) %>% summarise(mean_response=mean(value), mean_positive_response=mean(value[Response=="positive"]))

## compare diversity
df_tmp <- df
df_tmp$RBD2 <- df_tmp$RBD
df_tmp$RBD2[df_tmp$Response=="negative"] <- "Negatives"
df_tmp$RBD2 <- factor(df_tmp$RBD2, levels=c(unique(df_tmp$RBD), "Negatives"))
df_tmp <- df_tmp %>% group_by(group, RBD2) %>% summarise(N=n())
df_tmp_freq <- df_tmp
N_total <- sum(df_tmp_freq$N)
df_tmp_freq$p <- df_tmp_freq$N/N_total
df_tmp_freq %>% group_by(group) %>% summarise(hs=sum(-p*log(p)), hpi=1-(sum(N*(N-1))/(N_total*(N_total-1))))

# df_tmp_freq <- df_tmp_freq %>% filter(group=="Group 1")
# (hs <- sum(-df_tmp_freq$p*log(df_tmp_freq$p)))
# (hpi <- 1-(sum(df_tmp_freq$N*(df_tmp_freq$N-1))/(N_total*(N_total-1))))


# Two dimensional illustration of the "Response fitness"
colors_t <- scico(length(unique(df_boot_rst$parent_group2)), palette = 'batlow')
colors_t <- c("#000000", "#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00", "#a65628","#f781bf","#999999") # colorbrewer2

plot_2d <- function(df, x_var, y_var) {
	# x_var = "mean_pos_res"
	# y_var = "hpi"
	x_lab <- ifelse(grepl("pos", x_var), "Average %inhibition of positive responses", "Average %inhibition of responses")

	x_e <- paste0("e_", x_var)
	y_e <- paste0("e_", y_var)
	x_ci_l <- paste0("ci_", x_var, "_low")
	x_ci_h <- paste0("ci_", x_var, "_high")
	y_ci_l <- paste0("ci_", y_var, "_low")
	y_ci_h <- paste0("ci_", y_var, "_high")
	
	ggplot(df)+
		geom_segment(aes_string(x=x_e, xend=x_e, y=y_ci_l, yend=y_ci_h, color="parent_group2"), alpha=0.6, size=0.5)+
		geom_segment(aes_string(x=x_ci_l, xend=x_ci_h, y=y_e, yend=y_e, color="parent_group2"), alpha=0.6, size=0.5)+
		geom_point(aes_string(x=x_e, y=y_e, color="parent_group2"), alpha=0.8, shape=1, size=2, data=. %>% filter(open_circle))+
		geom_point(aes_string(x=x_e, y=y_e, color="parent_group2", fill="parent_group2"), alpha=0.8, shape=16, size=2, data=. %>% filter(!open_circle), show.legend=FALSE)+
		geom_text_repel(aes_string(x=x_e, y=y_e, color="parent_group2", label="group_rename"), bg.color = "white", bg.r = 0.05, size=2, segment.size=0.3, segment.alpha=0.9, segment.color="black", arrow = arrow(length = unit(0.01, "npc")), point.padding=0.28, box.padding=0.45, show.legend=FALSE, max.overlaps = Inf, force=10)+
		scale_fill_manual(name="Group", values=colors_t)+
		scale_color_manual(name="Group", values=colors_t)+
		ylab(expression("RBD cross-reactivity ("~pi~")"))+
		xlab(x_lab)+
		theme_bw()+
		theme(legend.position = "top")+
		guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
		NULL
}

unique(df_boot_rst$group)
df_boot_rst$open_circle <- FALSE
df_boot_rst$open_circle[grepl("Pre", df_boot_rst$group)] <- TRUE
df_boot_rst$open_circle[grepl("^D0", df_boot_rst$group)] <- TRUE
df_boot_rst$open_circle[grepl("Acute", df_boot_rst$group)] <- TRUE
df_boot_rst$open_circle[grepl("Neg", df_boot_rst$group)] <- TRUE

p_scov2_pos_res <- plot_2d(df_boot_rst %>% filter(RBDs=="SARS-CoV-2"), "mean_pos_res", "hpi")
ggsave("../results/2D_pos_res_hpi_scov2.pdf", width=8, height=6)
p_scov2_all_res <- plot_2d(df_boot_rst %>% filter(RBDs=="SARS-CoV-2"), "mean_res", "hpi")
ggsave("../results/2D_res_hpi_scov2.pdf", width=8, height=6)
p_sarbeco_pos_res <- plot_2d(df_boot_rst %>% filter(RBDs!="SARS-CoV-2"), "mean_pos_res", "hpi")
ggsave("../results/2D_pos_res_hpi_sarbeco.pdf", width=8, height=6)
p_sarbeco_all_res <- plot_2d(df_boot_rst %>% filter(RBDs!="SARS-CoV-2"), "mean_res", "hpi")
ggsave("../results/2D_res_hpi_sarbeco.pdf", width=8, height=6)


p_2d_out_pos_res <- (p_scov2_pos_res+ggtitle("A (SARS-CoV-2)"))/(p_sarbeco_pos_res+ggtitle("B (Sarbecoviruses)")) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("../results/2D_pos_res_hpi_combined.pdf", width = 8, height=10)

p_2d_out_all_res <- (p_scov2_all_res+ggtitle("A (SARS-CoV-2)"))/(p_sarbeco_all_res+ggtitle("B (Sarbecoviruses)")) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("../results/2D_all_res_hpi_combined.pdf", width = 8, height=10)


