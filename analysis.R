library(RColorBrewer)
library(gplots)
library(ggplot2)
library(plyr)
library(rstudioapi)
library(ggrepel)
library(dplyr)
library(multcomp)
library(tidyr)
library(stringr)
library(data.table)
library(plotly)
library(tidyr)
library(tidyselect)
library(readr)
library(uwot)
library(devtools)
library(tidyverse)
library(pheatmap)
library(limma)
library(DataCombine)
library(jcolors)
library(broom)


data.path = "M:/nas/21_Chris/Projects/Collaborations/1906_pSILAC_HBV_Kaechele_Protzer/main experiment/without LFQ normalization/analysis/R"
plots.path = "M:/nas/21_Chris/Projects/Collaborations/1906_pSILAC_HBV_Kaechele_Protzer/main experiment/without LFQ normalization/analysis/R/figures"


df <- as.data.frame(read_tsv(file.path(data.path, "proteinGroups.txt"), guess_max = 10000, col_names = TRUE))

colnames(df) <- gsub(" ", ".", colnames(df))

df.sel <- df %>%
  filter(is.na(Only.identified.by.site) & is.na(Reverse) & is.na(Potential.contaminant)) %>%
  select(starts_with("LFQ"), Gene.names, Majority.protein.IDs, id) #%>% select(sort(names(.)))

df.long <- df.sel %>%
  rename_all(funs(str_replace(., "LFQ.intensity.H.", "H_"))) %>%
  rename_all(funs(str_replace(., "LFQ.intensity.L.", "L_"))) %>%
  pivot_longer(., cols = c(starts_with("H"), starts_with("L")), names_to = "mschannel", values_to = "intensity") %>%
  mutate(log2_intensity = ifelse(is.finite(log2(intensity)), log2(intensity), NA_real_)) %>%
  separate(mschannel, into = c("mstag", "infection", "timepoint", "SILAC", "replicate"), sep="_") %>%
  mutate(mschannel = paste0(mstag, "_", infection, "_", timepoint, "_", replicate)) %>%
  mutate(msrun = paste0(infection, "_", timepoint, "_", replicate)) %>%
  select(-SILAC) %>%
  filter(timepoint != "0") %>%
  mutate(id.gene.names = str_trunc(paste0(id, "_", Gene.names), 20)) %>% 
  mutate(infection = factor(infection, levels = c("mock", "HBV", "HBVdX")),
         timepoint = factor(timepoint, levels = c("1", "2", "3", "4", "6", "8", "10"))) %>%
  arrange(mstag, infection, timepoint, replicate, id) %>%
  mutate(mschannel = factor(mschannel, levels = unique(mschannel)),
         msrun = factor(msrun, levels = unique(msrun)))


#boxplot of non-normalized data
ggplot(df.long, aes(x=mschannel, y=log2_intensity)) +
  geom_boxplot(size=0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))


#raw heatmap of non-normalized, imputed data
df.long.imput <- df.long %>%
  group_by(mschannel) %>%
  summarize(quant_mean = mean(log2_intensity, na.rm = TRUE),
            quant_sd = sd(log2_intensity, na.rm = TRUE)) %>%
  ungroup()

hm.df.long <- df.long %>%
  left_join(df.long.imput) %>%
  mutate(log2_intensity_imput = if_else(!is.na(log2_intensity), log2_intensity, 
                                                 rnorm(n(), mean = quant_mean - 1.8 * quant_sd, sd = 0.3 * quant_sd))) %>%
  select(-quant_mean, -quant_sd)

hm.df.wide <- hm.df.long %>% 
  select(id.gene.names, mschannel, log2_intensity_imput) %>% 
  pivot_wider(id_cols = "id.gene.names", names_from = "mschannel", values_from = "log2_intensity_imput")

row.names(hm.df.wide) <- hm.df.wide$id.gene.names
hm.df.wide <- select(hm.df.wide, -id.gene.names)

pheatmap(hm.df.wide, cluster_cols = F, show_rownames = F, fontsize_col = 5)


#raw heatmap of only row (L) normalized, imputed values
#impute
df.long.imput <- df.long %>%
  group_by(mschannel) %>%
  summarize(quant_mean = mean(log2_intensity, na.rm = TRUE),
            quant_sd = sd(log2_intensity, na.rm = TRUE)) %>%
  ungroup()

hm.df.long <- df.long %>%
  left_join(df.long.imput) %>%
  mutate(log2_intensity_imput = if_else(!is.na(log2_intensity), log2_intensity, 
                                        rnorm(n(), mean = quant_mean - 1.8 * quant_sd, sd = 0.3 * quant_sd))) %>%
  select(-quant_mean, -quant_sd)

#row normalize
sum.hm.row.norm <- hm.df.long %>% 
  filter(mstag == "L") %>% 
  group_by(id.gene.names) %>% 
  summarise(log2_intensity_row_median = median(log2_intensity, na.rm = T)) %>%
  ungroup() %>% 
  select(id.gene.names, log2_intensity_row_median)

hm.df.wide <- hm.df.long %>% 
  left_join(sum.hm.row.norm) %>% 
  mutate(log2_intensity_imput_row_norm = log2_intensity_imput - ifelse(!is.na(log2_intensity_row_median), log2_intensity_row_median, (hm.df.long %>% 
                                                                                                                            group_by(id.gene.names) %>% 
                                                                                                                            mutate(log2_intensity_row_median = median(log2_intensity_imput[mstag=="L"], na.rm = TRUE)) %>%
                                                                                                                            ungroup())$log2_intensity_row_median)) %>%
  select(id.gene.names, mschannel, log2_intensity_imput_row_norm) %>% 
  pivot_wider(id_cols = "id.gene.names", names_from = "mschannel", values_from = "log2_intensity_imput_row_norm")

row.names(hm.df.wide) <- hm.df.wide$id.gene.names
hm.df.wide <- select(hm.df.wide, -id.gene.names)

pheatmap(hm.df.wide, cluster_cols = F, show_rownames = F, fontsize_col = 5)



##row and col normalize within replicates

#caluclate row median within replicates
sum.row.norm.4rep <- df.long %>%
  filter(mstag == "L") %>%
  group_by(id.gene.names, infection, timepoint) %>%
  mutate(log2_intensity_row_median = median(log2_intensity, na.rm = TRUE)) %>%
  ungroup() %>%
  select(id.gene.names, infection, timepoint, log2_intensity_row_median)

#subtract row median within replicates
df.row.norm <- df.long %>%
  filter(mstag == "L") %>% 
  left_join(sum.row.norm.4rep) %>%
  mutate(log2_intensity_row_norm = log2_intensity - log2_intensity_row_median) %>%
  select(-log2_intensity_row_median)


#column normalize df within replicates

#calculate column normalization factors for each msrun
sum.rep.norm <- df.row.norm %>%
  group_by(msrun, infection, timepoint, replicate) %>%
  summarise(log2_intensity_rep_median = median(log2_intensity_row_norm, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(infection, timepoint) %>%
  mutate(log2_intensity_rep_median_shift = log2_intensity_rep_median - median(log2_intensity_rep_median, na.rm = TRUE)) %>%
  ungroup() %>% 
  select(msrun, log2_intensity_rep_median_shift)
  

#plot distribution of normalization factors
ggplot(sum.rep.norm) +
  geom_density(aes(x=log2_intensity_rep_median_shift))

#apply normalization factor
df.rep.norm <- df.long %>%
  left_join(sum.rep.norm) %>%
  mutate(log2_intensity_rep_norm = log2_intensity - log2_intensity_rep_median_shift) %>%
  select(-log2_intensity_rep_median_shift)


#calculate row median of replicate normalized values of L channel
sum.rep.row.norm <- df.rep.norm %>%
  filter(mstag=="L") %>%
  group_by(id.gene.names) %>%
  summarise(log2_intensity_row_median = median(log2_intensity_rep_norm, na.rm = TRUE)) %>%
  ungroup() %>%
  select(id.gene.names, log2_intensity_row_median)


#impute
sum.rep.norm.imput.df <- df.rep.norm %>%
  group_by(mschannel) %>%
  summarize(quant_mean = mean(log2_intensity_rep_norm, na.rm = TRUE),
            quant_sd = sd(log2_intensity_rep_norm, na.rm = TRUE)) %>%
  ungroup()

df.rep.norm.imput <- df.rep.norm %>%
  left_join(sum.rep.norm.imput.df) %>%
  mutate(log2_intensity_rep_norm_imput = if_else(!is.na(log2_intensity_rep_norm), log2_intensity_rep_norm, 
                                                 rnorm(n(), mean = quant_mean - 1.8 * quant_sd, sd = 0.3 * quant_sd))) %>%
  select(-quant_mean, -quant_sd)


##umap of imputed df
#row normalize

df.rep.row.norm.imput <- df.rep.norm.imput %>%
  left_join(sum.rep.row.norm) %>%
  mutate(log2_intensity_rep_norm_imput = log2_intensity_rep_norm_imput - 
           ifelse(!is.na(log2_intensity_row_median), log2_intensity_row_median, (df.rep.norm.imput %>% 
                                                                                   group_by(id.gene.names) %>% 
                                                                                   mutate(log2_intensity_rep_norm_imput_median = median(log2_intensity_rep_norm_imput[mstag == "L"], na.rm = TRUE)) %>%
                                                                                   ungroup())$log2_intensity_rep_norm_imput_median)) %>%
  select(-log2_intensity_row_median) %>%
  pivot_wider(., id_cols = "id.gene.names", names_from = "mschannel", values_from = "log2_intensity_rep_norm_imput")

# define umap as matrix which only allows numeric values. w/o umap would also try to cluster id and gene names and not only numerical values
umap.mtx <-  df.rep.row.norm.imput %>%
  select(starts_with("L_"),
         starts_with("H_")
  ) %>% as.matrix()

row.names(umap.mtx) <- df.rep.row.norm.imput$id.gene.names

#calculate umap
umap <- as_tibble(umap(umap.mtx))

#select data to color,... umap
df.rep.norm.imput.sel <- df.rep.norm.imput %>%
  filter(
    mstag == "L" &
      infection == "mock" &
      timepoint == "6" &
      replicate == "1"
  )

umap.sel <- umap %>%
  mutate(id.gene.names = df.rep.row.norm.imput$id.gene.names) %>%
  left_join(df.rep.norm.imput.sel)


#draw umap
ggplot(umap.sel, aes(x=V1, y=V2, 
                     color=log2_intensity_rep_norm
)) +
  geom_point() +
  geom_point(data=filter(umap.sel, str_detect(id.gene.names, "_HIST")), color="red") +
geom_label_repel(data=filter(umap.sel, str_detect(id.gene.names, "RPL")), aes(label = id.gene.names))



#heatmap of umap data
pheatmap(umap.mtx[,], cluster_cols = F, show_rownames = F, fontsize_col = 5)

#boxplot of umap data (already rep and row normalized)
df.rep.row.norm.imput.boxplot <- df.rep.row.norm.imput %>%
  pivot_longer(cols = c(starts_with("L"), starts_with("H")), names_to = "mschannel", values_to = "log2_intensity_rep_norm_imput")

ggplot(df.rep.row.norm.imput.boxplot, aes(x=mschannel, y=log2_intensity_rep_norm_imput)) +
  geom_boxplot(size=0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))

#boxplot of umap data w/o row norm, but rep norm
ggplot(df.rep.norm, aes(x=mschannel, y=log2_intensity_rep_norm)) +
  geom_boxplot(size=0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))



##select stable proteins (high L/H ratio, high L, low H intensity)

df.rep.norm.imput.mod <- df.rep.norm.imput %>%
  select(-intensity) %>%
  group_by(id.gene.names, msrun) %>%
  mutate(log2_LHratio = log2_intensity_rep_norm[mstag=="L"] - log2_intensity_rep_norm[mstag=="H"]) %>%
  ungroup()

sel.stab.prots <- df.rep.norm.imput.mod %>%
  group_by(id.gene.names) %>%
  filter(length(unique(msrun[!is.na(log2_intensity_rep_norm[mstag=="L"])])) >= 80 &
           median(log2_LHratio, na.rm = TRUE) >= 2) %>%
  ungroup() %>%
  select(id)

df.stab.prots <- df.rep.norm.imput.mod %>%
  filter(id %in% unique(sel.stab.prots$id))

length(unique(df.stab.prots$id.gene.names))


#draw heatmap of selected stab. proteins

df.stab.prots.imput.wide <- df.stab.prots %>%
  left_join(sum.rep.row.norm) %>% 
  mutate(log2_intensity_rep_row_norm_imput = log2_intensity_rep_norm_imput - log2_intensity_row_median) %>%
  select(-log2_intensity_row_median) %>%
  filter(mstag == "L") %>%
  select(id.gene.names, msrun, log2_intensity_rep_row_norm_imput) %>%
  pivot_wider(id_cols = "id.gene.names", names_from = "msrun", values_from = "log2_intensity_rep_row_norm_imput")

rownames(df.stab.prots.imput.wide) <- df.stab.prots.imput.wide$id.gene.names

method_dist = "euclidean"
method_hclust = "ward.D2"

data.dist <- dist(select(df.stab.prots.imput.wide, -id.gene.names), method = method_dist)
data.hclust <- hclust(data.dist, method = method_hclust)

breaksList = seq(-1.5, 1.5, by = 0.1)

pheatmap(select(df.stab.prots.imput.wide, -id.gene.names), color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                breaks = breaksList, cluster_rows = data.hclust, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE)

pheatmap(select(df.stab.prots.imput.wide, -id.gene.names),cluster_rows = data.hclust, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE)

cairo_pdf(file.path(plots_path, paste0("stab-prots_", method_dist, "_", method_hclust, ".pdf")), width = 9, height = 8)
print(p)
dev.off()



#plot stab. prots for all conditions
ggplot(data=df.stab.prots, aes(x=timepoint, y=log2_intensity_rep_norm, color=str_trunc(id.gene.names, 20), group=id)) +
  geom_point() +
  geom_smooth(aes(fill=str_trunc(id.gene.names, 20)), alpha=0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(~ mstag + infection)

#loop for each id.gene.names
id.gene.names.list <- unique(df.stab.prots$id.gene.names)

for (i in id.gene.names.list){
  ggplot(data=filter(df.stab.prots, id.gene.names == i), aes(x=timepoint, y=log2_intensity_rep_norm, color=id.gene.names)) +
    geom_point() +
    geom_smooth(aes(fill=str_trunc(id.gene.names, 10)), alpha=0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(~ mstag + infection) +
    scale_y_continuous(limits = c(20, 40)) +
    ggsave(file.path(plots.path, paste0("profile_plot_", i, ".pdf")), width = 11, height = 8)
}


#plot stab. prots for specific condition with several id.gene.names per page
df.stab.prots.sel <- df.stab.prots %>% 
  filter(mstag == "L" & infection == "HBVdX",
         id.gene.names %in% id.gene.names.list,
         #id.gene.names %in% id.gene.names.list[11:20],
         #id.gene.names %in% id.gene.names.list[21:30],
         #id.gene.names %in% id.gene.names.list[31:40],
         #id.gene.names %in% id.gene.names.list[41:49]
         )

ggplot(df.stab.prots.sel, aes(x=timepoint, y=log2_intensity_rep_norm, color=id.gene.names)) +
  geom_point() +
  geom_smooth(aes(fill=str_trunc(id.gene.names, 10)), alpha=0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  facet_wrap(~ id.gene.names) #+
  scale_y_continuous(limits = c(27, 38))
  
  
##final list of selected stable proteins
id.gene.names.list <- unique(df.stab.prots$id.gene.names)

stab.prots.notselected <- c("3744_HIST2H3A", "4087_GHDC", "5880_OXSM", "2191_HMGA2", "2615_DHODH", "5586_CARS2", "953_TOMM40", "1801_GLRX", "5944_FKBP11", "1399_HMGA1", "1733_GK;GK3P")

id.gene.names.df <- as.data.frame(id.gene.names.list)

colnames(id.gene.names.df) <- "id.gene.names"

stab.prots.final.list <- id.gene.names.df %>%
  filter(!(id.gene.names %in% stab.prots.notselected))

df.stab.prots <- df.rep.norm.imput %>%
  filter(id.gene.names %in% stab.prots.final.list$id.gene.names)

##calculate msrun shift based on stab proteins
sum.stab.prots <- df.stab.prots %>% 
  filter(mstag == "L") %>% 
  group_by(id.gene.names) %>%
  mutate(log2_intensity_rep_row_norm = log2_intensity_rep_norm - median(log2_intensity_rep_norm, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(msrun) %>% 
  summarise(log2_msrun_median_stab_prots = median(log2_intensity_rep_row_norm, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(log2_msrun_shift_stab_prots = log2_msrun_median_stab_prots - median(log2_msrun_median_stab_prots, na.rm = T)) %>% 
  select(msrun, log2_msrun_shift_stab_prots)
  
ggplot(sum.stab.prots) +
  geom_density(aes(x=log2_msrun_shift_stab_prots)) 
  
df.stab.prots.norm <- df.rep.norm.imput %>% 
  left_join(sum.stab.prots) %>% 
  mutate(log2_intensity_stabprots_norm_imput = log2_intensity_rep_norm_imput - log2_msrun_shift_stab_prots)
  
sum.stab.prots.norm.row.norm <- df.stab.prots.norm %>% 
  filter(mstag == "L") %>% 
  group_by(id.gene.names) %>% 
  summarise(log2_intensity_stabprots_norm_imput_row_median = median(log2_intensity_stabprots_norm_imput, na.rm = T)) %>% 
  ungroup() %>% 
  select(id.gene.names, log2_intensity_stabprots_norm_imput_row_median)
  
df.stab.prots.norm.row.norm <- df.stab.prots.norm %>% 
  left_join(sum.stab.prots.norm.row.norm) %>% 
  mutate(log2_intensity_stabprots_row_norm_imput = log2_intensity_stabprots_norm_imput - log2_intensity_stabprots_norm_imput_row_median)

df.stab.prots.norm.row.norm.wide <- df.stab.prots.norm.row.norm %>% 
  select(id.gene.names, mschannel, log2_intensity_stabprots_row_norm_imput) %>% 
  pivot_wider(id_cols = "id.gene.names", names_from = "mschannel", values_from = "log2_intensity_stabprots_row_norm_imput")

#make umap
umap.norm.mtx <-  df.stab.prots.norm.row.norm.wide %>%
  select(starts_with("L_"),
         starts_with("H_")
  ) %>% as.matrix()

row.names(umap.norm.mtx) <- df.stab.prots.norm.row.norm.wide$id.gene.names

umap <- as_tibble(umap(umap.norm.mtx))

#select data to color,... umap
df.rep.norm.imput.sel <- df.stab.prots.norm %>%
  filter(
    mstag == "H" &
      infection == "mock" &
      timepoint == "6" &
      replicate == "1"
  )

umap.sel <- umap %>%
  mutate(id.gene.names = df.rep.row.norm.imput$id.gene.names) %>%
  left_join(df.rep.norm.imput.sel)


#draw umap
ggplot(umap.sel, aes(x=V1, y=V2, 
                     color=log2_intensity_stabprots_norm_imput
)) +
  geom_point(size=0.5) +
  geom_point(data=filter(umap.sel, str_detect(id.gene.names, "_HIST")), color="red") +
  geom_point(data=filter(umap.sel, str_detect(id.gene.names, "_MRP")), color="green") +
  geom_point(data=filter(umap.sel, str_detect(id.gene.names, "_RPL")), color="darkgreen") +
  geom_label_repel(data=filter(umap.sel, str_detect(id.gene.names, "_HIST") & V1 < 9), aes(label = id.gene.names)) +
  scale_colour_jcolors_contin(palette = "pal12")

#heatmap of umap data
pheatmap(umap.norm.mtx[,], cluster_cols = F, show_rownames = F, fontsize_col = 5)


#boxplot of rep and stab_prot normalized data
ggplot(df.stab.prots.norm, aes(x=mschannel, y=log2_intensity_stabprots_norm_imput)) +
  geom_boxplot(size=0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))



### GLM

#create new conditions
df.4glm <- df.stab.prots.norm %>%
  mutate(mstag = factor(mstag, levels = c("L", "H")),
         time = as.integer(as.character(timepoint)),
         condition = paste0(mstag, "_", infection, "_", timepoint),
         after1 = factor(time >= 1),
         after2 = factor(time >= 2),
         after3 = factor(time >= 3),
         after4 = factor(time >= 4),
         after6 = factor(time >= 6),
         after8 = factor(time >= 8),
         after10 = factor(time >= 10))

#model.matrix

exp.des <- select(df.4glm, condition, mstag, infection, timepoint, starts_with("after")) %>% distinct()
exp.des.m <- model.matrix(~(mstag:(after2+after3+after4+after6+after8+after10))*infection + mstag + mstag:infection - infection,
                          data=exp.des)
rownames(exp.des.m) <- exp.des$condition
exp.des.m[,"infectionHBV"] <- exp.des.m[,"infectionHBV"] - exp.des.m[,"mstagH:infectionHBV"]
exp.des.m[,"infectionHBVdX"] <- exp.des.m[,"infectionHBVdX"] - exp.des.m[,"mstagH:infectionHBVdX"]
exp.des.m <- exp.des.m[, str_detect(colnames(exp.des.m), "after2FALSE|^infection\\w+$", negate = TRUE)]
colnames(exp.des.m) <- if_else(str_detect(colnames(exp.des.m), "^infection"),
                               str_c("mstagL:", colnames(exp.des.m)),
                               colnames(exp.des.m))
#det(exp.des.m)
#heatmap.2(exp.des.m, Rowv = F, Colv = F, margins = c(10,10))
pheatmap(exp.des.m, cluster_rows = F, cluster_cols = F, border_color = "grey")

#compute glm
df.glm <- df.4glm %>% #filter(id <= 100) %>%
  group_by(id, id.gene.names) %>% 
  #filter(n_distinct(mstag) >= 1 & n_distinct(timepoint) >= 5 & n_distinct(infection) >= 2) %>%
  do({
    pg_df <- .
    if (nrow(pg_df) > 0) {
      pg_glm <- glm(log2_intensity_stabprots_norm_imput ~ (mstag:(after2+after3+after4+after6+after8+after10))*infection + mstag + mstag:infection - infection, 
                    data = mutate(pg_df, id.gene.names = factor(id.gene.names)))
      tidy(pg_glm)
      #pg_glm <- glm.fit(exp.des.m[pg_df$condition,], pg_df$log2_intensity_rep_norm_imput)
      #pg_smry <- summary(pg_glm)
      #pg_smry_df <- as_tibble(pg_smry$coefficients, rownames="effect")
      #colnames(pg_smry_df) <- c("effect", "estimate", "std", "t_value", "p_value")
      #pg_smry_df
    } else {
      tibble()
    }
  })


###glm with H/L ratio instead of H intensity effect
#model.matrix

names(df.4glm)

df.4glm.HLratio <- df.4glm %>% 
  pivot_wider(id_cols = c("id.gene.names", "msrun"), names_from = "mstag",
              values_from = "log2_intensity_stabprots_norm_imput") %>%
  mutate(HL = H - L) %>% 
  pivot_longer(cols = c("L", "H", "HL"), names_to = "mstag", values_to = "log2_intensity_stabprots_norm_imput")

df.4glm.HLratio.mod <- df.4glm.HLratio %>% 
  filter(mstag != "H") %>% 
  left_join(select(df.4glm, starts_with("after"), msrun, id.gene.names, infection, timepoint)) %>% 
  arrange(mstag, infection, timepoint) %>% 
  mutate(mstag = factor(mstag, levels = c("L", "HL"))) %>% 
  mutate(condition = paste0(mstag, "_", infection, "_", timepoint)) %>% 
  distinct()

exp.des <- select(df.4glm.HLratio.mod, condition, mstag, infection, timepoint, starts_with("after")) %>% distinct()
exp.des.m <- model.matrix(~(mstag:(after2+after3+after4+after6+after8+after10))*infection + mstag + mstag:infection - infection,
                          data=exp.des)
rownames(exp.des.m) <- exp.des$condition
#exp.des.m[,"infectionHBV"] <- exp.des.m[,"infectionHBV"] - exp.des.m[,"mstagHL:infectionHBV"]
#exp.des.m[,"infectionHBVdX"] <- exp.des.m[,"infectionHBVdX"] - exp.des.m[,"mstagHL:infectionHBVdX"]
#exp.des.m <- exp.des.m[, str_detect(colnames(exp.des.m), "after2FALSE|^infection\\w+$", negate = TRUE)]
colnames(exp.des.m) <- if_else(str_detect(colnames(exp.des.m), "^infection"),
                               str_c("mstagL:", colnames(exp.des.m)),
                               colnames(exp.des.m))
#det(exp.des.m)
#heatmap.2(exp.des.m, Rowv = F, Colv = F, margins = c(10,10))
pheatmap(exp.des.m, cluster_rows = F, cluster_cols = F, border_color = "grey")


#compute glm
df.glm <- df.4glm.HLratio.mod %>% #filter(id <= 100) %>%
  group_by(id.gene.names) %>% 
  #filter(n_distinct(mstag) >= 1 & n_distinct(timepoint) >= 5 & n_distinct(infection) >= 2) %>%
  do({
    pg_df <- .
    if (nrow(pg_df) > 0) {
      pg_glm <- glm(log2_intensity_stabprots_norm_imput ~ (mstag:(after2+after3+after4+after6+after8+after10))*infection + mstag + mstag:infection - infection, 
                    data = mutate(pg_df, id.gene.names = factor(id.gene.names)))
      tidy(pg_glm)
      #pg_glm <- glm.fit(exp.des.m[pg_df$condition,], pg_df$log2_intensity_rep_norm_imput)
      #pg_smry <- summary(pg_glm)
      #pg_smry_df <- as_tibble(pg_smry$coefficients, rownames="effect")
      #colnames(pg_smry_df) <- c("effect", "estimate", "std", "t_value", "p_value")
      #pg_smry_df
    } else {
      tibble()
    }
  })
#####


#plotting all effects
df.glm <- df.glm %>%
  mutate(log10_pvalue = -log10(p.value))

#View(filter(df.glm, id==2431))

# df.glm.stats <- group_by(df.glm, term) %>%
#   summarise(estimate_mean = mean(estimate),
#             estimate_sd = sd(estimate),
#             estimate_median = median(estimate),
#             n_signif_up = sum(p.value <= 0.01 & estimate >= 0.1),
#             n_signif_down = sum(p.value <= 0.01 & estimate <= -0.1)) %>%
#   ungroup()

ggplot(df.glm, aes(x=estimate, y=log10_pvalue)) +
  geom_vline(xintercept = c(2,-2)) +
  geom_hline(yintercept = 2) +
  geom_point(alpha=0.1) +
  facet_wrap(~term, scales = "free") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 5))
#geom_label_repel(data=filter(df.glm.sel, abs(estimate) > 1.1 & log10_pvalue > 2), aes(label=id.gene.names))

#plotting selected effect
df.glm.sel <- df.glm %>%
  mutate(log10_pvalue = -log10(p.value)) %>%
  filter(term == "mstagL:after8TRUE:infectionHBV")

ggplot(df.glm.sel, aes(x=estimate, y=log10_pvalue)) +
  geom_vline(xintercept = c(2,-2)) +
  geom_hline(yintercept = 2) +
  geom_point(alpha=0.5) +
  geom_label_repel(data=filter(df.glm.sel, abs(estimate) > 2.5 & log10_pvalue > 2), aes(label=id.gene.names))


###annotate significant glm hits
df.glm.sig <- df.glm %>% 
  mutate(term = str_replace_all(term, ":", "-")) %>%
  #mutate(significant = ifelse((abs(estimate) > 2 & log10_pvalue > 2), TRUE, FALSE)) %>%
  mutate(significant_all = ifelse((abs(estimate) > 2 & log10_pvalue > 2), term, "")) %>% 
  #pivot_wider(id_cols = "id.gene.names", values_from = c("significant"), names_from = "term", names_prefix = "sig_") #%>% 
  left_join(select(df.long, id.gene.names, Majority.protein.IDs)) %>%
  distinct() %>% 
  select(id.gene.names, Majority.protein.IDs, significant_all) %>% 
  group_by(id.gene.names, Majority.protein.IDs) %>% 
  summarise(significant_all_merged = gsub('(;)+',';',paste0(significant_all, collapse = ";")))

write_tsv(df.glm.sig, file.path(data.path, "df.glm.sig_HL-intensity.txt"), col_names = T)


### merge GLM results with umap

umap.norm.mtx <-  df.stab.prots.norm.row.norm.wide %>%
  select(starts_with("L_"),
         starts_with("H_")
  ) %>% as.matrix()

row.names(umap.norm.mtx) <- df.stab.prots.norm.row.norm.wide$id.gene.names

umap <- as_tibble(umap(umap.norm.mtx))


### HC of umap results to define clusters for GO enrichments
umap.4HC <- umap %>%
  mutate(id.gene.names = df.rep.row.norm.imput$id.gene.names)

umap.4HC.mod <- select(umap.4HC, V1, V2)

row.names(umap.4HC.mod) <- umap.4HC$id.gene.names

hclust.results <- hclust(dist(umap.4HC.mod, method = "euclidean"), method = "ward.D2")
pheatmap(umap.4HC.mod, cluster_rows = hclust.results)

umap.HC.cluster <- cutree(hclust.results, k=40)
umap.HC.cluster.mod <- as.data.frame(umap.HC.cluster) %>% 
  rownames_to_column(., var = "id.gene.names") %>% 
  left_join(select(df.long, id.gene.names, Majority.protein.IDs)) %>%
  distinct()

write_tsv(umap.HC.cluster.mod, file.path(data.path, "umap.HC.cluster.mod.txt"), col_names = T)


# display cluster on umap with GO terms

umap.HC.cluster.mod_GOanot <- read_tsv(file.path(data.path, "umap.HC.cluster.mod_GOanot.txt"), col_names = T)

umap.sel <- umap %>%
  mutate(id.gene.names = df.rep.row.norm.imput$id.gene.names) %>%
  #left_join(df.glm) %>%
  left_join(select(umap.HC.cluster.mod_GOanot, id.gene.names, umap.HC.cluster, `GOBP name`, `GOCC name`, `GOMF name`))
  #filter(term == "mstagL:infectionHBV")

colnames(umap.sel) <- gsub(" ", ".", colnames(umap.sel))


GO.term <- "mitochondrion"
umap.cluster <- "14"

p <- ggplot(umap.sel, aes(x=V1, y=V2, 
                       color=as.integer(umap.HC.cluster), 
                       #color=pmin(3, pmax(-3, estimate)), 
                       #size=log10_pvalue
                       ), shape=1) +
    ggtitle(paste0(GO.term, "_cluster:", umap.cluster)) +
    geom_point(aes(text=id.gene.names)) +
    #scale_size_continuous(range = c(0.3, 3)) +
    #geom_point(data=filter(umap.sel, str_detect(id.gene.names, "_HIST")), color="red") +
    #geom_point(data=filter(umap.sel, str_detect(id.gene.names, "_MRP")), color="green") +
    #geom_point(data=filter(umap.sel, str_detect(id.gene.names, "_RPL")), color="darkgreen") +
    geom_label_repel(data=filter(umap.sel, str_detect(GOCC.name, GO.term) & umap.HC.cluster == umap.cluster), aes(label = id.gene.names)) +
    scale_colour_jcolors_contin(palette = "pal12") +
    theme_bw()

show(p)

cairo_pdf(file.path(plots.path, paste0("umap_glm_", umap.cluster, "_", GO.term, ".pdf")), width = 11.6, height = 8.2)
print(p)
dev.off()
###




###umap of L/H ratios
names(df.stab.prots.norm)
df.stab.prots.norm.LHratio <- df.stab.prots.norm %>% 
  pivot_wider(id_cols = c("id.gene.names", "msrun"), values_from = "log2_intensity_stabprots_norm_imput", names_from = "mstag") %>% 
  mutate(log2_intensity_stabprots_norm_imput_LHratio = H - L) %>% 
  group_by(id.gene.names) %>% 
  mutate(log2_intensity_stabprots_norm_imput_LHratio_row_norm = log2_intensity_stabprots_norm_imput_LHratio - median(log2_intensity_stabprots_norm_imput_LHratio)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = "id.gene.names", values_from = "log2_intensity_stabprots_norm_imput_LHratio", names_from = "msrun") 

umap.LHratio.mtx <- df.stab.prots.norm.LHratio %>% 
  select(starts_with("mock"),
         starts_with("HBV"),
         starts_with("HBVdX")) %>% 
  as.matrix()

row.names(umap.LHratio.mtx) <- df.stab.prots.norm.LHratio$id.gene.names

umap <- as_tibble(umap(umap.LHratio.mtx))

ggplot(umap, aes(x=V1, y=V2)) +
  geom_point()
###


###umap of H intensity only
df.stab.prots.norm.H <- df.stab.prots.norm %>% 
  filter(mstag == "H") %>% 
  group_by(id.gene.names) %>% 
  mutate(log2_intensity_stabprots_norm_imput = log2_intensity_stabprots_norm_imput - median(log2_intensity_stabprots_norm_imput, na.rm = T)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c("id.gene.names", "mschannel"), values_from = "log2_intensity_stabprots_norm_imput", names_from = "mschannel")
  
umap.H.mtx <- df.stab.prots.norm.H %>% 
  select(starts_with("H_")) %>% 
  as.matrix()

row.names(umap.LHratio.mtx) <- df.stab.prots.norm.H$id.gene.names

umap <- as_tibble(umap(umap.LHratio.mtx))
###



##ploty umap+glm for single glm effects

umap.sel <- umap %>%
  mutate(id.gene.names = df.rep.row.norm.imput$id.gene.names) %>%
  left_join(df.glm) %>%
  filter(term == "mstagL:infectionHBV")

ggplotly(
ggplot(umap.sel, aes(x=V1, y=V2, 
                     color=estimate, 
                     #color=pmin(3, pmax(-3, estimate)), 
                     size=log10_pvalue)) +
  ggtitle(paste0(umap.sel$term)) +
  geom_point(aes(text=id.gene.names)) +
  scale_size_continuous(range = c(0.3, 3)) +
  #geom_point(data=filter(umap.sel, str_detect(id.gene.names, "_HIST")), color="red") +
  #geom_point(data=filter(umap.sel, str_detect(id.gene.names, "_MRP")), color="green") +
  #geom_point(data=filter(umap.sel, str_detect(id.gene.names, "_RPL")), color="darkgreen") +
  #geom_label_repel(data=filter(umap.sel, str_detect(id.gene.names, "_HIST") & V1 < 9), aes(label = id.gene.names)) +
  scale_colour_jcolors_contin(palette = "pal12") +
  theme_bw(),
  tooltip = "text"
)



## plot umap+glm in loop for all glm effects
for (i in unique(df.glm$term)){

umap.sel <- umap %>%
  mutate(id.gene.names = df.rep.row.norm.imput$id.gene.names) %>%
  left_join(df.glm) %>%
  filter(term == i)

cairo_pdf(file.path(plots.path, paste0("umap_glm-HLratio_", str_replace_all(i, ":", "-"), ".pdf")), width = 11.6, height = 8.2)

p<-ggplot(umap.sel, aes(x=V1, y=V2, 
                     color=estimate,
                     #color=pmin(2, pmax(-2, estimate)), 
                     size=log10_pvalue)) +
  ggtitle(paste0(i)) +
  geom_point() +
  scale_size_continuous(range = c(0.3, 3)) +
  #geom_point(data=filter(umap.sel, str_detect(id.gene.names, "_HIST")), color="red") +
  #geom_point(data=filter(umap.sel, str_detect(id.gene.names, "_MRP")), color="green") +
  #geom_point(data=filter(umap.sel, str_detect(id.gene.names, "_RPL")), color="darkgreen") +
  #geom_label_repel(data=filter(umap.sel, str_detect(id.gene.names, "_HIST") & V1 < 9), aes(label = id.gene.names)) +
  scale_colour_jcolors_contin(palette = "pal12") +
  theme_bw()

print(p)
dev.off()
}

###normalize stabprot/rep normalized intensities that all 3 infection conditions are starting from the same intensity at timepoint 1d
#on rep norm df by defining a set of proteins which are stable between infections
names(df.rep.norm.imput)
df.rep.norm.imput.tp1d.stats.pre <- df.rep.norm.imput %>% 
  filter(mstag == "L", timepoint == "1") %>%
  group_by(infection, id.gene.names) %>% 
  summarise(log2_intensity = median(log2_intensity_rep_norm_imput)) %>%
  ungroup()

df.rep.norm.imput.tp1d.stats.pairs <- left_join(
  df.rep.norm.imput.tp1d.stats.pre,
  df.rep.norm.imput.tp1d.stats.pre,
  by=c("id.gene.names")) %>%
  filter(as.character(infection.x) < as.character(infection.y)) %>%
  mutate(log2_intensity_delta = log2_intensity.x - log2_intensity.y) %>%
  group_by(infection.x, infection.y) %>%
  mutate(delta_percent = percent_rank(log2_intensity_delta)) %>% 
  ungroup()

df.rep.norm.imput.tp1d.stats.genes <- df.rep.norm.imput.tp1d.stats.pairs %>%
  group_by(id.gene.names) %>%
  summarise(low_percent = min(delta_percent),
            high_percent = max(delta_percent)) %>%
  ungroup() %>%
  mutate(is_used = low_percent >= 0.25 & high_percent <= 0.75)

df.rep.norm.imput.tp1d.stats <- df.rep.norm.imput.tp1d.stats.genes %>% 
  filter(is_used == TRUE) %>% 
  left_join(df.rep.norm.imput.tp1d.stats.pre) %>% 
  group_by(id.gene.names) %>% 
  mutate(log2_intensity_row_norm = log2_intensity - median(log2_intensity, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(infection) %>% 
  summarise(log2_intensity_row_norm_median = mean(log2_intensity_row_norm)) %>% 
  ungroup() %>% 
  mutate(log2_intensity_row_norm_median_shift = log2_intensity_row_norm_median - median(log2_intensity_row_norm_median))
  
df.rep.norm.imput.tp1d.norm <- df.rep.norm.imput %>% 
  left_join(df.rep.norm.imput.tp1d.stats) %>% 
  mutate(log2_intensity_rep_norm_imput_tp1d_norm = log2_intensity_rep_norm_imput - log2_intensity_row_norm_median_shift) %>% 
  mutate(is.detected = ifelse(is.na(log2_intensity), FALSE, TRUE))





### line plot of proteins of interest using norm. intensities

infection_palette = c("mock" = "gray", "HBV" = "firebrick", "HBVdX" = "orange")
is.detected_palette = c("TRUE" = 16, "FALSE" = 1)

sel_data.df <- filter(df.rep.norm.imput.tp1d.norm, str_detect(Gene.names, "^APP")) %>% 
  mutate(timepoint = as.numeric(as.character(timepoint)))

ggplot(sel_data.df, aes(x = timepoint, y = log2_intensity_rep_norm_imput_tp1d_norm,
                        color=infection, fill=infection, group=infection)) +
  geom_smooth(alpha=0.5) +
  geom_point(data=sel_data.df, aes(x = timepoint, y = log2_intensity_rep_norm_imput_tp1d_norm,
                              color=infection, fill=infection, shape=is.detected)) +
  facet_grid(mstag ~ id.gene.names, scales = "free") +
  scale_color_manual(values = infection_palette) +
  scale_fill_manual(values = infection_palette) +
  scale_shape_manual(values = is.detected_palette) +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 6, 8, 10)) +
  theme_bw()

###


### line plot of proteins of interest using LH ratios

df.stab.prots.norm.LHratio.lineplot <- df.stab.prots.norm %>% 
  pivot_wider(id_cols = c("id.gene.names", "msrun", "timepoint", "infection"), values_from = "log2_intensity_stabprots_norm_imput", names_from = "mstag") %>% 
  mutate(log2_intensity_stabprots_norm_imput_LHratio = H - L)

infection_palette = c("mock" = "gray", "HBV" = "firebrick", "HBVdX" = "orange")

sel_data.df <- filter(df.stab.prots.norm.LHratio.lineplot, str_detect(id.gene.names, "SQSTM1")) %>% 
  mutate(timepoint = as.numeric(as.character(timepoint)))

ggplot(sel_data.df, aes(x = as.integer(timepoint), y = log2_intensity_stabprots_norm_imput_LHratio,
                        color=infection, fill=infection)) +
  geom_smooth(alpha=0.5) +
  geom_point() +
  facet_grid(~id.gene.names, scales = "free_y") +
  scale_color_manual(values = infection_palette) +
  scale_fill_manual(values = infection_palette) +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 6, 8, 10)) +
  theme_bw()

###

