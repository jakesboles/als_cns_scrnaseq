library(Seurat)
library(scCustomize)
library(tidyverse)
library(lme4)
library(emmeans)
library(multcomp)
library(UCell)
# library(org.Hs.eg.db)
library(msigdbr)
library(clusterProfiler)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

cell <- "Oligodendrocyte"

results_dir <- paste0("results/wgcna_consensus/", cell, "/")

scores <- read.csv(paste0(results_dir, "module_scores_ucell.csv"))

scores <- scores %>% 
  mutate(id = str_split_i(orig.ident, "_", i = 1))

cols <- colnames(scores)[str_detect(colnames(scores), "UCell_kNN")]
color <- str_remove_all(cols, "_UCell_kNN")

stats <- list()
fit <- list()
emm <- list()

for (i in seq_along(cols)){
  
  message(color[i])
  
  df <- scores %>%
    dplyr::rename("active_col" = cols[i])
  
  fit[[i]] <- lme4::lmer(active_col ~ group * tissue + (1|id/orig.ident),
                    data = df)
  suppressMessages({
    print(joint_tests(fit[[i]]))
    
    emm[[i]] <- emmeans(fit[[i]], pairwise ~ group | tissue)
  })
  
  stats[[i]] <- multcomp::cld(emm[[i]], Letters = letters) %>%
    mutate(.group = str_remove_all(.group, " ")) %>%
    mutate(module = color[i])
}

saveRDS(fit, 
        file = paste0(results_dir, "lmer_stats.rds"))
saveRDS(emm,
        file = paste0(results_dir, "lmer_emm.rds"))

stats_df <- list_rbind(stats)

stats_df %>%
  mutate(group = factor(group,
                        levels = c("Control", "sALS", "C9orf72")),
         tissue = factor(tissue,
                         levels = c("Motor cortex", "Cervical spinal cord"))) %>%
  ggplot(aes(x = tissue,
             y = emmean)) + 
  geom_crossbar(aes(ymin = asymp.LCL,
                    ymax = asymp.UCL,
                    fill = group),
                position = position_dodge(width = 1)) + 
  geom_text(aes(label = .group,
                y = (asymp.UCL + emmean) / 2,
                group = group),
            position = position_dodge(width = 1)) + 
  scale_fill_manual(values = c("dodgerblue1", "magenta1", "chartreuse1")) + 
  facet_wrap(. ~ module,
             scales = "free_y") +
  # ggtitle(cell_type) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"))
ggsave(filename = paste0(results_dir, "lmer_stats_overview.png"),
       units = "in", dpi = 600,
       height = length(cols)*.7,
       width = length(cols)*.7)

# GSEA on modules ---------------------------------------------------------

m_t2g <- msigdbr(species = "Homo sapiens",
                 category = "C2")
m_t2g <- m_t2g %>%
  filter(gs_subcollection %in% c("CP:BIOCARTA", "CP:KEGG_LEGACY", "CP:REACTOME", "CP:WIKIPATHWAYS", "CP:PID")) %>%
  dplyr::select(c(gs_name, gene_symbol))
# unique(m_t2g$gs_name)

go_t2g <- msigdbr(species = "Homo sapiens",
                  category = "C5")
go_t2g <- go_t2g %>%
  # filter(gs_subcollection %in% c("GO:BP", "GO:CC") %>%
  dplyr::select(c(gs_name, gene_symbol))

t2g <- rbind(go_t2g, m_t2g)

modules <- read.csv(paste0(results_dir, "modules.csv"))

background <- modules$gene_name %>% unique()

colors <- unique(modules$module)
colors <- colors[colors != "grey"]

for (j in seq_along(colors)){
  
  message(paste0(str_to_title(colors[j]), " module"))
  
  genes <- modules %>%
    filter(module == colors[j]) %>%
    pull(gene_name)
  
  em <- enricher(genes,
                 TERM2GENE = t2g,
                 universe = background)
  
  if (em@result$p.adjust %>% min() < 0.05){
    p <- dotplot(em,
                 showCategory = 20,
                 x = "FoldEnrichment",
                 color = "p.adjust",
                 size = "GeneRatio")
    ggsave(p,
           filename = paste0(results_dir, colors[j], ".png"),
           units = "in", dpi = 600,
           height = 10, width = 8)
  }
  
  write.csv(as.data.frame(em@result),
            file = paste0(results_dir, colors[j], "_gsea.csv"))
}
