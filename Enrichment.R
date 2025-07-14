library(tidyverse)
library(ggrepel)
library(ggprism)
library(patchwork)
library(ComplexUpset)
library(gprofiler2)

panel <- data.table::fread("panel_info.tsv") %>% mutate(Panel=gsub("_II","",Panel)) %>% filter(!duplicated(Target))
panel$Gene %>% map(~ str_split(.,";")) %>% unlist %>% unique %>% as_tibble %>% data.table::fwrite("background_set_olink.txt", col.names = F)

dat <- data.table::fread("step1_bon_sig_pro.txt", header = F)
excl_pro <- c('ERI1','ACTA2','TRDMT1','ARFIP1','MST1','AARSD1','CD164L2','GIPC2','CA11','TEF','DOK2','TNFRSF6B','MYO9B','ACP1','GSTP1','TIMP4','EFEMP1','ABO','FKBP1B','CEACAM19')
dat <- dat %>% filter(!V1 %in% excl_pro)
dat %>% data.table::fwrite("step1_bon_sig_pro_excel_stringdb.txt", col.names = F)

plot <- data.table::fread("gProfiler_hsapiens_2025-3-14 10-08-45__intersections.csv") %>% filter(highlighted == 'TRUE') %>% mutate(term_name = Hmisc::capitalize(term_name))

plot <- plot %>%
  mutate(gene_list = str_split(as.character(intersections), "[,]"))

gene_term_assoc <- plot %>%
  select(source, term_name, negative_log10_of_adjusted_p_value, gene_list) %>%
  unnest(gene_list) %>%
  rename(gene = gene_list)

trait_assoc <- data.table::fread("result/step1_MR.csv") %>% filter(P_BON < 0.05) %>% select(outcome, exposure)
names(trait_assoc) <- c("trait","gene")

traits <- unique(trait_assoc$trait)
num_traits <- length(traits)

plot$source <- factor(plot$source)
source_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(min(9, 3), "Set2"))(length(unique(plot$source))),
  unique(plot$source)
)

plot <- plot %>%
  arrange(desc(source), negative_log10_of_adjusted_p_value) %>%
  mutate(term_name = factor(term_name, levels = unique(term_name)))

gene_term_assoc$term_name <- factor(gene_term_assoc$term_name, 
                                    levels = levels(plot$term_name))

all_genes <- gene_term_assoc %>%
  group_by(gene) %>%
  summarize(first_appearance = min(which(levels(term_name) == term_name[1]))) %>%
  arrange(first_appearance) %>%
  pull(gene)

p1 <- ggplot(plot, aes(x = negative_log10_of_adjusted_p_value, y = term_name, fill = source)) +
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = source_colors) +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 3)) +
  labs(x = "-log10(FDR q)", y = NULL, fill = NULL) +
  theme_prism(base_fontface = "plain", base_size = 4) +
  theme(axis.text.y = element_text(size = 10),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.5),
        legend.position = "top",
        axis.text.x = element_text(color = "black", size = 10, margin = margin(t = 3)),
        axis.title = element_text(size = 12, colour = "black"),
        axis.ticks.length = unit(rep(0.25, 1),'lines'),
        legend.text = element_text(size = 10),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"))

create_polygon_data <- function(gene_term_assoc, trait_assoc, all_genes, term_levels) {
  polygons <- data.frame()
  
  gene_term_combos <- gene_term_assoc %>%
    select(gene, term_name) %>%
    distinct()
  
  for (i in 1:nrow(gene_term_combos)) {
    gene <- gene_term_combos$gene[i]
    term <- as.character(gene_term_combos$term_name[i])
    
    gene_traits <- trait_assoc %>%
      filter(gene == !!gene) %>%
      pull(trait)
    
    gene_pos <- which(all_genes == gene)
    term_pos <- which(term_levels == term)
    
    if (length(gene_traits) == 0) {
      polygons <- bind_rows(polygons, data.frame(
        gene = gene,
        term_name = term,
        trait = "Unknown",
        x = gene_pos + c(-0.45, 0.45, 0.45, -0.45),
        y = term_pos + c(-0.45, -0.45, 0.45, 0.45),
        id = paste(gene, term, "Unknown", sep = "_"),
        stringsAsFactors = FALSE
      ))
    }
    else if (length(gene_traits) == 1) {
      polygons <- bind_rows(polygons, data.frame(
        gene = gene,
        term_name = term,
        trait = gene_traits[1],
        x = gene_pos + c(-0.45, 0.45, 0.45, -0.45),
        y = term_pos + c(-0.45, -0.45, 0.45, 0.45),
        id = paste(gene, term, gene_traits[1], sep = "_"),
        stringsAsFactors = FALSE
      ))
    }
    else if (length(gene_traits) == 2) {
      polygons <- bind_rows(polygons, data.frame(
        gene = gene,
        term_name = term,
        trait = gene_traits[1],
        x = gene_pos + c(-0.45, 0.45, 0.45),
        y = term_pos + c(-0.45, -0.45, 0.45),
        id = paste(gene, term, gene_traits[1], sep = "_"),
        stringsAsFactors = FALSE
      ))
      
      polygons <- bind_rows(polygons, data.frame(
        gene = gene,
        term_name = term,
        trait = gene_traits[2],
        x = gene_pos + c(-0.45, 0.45, -0.45),
        y = term_pos + c(-0.45, 0.45, 0.45),
        id = paste(gene, term, gene_traits[2], sep = "_"),
        stringsAsFactors = FALSE
      ))
    }
    else {
      n_traits <- length(gene_traits)
      width <- 0.9 / n_traits
      
      for (j in 1:n_traits) {
        left <- -0.45 + (j-1) * width
        right <- -0.45 + j * width
        
        polygons <- bind_rows(polygons, data.frame(
          gene = gene,
          term_name = term,
          trait = gene_traits[j],
          x = gene_pos + c(left, right, right, left),
          y = term_pos + c(-0.45, -0.45, 0.45, 0.45),
          id = paste(gene, term, gene_traits[j], j, sep = "_"),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(polygons)
}
polygon_data <- create_polygon_data(
  gene_term_assoc, 
  trait_assoc, 
  all_genes, 
  levels(plot$term_name)
)

p2 <- polygon_data %>% 
  mutate(trait = factor(trait, levels = c('kdm_delta','phenoage_delta','frailty_score','healthspan','ltl'),
                        labels = c('KDM-BA acceleration','PhenoAge acceleration',
                                   'Frailty index','Healthspan','Leukocyte telomere length'))) %>%
  ggplot(aes(x = x, y = y, fill = trait, group = id)) +
  geom_polygon(color = "white", size = 0.2) +
  #scale_fill_manual(values = trait_colors, drop = F) +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(breaks = 1:length(all_genes), labels = all_genes,
                     expand = c(0, 0), position = "top") +
  scale_y_continuous(breaks = 1:length(levels(plot$term_name)), 
                     labels = levels(plot$term_name), expand = c(0.01, 0.01)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_prism(base_fontface = "plain", base_size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 10),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(0.5,"cm"),
        panel.grid.major.x = element_line(color="grey",linewidth=0.5),
        panel.grid.major.y = element_line(color="grey",linewidth=0.5))

combined_plot <- p1 + p2 + plot_layout(widths = c(1, 4))

ggsave("PLOT3.pdf", plot = combined_plot, width = 18, height = 6.4)


