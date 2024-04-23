## -----------------------------------------------------------------------------
# Code accompanying the publication: "cfSPI: NGS-based Aspergillus detection in plasma and lung lavage of children with invasive pulmonary aspergillosis"
# A.E.Wesdorp, L. Rotte et al. 

# Please direct any questions or comments to a.e.wesdorp@umcutrecht.nl


## -----------------------------------------------------------------------------
source("../functions.R")
source("../cols.R")
source("../dirs.R")

## -----------------------------------------------------------------------------
lst_fisher_genus <- readRDS(paste0(INT_SPL, "lst_fisher_genus.Rdata"))
lst_fisher_species <- readRDS(paste0(INT_SPL, "lst_fisher_species.Rdata"))
Asp_genus_RPM <- readRDS(paste0(INT_SPL, "Asp_genus_RPM.Rdata"))
Asp_species_RPM <- readRDS(paste0(INT_SPL, "Asp_species_RPM.Rdata"))


## -----------------------------------------------------------------------------
ann_NA <- data.frame(sample_short = c("A15","A14"), p.value = c(125,125), label = "not available",
                       origin = factor("Plasma ss-cfDNA", levels = c("BAL ss-cfDNA", "Plasma ss-cfDNA")))

ann_NA_2 <- data.frame(sample_short = c("A15","A14","A15","A14"), p.value = 150, label = "not available",
                       origin = factor("Plasma ss-cfDNA", levels = c("BAL ss-cfDNA", "Plasma ss-cfDNA")), 
                      db = c("dREM.260", "dREM.260", "cRE.21", "cRE.21"))


## -----------------------------------------------------------------------------
dt_fisher_species <- lst_fisher_species[["EPRSc2"]] %>% 
    melt(id.vars = c('taxName','count','db','threshold','sample','control_type','origin'), 
                 variable.name = "test_type", value.name = "p.value") %>% 
        filter(grepl(sample, pattern = "asp$")) %>% 
        mutate(sample = factor(x = sample, levels = c(rev(unique(as.factor(lst_fisher_species[[1]]$sample)))))) %>% 
        mutate(sample_short = str_sub(sample, 1, 3)) %>% 
        filter(!is.na(p.value)) %>%  
        filter(sample_short %nin% c("A11", "A12", "A13", "A16")) %>% 
        mutate(shape_test = ifelse(taxName == "Aspergillus fumigatus", 
                                 yes = taxName, no = "Other species")) %>% 
        mutate(col_test = ifelse(test_type == "fisher.t.int" & p.value < 0.001, 
                                 yes = "p < 0.001, internal control",
                          ifelse(test_type == "fisher.t.ext" & p.value < 0.001, 
                                 yes = "p < 0.001, external controls", 
                          ifelse(test_type == "fisher.t.int" & p.value > 0.001, 
                                 yes = "p > 0.001, internal control", 
                          ifelse(test_type == "fisher.t.ext" & p.value > 0.001, 
                                 yes = "p > 0.001, external controls", no = NA))))) %>% 
        filter(grepl(col_test, pattern = "external control")) %>% 
        mutate(p.value = -log10(p.value)) 

plt_fisher_species <- dt_fisher_species %>% 
        ggplot(aes(x = p.value, y = sample_short, col = col_test)) + 
            geom_hline(yintercept = 2.5, linetype = 1 , color = "black", linewidth = 1) + 
            geom_hline(yintercept = seq(1.5,6.5,1), linetype = 1 , color = "black", linewidth = 0.15) + 
            geom_vline(xintercept = -log10(0.001), linetype="dotted", color = "black", linewidth=1.5) + 
            geom_point(aes(shape = shape_test), size = 4, alpha = 0.5) + 
            geom_point(aes(shape = shape_test), data = . %>% filter(grepl(col_test, pattern = "p < 0.001")), size = 4, alpha = 1) +
            facet_grid(cols = vars(origin), 
                       drop = TRUE, scales = "free_y", space = "free") + 
            ylab("") + xlab("mean p-value Fisher's exact\n-log10") + 
            theme_bw() + labs(color = "Fisher's exact (mean)", shape = "Species") + 
            scale_color_manual(values = c("p < 0.001, internal control" = colors_mc[6], 
                                          "p < 0.001, external controls" = colors_mc[5], 
                                          "p > 0.001, internal control" = "darkgrey",
                                          "p > 0.001, external controls" = "darkgrey")) + 
            scale_shape_manual(values = c(17,15)) + 
            scale_y_discrete(limits=rev) +
            guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) + 
            geom_text(data = ann_NA, label = "not available", col = "black", size = 5, fontface="italic") +
            geom_text(data = dt_fisher_species %>% filter(p.value == Inf) %>% 
                      mutate(p.value = 225), label = "Inf", col = colors_mc[5], size = 5) + xlim(0,250) + 
            ggtitle("Aspergillus spp. identification", subtitle = "cRE.21; CT=0.4")


## -----------------------------------------------------------------------------
df_fisher_genus <- lst_fisher_genus[["EPRSFv64MCAspDM"]] %>% 
    melt(id.vars = c('db','sample','Asp_genus_count','total_read_count','control_type','origin'), 
                 variable.name = "test_type", value.name = "p.value") %>%   
        filter(grepl(sample, pattern = "asp$")) %>% 
        mutate(sample = factor(x = sample, levels = c(rev(unique(as.factor(lst_fisher_genus[[1]]$sample)))))) %>% 
        mutate(sample_short = str_sub(sample, 1, 3)) %>% 
        filter(!is.na(p.value)) %>%  
        filter(sample_short %nin% c("A11", "A12", "A13", "A16")) %>% 
        
        mutate(col_test = ifelse(test_type == "fisher.t.int" & p.value < 0.001, 
                                 yes = "p < 0.001, internal control",
                          ifelse(test_type == "fisher.t.ext" & p.value < 0.001, 
                                 yes = "p < 0.001, external controls", 
                          ifelse(test_type == "fisher.t.int" & p.value > 0.001, 
                                 yes = "p > 0.001, internal control", 
                          ifelse(test_type == "fisher.t.ext" & p.value > 0.001, 
                                 yes = "p > 0.001, external controls", no = NA))))) %>% 
        filter(grepl(col_test, pattern = "external control")) %>% 
        mutate(p.value = -log10(p.value)) %>% 
        mutate(shape_test = "Aspergillus")  

plt_fisher_genus <- df_fisher_genus %>% 
        ggplot(aes(x = p.value, y = sample_short, col = col_test)) + 
            geom_hline(yintercept = 2.5, linetype = 1 , color = "black", linewidth = 0.4) + 
            geom_hline(yintercept = seq(1.5,6.5,1), linetype = 1 , color = "black", linewidth = 0.15) + 
            geom_vline(xintercept = -log10(0.001), linetype="dotted", color = "black", linewidth=1.5) + 
            geom_point(aes(shape = shape_test), size = 4, alpha = 0.5) + 
            geom_point(aes(shape = shape_test), data = . %>% filter(grepl(col_test, pattern = "p < 0.001")), size = 4, alpha = 1) +
            facet_grid(cols = vars(origin), 
                       drop = TRUE, scales = "free_y", space = "free") + 
            ylab("") + xlab("mean p-value Fisher's exact\n-log10") + 
            theme_bw() + labs(color = "Fisher's exact (mean)", shape = "Genus") + 
            scale_color_manual(values = c("p < 0.001, internal control" = colors_mc[6], 
                                          "p < 0.001, external controls" = colors_mc[5], 
                                          "p > 0.001, internal control" = "darkgrey",
                                          "p > 0.001, external controls" = "darkgrey")) + 
            scale_y_discrete(limits=rev) +
            geom_text(data = df_fisher_genus %>% filter(p.value == Inf) %>% 
                      mutate(p.value = 225), label = "Inf", col = colors_mc[5], size = 5) +
            geom_text(data = ann_NA, label = "not available", col = "black", size = 5, fontface="italic") +
            xlim(0,250) + 
            ggtitle("Aspergillus (G) identification", subtitle = "dREM.260; CT=0.9") + 
            guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))


## -----------------------------------------------------------------------------
species_RPM_plots <- Asp_species_RPM %>% 
        filter(grepl(threshold, pattern = str_extract(0.4, pattern = "[0-9].[0-9]"))) %>% 
        filter(taxName == "Aspergillus fumigatus") %>% 
        melt(varnames = c("db", "sample"), id.vars = c("taxName", "db", "threshold", "sample"), value.name = "value") %>%  
        mutate(sample_short = str_sub(sample, 1, 3)) %>% 
        filter(grepl(sample, pattern = "asp$")) %>% 
        filter(db == "EPRSc2") %>% 
        mutate(sample = ifelse(grepl(sample, pattern = "[0-9]B[a-z]"), yes = "BAL", 
                ifelse(grepl(sample, pattern = "[0-9]P[a-z]"), yes = "plasma", no = "control"))) %>% 
        filter(sample_short %nin% c("A11", "A12", "A13", "A16")) %>% 
        mutate(value = value*10^6) %>% 
        mutate(RPM = "fractional abundance") %>% 
        mutate(origin = ifelse(sample == "plasma", yes = "Plasma ss-cfDNA",
                        ifelse(sample == "BAL", yes = "BAL ss-cfDNA", no = NA))) %>% 
        ggplot(aes(x = value, y = sample_short)) +
            geom_segment(aes(x=0, xend=value, y=sample_short, yend=sample_short, color = RPM), size = 2, alpha = 0.7) +
            geom_point(aes(color = RPM), size = 4,shape = 17) + 
            scale_color_manual(values = colors_mc[[4]],
                              guide = guide_legend(override.aes = list(shape = c(NA)))) + 
            facet_grid(cols = vars(origin)) + 
            scale_y_discrete(limits=rev) + 
            xlab("RPM Aspergillus fumigatus") + ylab("") + theme_bw() +
            geom_hline(yintercept = 2.5, linetype = 1 , color = "black", size = 0.4) +
            geom_hline(yintercept = seq(1.5,6.5,1), linetype = 1 , color = "black", size = 0.15) + 
            labs(color = "") + xlim(0,60) + 
            geom_text(data = ann_NA_2 %>% mutate(value = 30), label = "not available", col = "black", size = 5, fontface="italic")


## -----------------------------------------------------------------------------
genus_RPM_plots <- list()
for (fav_db in dbs_sel_min){
    if (fav_db == "EPRSc2"){print(fav_db); tmp_genus_threshold = "conf0.8"}
    if (fav_db == "EPRSFv64MCAspDM"){print(fav_db); tmp_genus_threshold = "conf0.9"}
                
    print(tmp_genus_threshold)
    genus_RPM_plots[[fav_db]] <- Asp_genus_RPM[[tmp_genus_threshold]] %>% 
        melt(varnames = c("db", "sample"), value.name = "value") %>%  
        mutate(sample_short = str_sub(sample, 1, 3)) %>% 
        filter(grepl(sample, pattern = "asp$")) %>% 
        filter(db == fav_db) %>% 
        mutate(sample = ifelse(grepl(sample, pattern = "[0-9]B[a-z]"), yes = "BAL", 
                ifelse(grepl(sample, pattern = "[0-9]P[a-z]"), yes = "plasma", no = "control"))) %>% 
        filter(sample_short %nin% c("A11", "A12", "A13", "A16")) %>% 
        mutate(RPM = "fractional abundance") %>% 
        mutate(origin = ifelse(sample == "plasma", yes = "Plasma ss-cfDNA",
                        ifelse(sample == "BAL", yes = "BAL ss-cfDNA", no = NA))) %>% 
        ggplot(aes(x = value, y = sample_short)) +
            geom_segment(aes(x=0, xend=value, y=sample_short, yend=sample_short, color = RPM), size = 2, alpha = 0.7) +
            geom_point(aes(color = RPM), size = 4,shape = 19) + 
            scale_color_manual(values = colors_mc[[4]],
                              guide = guide_legend(override.aes = list(shape = c(NA)))) + 
            facet_grid(cols = vars(origin)) + 
            scale_y_discrete(limits=rev) + 
            xlab("RPM Aspergillus (G)") + ylab("") + theme_bw() +
            geom_hline(yintercept = 2.5, linetype = 1 , color = "black", size = 0.4) +
            geom_hline(yintercept = seq(1.5,6.5,1), linetype = 1 , color = "black", size = 0.15) + 
            labs(color = "") + xlim(0,60) + 
            geom_text(data = ann_NA_2 %>% mutate(value = 30), label = "not available", col = "black", size = 5, fontface="italic")
}


## -----------------------------------------------------------------------------
options(repr.plot.width=20, repr.plot.height=10) 
layout <- "AC
           BD"
Fig6 <- plt_fisher_species + 
    species_RPM_plots + 
    plt_fisher_genus + 
    genus_RPM_plots[["EPRSFv64MCAspDM"]] + 
        plot_layout(ncol = 2,design = layout, widths = c(3,3,2)) &
        plot_annotation(tag_levels = 'a') &
        theme(legend.position='right', legend.justification='top', legend.direction = 'vertical',
            legend.key.size = unit(0.75, 'cm'), legend.key.height = unit(0.75, 'cm'), legend.key.width = unit(0.75, 'cm'), 
            text = element_text(size = 17), legend.text = element_text(size = 15, colour = "black"),
            plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
            plot.subtitle = element_text(size = 17, hjust = 0.5),
            strip.background=element_rect(fill="white", color = "white"), 
            panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()) 

Fig6

ggsave("../../output/figures/Fig6_R.png", 
       Fig6, width = 20, height = 10)
ggsave("../../output/figures/Fig6_R.pdf", 
       Fig6, width = 20, height = 10)


## -----------------------------------------------------------------------------
file.exists("Rplots.pdf")
file.remove("Rplots.pdf")