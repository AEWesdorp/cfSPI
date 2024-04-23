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


## -----------------------------------------------------------------------------
genus_db="EPRSFv64MCAspDM"
int_ctrls_fisher_genus <- lst_fisher_genus[[genus_db]] %>% 
    melt(id.vars = c('db','sample','Asp_genus_count','total_read_count','control_type','origin'), 
                 variable.name = "test_type", value.name = "p.value") %>%  
        mutate(taxName = "Aspergillus (G)") %>% 
        filter(test_type == "fisher.t.ext") %>% 
        filter(grepl(sample, pattern = "ctrl$") | grepl(sample, pattern = "ctrl[2-3]")) %>%
        mutate(sample = factor(x = sample, levels = c(rev(unique(as.factor(lst_fisher_genus[[genus_db]]$sample)))))) %>% 
        mutate(sample_short = str_sub(sample, 1, 3)) %>% 
        filter(!is.na(p.value)) %>%  
        filter(sample_short %in% c("A01", "A02", "A03", "A04", "A05")) %>% 
        mutate(sample_short = ifelse(!grepl(sample, pattern = "Pctrl$"), 
                                    yes = paste0(sample_short, sub(sample, pattern = ".*Pctrl", replacement = "-")), 
                                     no = sample_short)) %>% 
        mutate(col_test = ifelse(test_type == "fisher.t.int" & p.value < 0.001, 
                                 yes = "p < 0.001, internal controls",
                          ifelse(test_type == "fisher.t.ext" & p.value < 0.001, 
                                 yes = "p < 0.001, external controls", no = "p > 0.001"))) %>%
        mutate(p.value = -log10(p.value)) %>% 
        ggplot(aes(x = p.value, y = sample_short, col = col_test, shape = taxName)) + 
            geom_vline(xintercept = -log10(0.001), linetype="dotted", color = "black", size=1.5) + 
            geom_hline(yintercept = seq(3.5,6.5,1), linetype = 1 , color = "black", size = 0.15) + 
            geom_point(size = 4, alpha = 0.5) + 
            geom_point(data = . %>% filter(col_test != "p > 0.001"), size = 4, alpha = 0.5) +
            facet_grid(rows = vars(origin), 
                       drop = FALSE, scales = "free_y", space = "free") + 
            ylab("") + xlab("mean p-value Fisher's exact\n-log10") + 
            theme_bw() + labs(color = "", shape = "") + guides(col = "none") +
            scale_color_manual(values = c("p < 0.001, internal control" = colors_mc[6], 
                                          "p < 0.001, external controls" = colors_mc[5], 
                                          "p > 0.001, internal control" = alpha("darkgrey", alpha = 0.2),
                                          "p > 0.001, external controls" = alpha("darkgrey", alpha = 0.2))) + 
            ggtitle("Aspergillus genus level", subtitle = "dREM.260; CT=0.9") + 
            scale_y_discrete(limits=rev) 
#int_ctrls_fisher_genus


## -----------------------------------------------------------------------------
lst_fisher_genus %>% names()


## -----------------------------------------------------------------------------
ext_ctrls_fisher_genus <- lst_fisher_genus[['EPRSFv64MCAspDM']] %>% 
    melt(id.vars = c('db','sample','Asp_genus_count','total_read_count','control_type','origin'), 
                 variable.name = "test_type", value.name = "p.value") %>% 
        mutate(db = ifelse(db == "EPRSc2", yes = "cRE.21", no = as.character(db))) %>% 
        mutate(db = ifelse(db == "EPRSFv64MCAspDM", yes = "dREM.260", no = as.character(db))) %>% 
        mutate(db = factor(db, levels = c("dREM.260", "cRE.21"))) %>% 
        mutate(taxName = "Aspergillus (G)") %>% 
        filter(test_type == "fisher.t.ext") %>% 
        filter(grepl(sample, pattern = "^H")) %>% 
        filter(grepl(sample, pattern = "ctrl$")) %>% 
        mutate(sample = factor(x = sample, levels = c(rev(unique(as.factor(.$sample)))))) %>% 
        mutate(sample_short = str_sub(sample, 1, 3)) %>% 
        filter(!is.na(p.value)) %>%  
        mutate(col_test = ifelse(test_type == "fisher.t.int" & p.value < 0.001, 
                                 yes = "p < 0.001, internal control",
                          ifelse(test_type == "fisher.t.ext" & p.value < 0.001, 
                                 yes = "p < 0.001, external controls", 
                          ifelse(test_type == "fisher.t.int" & p.value > 0.001, 
                                 yes = "p > 0.001, internal control", 
                          ifelse(test_type == "fisher.t.ext" & p.value > 0.001, 
                                 yes = "p > 0.001, external controls", no = NA))))) %>% 
        mutate(p.value = -log10(p.value)) %>% 
        ggplot(aes(x = p.value, y = sample_short, col = col_test, shape = taxName)) + 
            geom_rect(data = . %>% filter(db == "cRE.21"), aes(fill = "red"),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
            geom_vline(xintercept = -log10(0.001), linetype="dotted", color = "black", size=1.5) + 
            geom_point(size = 4, alpha = 0.5) + 
            geom_point(data = . %>% filter(col_test != "p > 0.001"), size = 4, alpha = 0.5) +
            facet_grid(rows = vars(origin), 
                       drop = FALSE, scales = "free_y", space = "free") + 
            ylab("") + xlab("mean p-value Fisher's exact\n-log10") + 
            theme_bw() + labs(color = "", shape = "") + guides(fill = "none", color = "none") +
            scale_color_manual(values = c("p < 0.001, internal control" = colors_mc[6], 
                                          "p < 0.001, external controls" = colors_mc[5], 
                                          "p > 0.001, internal control" = alpha("darkgrey", alpha = 0.2),
                                          "p > 0.001, external controls" = alpha("darkgrey", alpha = 0.2))) + 
            ggtitle("Aspergillus genus level") + 
            scale_y_discrete(limits=rev) 


## -----------------------------------------------------------------------------
species_db="EPRSc2"
int_ctrls_fisher_species <- lst_fisher_species[[species_db]] %>% 
    melt(id.vars = c('taxName','count','db','threshold','sample','control_type','origin'), 
                 variable.name = "test_type", value.name = "p.value") %>% 
        filter(test_type == "fisher.t.ext") %>% 
        filter(grepl(sample, pattern = "ctrl$") | grepl(sample, pattern = "ctrl[2-3]")) %>%
        mutate(sample = factor(x = sample, levels = c(rev(unique(as.factor(lst_fisher_species[[species_db]]$sample)))))) %>% 
        mutate(sample_short = str_sub(sample, 1, 3)) %>% 
        filter(!is.na(p.value)) %>%  
        filter(sample_short %in% c("A01", "A02", "A03", "A04", "A05")) %>% 
        mutate(sample_short = ifelse(!grepl(sample, pattern = "Pctrl$"), 
                                    yes = paste0(sample_short, sub(sample, pattern = ".*Pctrl", replacement = "-")), 
                                     no = sample_short)) %>% 
        mutate(shape_test = taxName) %>% 
        mutate(shape_test = ifelse(p.value < 0.001, 
                                 yes = taxName, no = "Other Aspergillus spp.")) %>% 
        mutate(col_test = ifelse(test_type == "fisher.t.int" & p.value < 0.001, 
                                 yes = "p < 0.001, internal control",
                          ifelse(test_type == "fisher.t.ext" & p.value < 0.001, 
                                 yes = "p < 0.001, external controls", 
                          ifelse(test_type == "fisher.t.int" & p.value > 0.001, 
                                 yes = "p > 0.001, internal control", 
                          ifelse(test_type == "fisher.t.ext" & p.value > 0.001, 
                                 yes = "p > 0.001, external controls", no = NA))))) %>% 
        mutate(p.value = -log10(p.value)) %>% 
        ggplot(aes(x = p.value, y = sample_short, col = col_test, shape = shape_test)) + 
            geom_vline(xintercept = -log10(0.001), linetype="dotted", color = "black", size=1.5) + 
            geom_hline(yintercept = seq(3.5,6.5,1), linetype = 1 , color = "black", size = 0.15) + 
            geom_point(size = 4, alpha = 0.5) +
            geom_point(data = . %>% filter(col_test != "p > 0.001"), size = 4, alpha = 0.5) +
            facet_grid(rows = vars(origin), 
                       drop = FALSE, scales = "free_y", space = "free") + 
            ylab("") + xlab("mean p-value Fisher's exact\n-log10") + 
            theme_bw() + labs(color = "", shape = "") + #guides(col = FALSE, shape = FALSE) + 
            scale_color_manual(values = c("p < 0.001, internal control" = colors_mc[6], 
                                          "p < 0.001, external controls" = colors_mc[5], 
                                          "p > 0.001, internal control" = alpha("darkgrey", alpha = 0.2),
                                          "p > 0.001, external controls" = alpha("darkgrey", alpha = 0.2))) + 
            scale_shape_manual(values = c("Aspergillus fumigatus" = 17, 
                                        "Aspergillus clavatus" = 18, 
                                        "Aspergillus glaucus" = 25, 
                                        "Other Aspergillus spp." = 15), drop = FALSE) + 
            ggtitle("Aspergillus species level", subtitle = "cRE.21; CT=0.4") + 
            scale_y_discrete(limits=rev) + 
            guides(col = guide_legend(override.aes = list(order = 1)),
                    shape = guide_legend(override.aes = list(order = 2))) 


## -----------------------------------------------------------------------------
species_db="EPRSc2" 
ext_ctrls_fisher_species <- lst_fisher_species[["EPRSc2"]] %>% 
    melt(id.vars = c('taxName','count','db','threshold','sample','control_type','origin'), 
                 variable.name = "test_type", value.name = "p.value") %>% 
        mutate(db = ifelse(db == "EPRSc2", yes = "cRE.21", no = as.character(db))) %>% 
        mutate(db = ifelse(db == "EPRSFv64MCAspDM", yes = "dREM.260", no = as.character(db))) %>% 
        filter(test_type == "fisher.t.ext") %>% 
        filter(grepl(sample, pattern = "^H")) %>% 
        filter(grepl(sample, pattern = "ctrl$")) %>% 
        mutate(sample = factor(x = sample, levels = c(rev(unique(as.factor(lst_fisher_species[[species_db]]$sample)))))) %>% 
        mutate(sample_short = str_sub(sample, 1, 3)) %>% 
        filter(!is.na(p.value)) %>%  
        mutate(shape_test = taxName) %>% 
        mutate(shape_test = ifelse(p.value < 0.001, 
                                 yes = taxName, no = "Aspergillus spp.")) %>% 
        mutate(col_test = ifelse(test_type == "fisher.t.int" & p.value < 0.001, 
                                 yes = "p < 0.001, internal control",
                          ifelse(test_type == "fisher.t.ext" & p.value < 0.001, 
                                 yes = "p < 0.001, external controls", 
                          ifelse(test_type == "fisher.t.int" & p.value > 0.001, 
                                 yes = "p > 0.001, internal control", 
                          ifelse(test_type == "fisher.t.ext" & p.value > 0.001, 
                                 yes = "p > 0.001, external controls", no = NA))))) %>% 
        mutate(p.value = -log10(p.value)) %>% 
        ggplot(aes(x = p.value, y = sample_short, col = col_test, shape = shape_test)) + 
            geom_rect(data = . %>% filter(db == "dREM.260") %>% filter(taxName == "Aspergillus fumigatus"), aes(fill = "blue"),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.02) +
            geom_vline(xintercept = -log10(0.001), linetype="dotted", color = "black", size=1.5) + 
            geom_hline(yintercept = seq(1.5,10.5,1), linetype = 1 , color = "black", size = 0.15) + 
            geom_point(size = 4, alpha = 0.5) +
            geom_point(data = . %>% filter(col_test != "p > 0.001"), size = 4, alpha = 0.5) +
            facet_grid(rows = vars(origin), 
                       drop = FALSE, scales = "free_y", space = "free") + 
            ylab("") + xlab("mean p-value Fisher's exact\n-log10") + 
            guides(col = FALSE) + 
            theme_bw() + labs(color = "", shape = "") + guides(fill = "none", col = "none") +
            scale_color_manual(values = c("p < 0.001, internal control" = colors_mc[6], 
                                          "p < 0.001, external controls" = colors_mc[5], 
                                          "p > 0.001, internal control" = alpha("darkgrey", alpha = 0.2),
                                          "p > 0.001, external controls" = alpha("darkgrey", alpha = 0.2))) + 
            scale_shape_manual(values = c("Aspergillus fumigatus" = 17, 
                                        "Aspergillus clavatus" = 18, 
                                        "Aspergillus glaucus" = 25, 
                                        "Aspergillus spp." = 15), drop = FALSE) + 
            ggtitle("Aspergillus species level") + 
            scale_y_discrete(limits=rev) 


## -----------------------------------------------------------------------------
options(repr.plot.width=20, repr.plot.height=16) 
layout <- "AC
           BD"
SuplFig13 <- 
    (int_ctrls_fisher_species + labs(tag = 'a')) + 
    (ext_ctrls_fisher_species + labs(tag = 'b')) + 
    (int_ctrls_fisher_genus + labs(tag = 'c')) + 
    (ext_ctrls_fisher_genus + labs(tag = 'd')) + 
        plot_layout(nrow = 2, heights = c(1.3,4), design = layout) & 
        guides(size = 30) & 
        theme(legend.position='right', legend.justification='top', legend.direction = 'vertical',
            legend.key.size = unit(0.75, 'cm'), legend.key.height = unit(0.75, 'cm'), legend.key.width = unit(0.75, 'cm'), 
            text = element_text(size = 17), legend.text = element_text(size = 15, colour = "black"),
            plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
            plot.subtitle = element_text(size = 17, hjust = 0.5), 
            plot.tag = element_text(face = 'bold', size = 20), 
            strip.background = element_rect(fill="white", color = "white"),
            panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()) &
            xlim(0,20)

SuplFig13

ggsave("../../output/suppl_figures/SuplFig13_R.png", 
       SuplFig13, width = 20, height = 16)
ggsave("../../output/suppl_figures/SuplFig13_R.pdf", 
       SuplFig13, width = 20, height = 16)


## -----------------------------------------------------------------------------
file.exists("Rplots.pdf")
file.remove("Rplots.pdf")