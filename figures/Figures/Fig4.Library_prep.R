## -----------------------------------------------------------------------------
# Code accompanying the publication: "cfSPI: NGS-based Aspergillus detection in plasma and lung lavage of children with invasive pulmonary aspergillosis"
# A.E.Wesdorp, L. Rotte et al. 

# Please direct any questions or comments to a.e.wesdorp@umcutrecht.nl


## -----------------------------------------------------------------------------
source("../functions.R")
source("../cols.R")
source("../dirs.R")


## -----------------------------------------------------------------------------
Fungi_RPM <- readRDS(paste0(INT_SPL, "Fungi_RPM.Rdata"))
Asp_genus_RPM <- readRDS(paste0(INT_SPL, "Asp_genus_RPM.Rdata"))
Asp_fum_RPM <- readRDS(paste0(INT_SPL, "Asp_fum_RPM.Rdata"))

Fungi_count <- readRDS(paste0(INT_SPL, "Fungi_count.Rdata"))
Asp_genus_count <- readRDS(paste0(INT_SPL, "Asp_genus_count.Rdata"))
Asp_fum_count <- readRDS(paste0(INT_SPL, "Asp_fum_count.Rdata"))
total_read_count <- readRDS(paste0(INT_SPL, "total_read_count.Rdata"))


## -----------------------------------------------------------------------------
Fungal_RPM_df <- Fungi_RPM[[wright_threshold]] %>% 
    melt(varnames = c("db", "sample_id"), value.name = "value") %>%  
    mutate(name = "Fungi") %>% 
    filter(db %in% dbs_sel_min)
Asp_genus_RPM_df <- Asp_genus_RPM[[wright_threshold]] %>% 
    melt(varnames = c("db", "sample_id"), value.name = "value") %>%  
    mutate(name = "Aspergillus") %>% 
    filter(db %in% dbs_sel_min)

df_RPM <- rbind(Fungal_RPM_df) %>%
    mutate(origin = ifelse(grepl(sample_id, pattern = "^A") & grepl(sample_id, pattern = "asp"), yes = "IPA samples",
                    ifelse(grepl(sample_id, pattern = "^A") & grepl(sample_id, pattern = "ctrl"), yes = "Internal controls",
                    ifelse(grepl(sample_id, pattern = "^H"), yes = "External controls", 
                    ifelse(grepl(sample_id, pattern = "^N"), yes = "controls", no = NA)))))  %>% 

    mutate(sample_short = str_replace(sample_id, str_sub(sample_id, 1, 3), "")) %>% 
    mutate(sample = ifelse(grepl(sample_id, pattern = "[0-9]B[a-z]"), yes = "BAL", 
                  ifelse(grepl(sample_id, pattern = "[0-9]P[a-z]"), yes = "plasma", no = "control")))  %>% 
    mutate(type = ifelse(grepl(sample_id, pattern = "K$"), yes = "ds-cfDNA", 
                  ifelse(grepl(sample_id, pattern = "P$"), yes = "ds-wcDNA", no = "ss-cfDNA"))) %>% 
    mutate(sample_type = paste(sample, type))  %>% 
    mutate(origin = fct_relevel(origin, c("IPA samples","Internal controls","External controls", "controls")))  %>% 
    mutate(sample_type = fct_relevel(sample_type, c("plasma ds-cfDNA", "plasma ss-cfDNA", 
                                                    "BAL ds-cfDNA", "BAL ss-cfDNA", 
                                                    "BAL ds-wcDNA")))  %>% 
    mutate(type = fct_relevel(type, c("ss-cfDNA", "ds-cfDNA", "ds-wcDNA"))) %>% 
    filter(origin != "controls")  %>% 
    filter(origin != "Internal controls")  %>% 
    mutate(sample_sh = substr(sample_id, 1, 3))  %>% 
    mutate(sample_sh = ifelse(grepl(sample_id, pattern = "VAL"), paste0(sample_sh, "val"), sample_sh)) %>% 
    filter(sample_id != "A05PaspVAL") %>% 
    mutate(sample_short = str_sub(sample_id, 1, 3)) %>% 
    filter(sample_short %nin% c("A11", "A12", "A13", "A16"))


## -----------------------------------------------------------------------------
my_comp_prep <- list(c("plasma ds-cfDNA", "plasma ss-cfDNA"),
                c("BAL ds-cfDNA", "BAL ss-cfDNA"))
my_comp_specimen <- list(c("BAL ss-cfDNA", "plasma ss-cfDNA"), 
                 c("BAL ds-wcDNA", "plasma ss-cfDNA"),                  
                c("BAL ds-wcDNA", "BAL ss-cfDNA"))

stat.test <- compare_means(as.data.frame(df_RPM), 
                           formula = value~sample_type, group.by = c("db","name", "origin"), 
                           method = "wilcox.test",
                           p.adjust.method = "bonferroni") %>% 
    mutate(db = ifelse(db == "EPRSc2", yes = "cRE.21", no = as.character(db))) %>% 
    mutate(db = ifelse(db == "EPRSFv64MCAspDM", yes = "dREM.260", no = as.character(db)))
stat.test$plot_pvalue <- round(stat.test$p, digits = 3)

stat.test_prep <- stat.test %>%
    filter(origin == "IPA samples") %>% 
    filter(group1 %in% my_comp_prep[[1]] & group2 %in% my_comp_prep[[1]] | 
           group1 %in% my_comp_prep[[2]] & group2 %in% my_comp_prep[[2]] ) %>% 
    filter(db %in% c("cRE.21")) 

stat.test_specimen <- stat.test %>%
    filter(group1 %in% my_comp_specimen[[1]] & group2 %in% my_comp_specimen[[1]] | 
           group1 %in% my_comp_specimen[[2]] & group2 %in% my_comp_specimen[[2]] | 
           group1 %in% my_comp_specimen[[3]] & group2 %in% my_comp_specimen[[3]] ) %>% 
    filter(db %in% c("cRE.21")) 

fig_RPM_prep <- df_RPM %>% 
    filter(origin == "IPA samples") %>% 
    filter(db %in% dbs_sel_min) %>%
    filter(!grepl(sample_type, pattern = "wcDNA")) %>% 
    mutate(db = ifelse(db == "EPRSc2", yes = "cRE.21", no = as.character(db))) %>% 
    mutate(db = ifelse(db == "EPRSFv64MCAspDM", yes = "dREM.260", no = as.character(db))) %>% 
    filter(db == "cRE.21") %>% 
    ggplot(., aes(x=sample_type, y=value, color = type)) + 
        geom_line(aes(x = sample_type, y = value, group = sample_sh), color = "grey", linetype = "dashed") + 
        geom_boxplot(aes(x = sample_type, y=value, fill = type), lwd=0.8, fill = "white") + 
        geom_point(aes(x = sample_type, y=value, col = type), size = 2) + 
        stat_pvalue_manual(stat.test_prep, y.position = rep(87.5), coord.flip = TRUE,
                           label = "p.signif",remove.bracket = FALSE) +
        facet_grid(rows = vars(origin), cols = vars(db), space = "free", scales = "free") + 
        scale_color_manual(values = c("ss-cfDNA" = colors_mc[6], 
                                      "ds-cfDNA" = colors_mc[5], 
                                      "ds-wcDNA" = colors_mc[4])) + 
        theme_bw() +
        theme(strip.background=element_rect(fill="white", color = "white"), 
              strip.text.y.right = element_text(), 
              axis.title.y=element_blank(), legend.title=element_blank()) + 
        ylab("RPM fungi (K)") +
        scale_y_continuous(limits=c(-10, 100), breaks = seq(0,100,50)) +
        geom_text(aes(label=paste0("n=", after_stat(count))), y=-12, stat='count', size=5, hjust = 0) +
        coord_flip()+ ggtitle("ss-cfDNA vs ds-cfDNA library")

fig_RPM_specimen <- df_RPM %>% 
    filter(db %in% dbs_sel_min) %>%
    filter(!grepl(sample_type, pattern = "ds-cfDNA")) %>% 
    mutate(db = ifelse(db == "EPRSc2", yes = "cRE.21", no = as.character(db))) %>% 
    mutate(db = ifelse(db == "EPRSFv64MCAspDM", yes = "dREM.260", no = as.character(db))) %>%      
    filter(db == "cRE.21") %>% 
    ggplot(., aes(x=sample_type, y=value, color = type)) + 
        geom_line(aes(x = sample_type, y = value, group = sample_sh), color = "grey", linetype = "dashed") + 
        geom_boxplot(aes(x = sample_type, y=value, fill = type), lwd=0.8, fill = "white") + 
        geom_point(aes(x = sample_type, y=value, col = type), size = 2) + 
        geom_text(aes(label=paste0("n=", after_stat(count))), y=-12, stat='count', size=5, hjust = 0) +
        scale_color_manual(values = c("ss-cfDNA" = colors_mc[6], 
                                      "ds-cfDNA" = colors_mc[5], 
                                      "ds-wcDNA" = colors_mc[4])) + 
        stat_pvalue_manual(stat.test_specimen, y.position = rep(c(75,87.5,100),2), coord.flip = TRUE,
                           label = "p.signif",remove.bracket = FALSE) +
        facet_grid(rows = vars(origin), cols = vars(db), space = "free", scales = "free") + 
        scale_alpha_manual(values = c(0.3, 1)) +
        theme_bw() +
        theme(strip.background=element_rect(fill="white", color = "white"), 
              strip.text.y.right = element_text(), 
              axis.title.y=element_blank(), legend.title=element_blank()) + 
        ylab("RPM fungi (K)") +
        scale_y_continuous(limits=c(-10, 100), breaks = seq(0,100,50))+
        coord_flip() + ggtitle("wcDNA vs ss-cfDNA sample")

options(repr.plot.width=10, repr.plot.height=10)
Fig4 <- 
    fig_RPM_prep + fig_RPM_specimen + 
        plot_layout(ncol = 1, heights = c(1,1.6)) &
        plot_annotation(tag_levels = 'a') &
        guides(size = 30) & 
        theme(legend.position='right', legend.justification='center', legend.direction = 'vertical',
            legend.key.size = unit(0.75, 'cm'), legend.key.height = unit(0.75, 'cm'), legend.key.width = unit(0.75, 'cm'), 
            text = element_text(size = 17), legend.text = element_text(size = 15, colour = "black"),
            plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
            plot.subtitle = element_text(size = 17, hjust = 0.5))

Fig4

ggsave("../../output/figures/Fig4_R.png", 
       Fig4, width = 10, height = 10)
ggsave("../../output/figures/Fig4_R.pdf", 
       Fig4, width = 10, height = 10)


## -----------------------------------------------------------------------------
file.exists("Rplots.pdf")
file.remove("Rplots.pdf")