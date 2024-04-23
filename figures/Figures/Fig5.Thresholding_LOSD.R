## -----------------------------------------------------------------------------
# Code accompanying the publication: "cfSPI: NGS-based Aspergillus detection in plasma and lung lavage of children with invasive pulmonary aspergillosis"
# A.E.Wesdorp, L. Rotte et al. 

# Please direct any questions or comments to a.e.wesdorp@umcutrecht.nl


## -----------------------------------------------------------------------------
source("../functions.R")
source("../cols.R")
source("../dirs.R")


## -----------------------------------------------------------------------------
SPL_AspG_RPM <- readRDS(paste0(INT_SPL, "Asp_genus_RPM.Rdata"))
SPL_AspS_RPM <- readRDS(paste0(INT_SPL, "Asp_species_RPM.Rdata"))
lst_perc_Asp <- readRDS(paste0(INT_ReSeq, "lst_perc_Asp.Rdata"))

classification_table <- read.csv(paste0(INT_ReSeq, "classification_table.csv"))

LOSD_genus_table_withNoise <- readRDS(paste0(INT_LOSD, "LOSD_genus_table_withNoise.Rdata"))
LOSD_species_table_withNoise <- readRDS(paste0(INT_LOSD, "LOSD_species_table_withNoise.Rdata"))


## -----------------------------------------------------------------------------
if(exists("RPM_AspG")){rm("RPM_AspG")}
for (n in names(SPL_AspG_RPM)){
    tmp_RPM_AspG <- SPL_AspG_RPM[[n]] %>% 
        melt(varnames = c("db", "sample_id"), value.name = "value") %>%  
        mutate(thrshld = gsub(n, pattern = "conf", replacement = "")) %>% 
        filter(db %in% dbs_sel_min) %>%
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
        filter(origin != "controls") %>% 
        filter(origin != "Internal controls")  %>% 
        mutate(sample_sh = substr(sample_id, 1, 3))  %>% 
        mutate(sample_sh = ifelse(grepl(sample_id, pattern = "VAL"), paste0(sample_sh, "val"), sample_sh)) 
    if(exists("RPM_AspG")){RPM_AspG <- rbind(RPM_AspG, tmp_RPM_AspG)}
    else{RPM_AspG <- tmp_RPM_AspG}
}


## -----------------------------------------------------------------------------
ext_cntrl_AspG <- RPM_AspG %>% 
    mutate(db = ifelse(db == "EPRSFv64MCAspDM", yes = "dREM.260", no = as.character(db))) %>% 
    filter(db == "dREM.260") %>% 
    filter(origin == "External controls") %>% filter(type == "ss-cfDNA") %>% 
    ggplot(aes(x = thrshld, y = value, fill = sample_type)) + 
        #geom_point(position = position_jitterdodge()) +
        geom_boxplot(alpha=0.6) + 
        facet_grid(rows = vars(sample_type)) + 
        ylab("RPM Aspergillus (G)") + xlab("Confidence Threshold") + 
        ggtitle(label = "Genus level", subtitle = "External Controls; dREM.260") + 
        scale_fill_manual(values = c("BAL ss-cfDNA" = colors_mc[3], "plasma ss-cfDNA" = colors_mc[6])) +
        scale_color_manual(values = c("BAL ss-cfDNA" = "darkgrey", "plasma ss-cfDNA" = "black")) +      
        labs(col = "", fill = "") + guides(fill = "none", col = "none")


## -----------------------------------------------------------------------------
RPM_AspS <- SPL_AspS_RPM %>% 
    rename(sample = "sample_id") %>% 
    rename(RPM = "value") %>% 
        mutate(thrshld = gsub(threshold, pattern = "conf", replacement = "")) %>% 
        filter(db %in% dbs_sel_min) %>%
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
        filter(origin != "controls") %>% 
        filter(origin != "Internal controls")  %>% 
        mutate(sample_sh = substr(sample_id, 1, 3))  %>% 
        mutate(sample_sh = ifelse(grepl(sample_id, pattern = "VAL"), paste0(sample_sh, "val"), sample_sh)) 


## -----------------------------------------------------------------------------
ext_cntrl_AspS <- RPM_AspS %>% 
    mutate(db = ifelse(db == "EPRSc2", yes = "cRE.21", no = as.character(db))) %>% 
    filter(db == "cRE.21") %>% 
    mutate(RPM = value*10^6) %>% 
    filter(origin == "External controls") %>% filter(type == "ss-cfDNA") %>% 
    ggplot(aes(x = thrshld, y = RPM, fill = sample_type)) + 
        #geom_point(position = position_jitterdodge()) +
        geom_boxplot(alpha=0.6) + 
        facet_grid(rows = vars(sample_type)) + 
        ylab("RPM Aspergillus (S)") + xlab("Confidence Threshold") + 
        ggtitle(label = "Species level", "External Controls; cRE.21") + 
        scale_fill_manual(values = c("BAL ss-cfDNA" = colors_mc[3], "plasma ss-cfDNA" = colors_mc[6])) +
        scale_color_manual(values = c("BAL ss-cfDNA" = "darkgrey", "plasma ss-cfDNA" = "black")) +      
        labs(col = "", fill = "") + guides(fill = "none", col = "none")


## -----------------------------------------------------------------------------
plt_LOSD <- list()
for (j in 1:2){
    plt_LOSD[[j]] <- list(LOSD_species_table_withNoise, LOSD_genus_table_withNoise)[[j]] %>% 
        as.data.frame() %>% 
        melt(id.vars = colnames(.)[!grepl(colnames(.), pattern = "LOSD")], value.name = "LOSD", variable.name = "control_set") %>%
        mutate(M_reads = as.numeric(gsub(control_set, pattern = "[^0-9.-]+", replacement = ""))) %>% 
        filter(M_reads == 70) %>% 
        select(-strain) %>% 
        mutate(sample = ifelse(grepl(control_set, pattern = "_Bctrls"), yes = 'BAL ss-cfDNA', 
                        ifelse(grepl(control_set, pattern = "_Pctrls"), yes = 'plasma ss-cfDNA', NA))) %>% 
        mutate(sample = factor(sample, levels = c('plasma ss-cfDNA', 'BAL ss-cfDNA'))) %>% 
        group_by(LOSD, db, threshold, sample) %>% 
        summarise(n = n(), .groups = "keep") %>% 
        mutate(MPM = factor(LOSD, levels = rev(c(0.25,0.5,1,2,4,8,16,32,64,128,256,512,1024,2048,4096,">4096"))))  %>% 
        mutate(db = ifelse(db == "EPRSc2", yes = "cRE.21", no = as.character(db))) %>% 
        mutate(db = ifelse(db == "EPRSFv64MCAspDM", yes = "dREM.260", no = as.character(db))) %>% 
        filter(db == c(c("cRE.21", "dREM.260")[j])) %>% 
        ggplot(aes(x = threshold, y = n, fill = MPM)) + 
            geom_col(position = "fill") +
            scale_fill_manual(drop = FALSE, na.value = "black", 
                              breaks = rev(c(0.25,0.5,1,2,4,8,16,32,64,128,256,512,1024,2048,4096,">4096")), 
                              values = c(colors_mc[7], browns(5)[1:4], reds(7)[2:7], blues(6)[1:5])) + 
            guides(color = guide_legend(override.aes = list(fill = NA)),
            linetype = guide_legend(override.aes = list(fill = NA))) +
            #scale_color_manual(values = c("BAL ss-cfDNA" = "darkgrey", "plasma ss-cfDNA" = "black")) +
            scale_y_continuous(breaks=c(0,seq(5, 55, 10))) +
            theme_bw() + theme(legend.key = element_rect(fill = "white")) +
            ylab("Fraction Aspergillus simulations") + 
            ggtitle(label = "", subtitle = "LOSD, 70M reads") + 
            facet_grid(cols = vars(db), rows = vars(sample), scales = "free_x") + xlab("Confidence Threshold") + 
            labs(fill = "Molecules\nPer\nMillion (MPM)") 
}


## -----------------------------------------------------------------------------
options(repr.plot.width=20, repr.plot.height=14)
layout= "
    AC
    BD"
Fig5 <- 
    (ext_cntrl_AspS + labs(tag = 'a')) +  
    (plt_LOSD[[1]] + labs(tag = 'b') + ggtitle("Species level", subtitle = "Limit of Significant Detection; cRE.21")) +
    (ext_cntrl_AspG + labs(tag = 'c')) +  
    (plt_LOSD[[2]] + labs(tag = 'd') + ggtitle("Genus level", subtitle = "Limit of Significant Detection; dREM.260")) +
        plot_layout(ncol = 2, guides = "collect", design = layout, heights = c(1,2)) &
        guides(size = 30) & 
        theme_bw() &
        theme(legend.position='right', legend.justification='bottom', legend.direction = 'vertical',
            legend.key.size = unit(0.75, 'cm'), legend.key.height = unit(0.75, 'cm'), legend.key.width = unit(0.75, 'cm'), 
            text = element_text(size = 17), legend.text = element_text(size = 15, colour = "black"),
            plot.title = element_text(size = 20, face = "bold"), 
            plot.subtitle = element_text(size = 17, hjust = 0.5), 
            strip.text.x = element_blank(), 
            strip.background = element_blank())

Fig5

ggsave("../../output/figures/Fig5_R.png", 
       Fig5, width = 20, height = 14)
ggsave("../../output/figures/Fig5_R.pdf", 
       Fig5, width = 20, height = 14)


## -----------------------------------------------------------------------------
file.exists("Rplots.pdf")
file.remove("Rplots.pdf")