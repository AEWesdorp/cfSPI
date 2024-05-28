## -----------------------------------------------------------------------------
DB_DIR="../../databases/"


## -----------------------------------------------------------------------------
SIM_DIR="../../output/cfspi/simulations/"
SIM_STATS_DIR=paste0(SIM_DIR, "results/stats/")
SIM_K2_REP=paste0(SIM_DIR,"results/kraken2_report/after_host_mapping/")
SIM_K2_OUT=paste0(SIM_DIR, "results/kraken2_output/after_host_mapping/")


## -----------------------------------------------------------------------------
SPL_DIR="../../output/cfspi/samples/"
SPL_STATS_DIR=paste0(SPL_DIR, "results/stats/")
SPL_K2_REP=paste0(SPL_DIR,"results/kraken2_report/after_host_mapping/")
SPL_K2_OUT=paste0(SPL_DIR, "results/kraken2_output/after_host_mapping/") 
SPL_K2pb_REP=paste0(SPL_DIR, "results/kraken2_standardDB_report/") 


## -----------------------------------------------------------------------------
SPL_Asp_map="../../output/mapAspergillus/"


## -----------------------------------------------------------------------------
INT_database_stats="../../output/preprocessing/DBstats/"
INT_SPL="../../output/preprocessing/samples/"
INT_ReSeq="../../output/preprocessing/simulations/"
INT_LOSD="../../output/preprocessing/LoSD/"


## -----------------------------------------------------------------------------
dbs <- c('RS_minusT2T','RS','EPRSFv46','EPRSFv46DM','EPRSc2','EPRSFv46MCAspDM','EPRSFv64','EPRSFv64DM','EPRSFv64MCAspDM')
dbs_minT2T <- c('RS','EPRSFv46DM','EPRSc2','EPRSFv46MCAspDM','EPRSFv64DM','EPRSFv64MCAspDM')
dbs_mut <- c('RS w/o CHM13v2','RS','EPRSFv46','EPRSFv46DM','EPRSc2','EPRSFv46MCAspDM','EPRSFv64','EPRSFv64DM','EPRSFv64MCAspDM')
dbs_sel <- c('RS','EPRSc2','EPRSFv64MCAspDM')
dbs_sel_min <- c('EPRSc2','EPRSFv64MCAspDM')

dbs_decon = c('EPRSFv46DM','EPRSc2')
dbs_aug = c('RS','EPRSFv46DM','EPRSFv46MCAspDM','EPRSFv64DM','EPRSFv64MCAspDM')


## -----------------------------------------------------------------------------
output_fig_suppl <- "../../output/suppl_figures/"
output_fig <- "../../output/figures/"


## -----------------------------------------------------------------------------
Karius_AL <- "../resources/Karius_Aspergillus_list.txt"


## -----------------------------------------------------------------------------
fav_threshold="conf0.0"
wright_threshold="conf0.8"