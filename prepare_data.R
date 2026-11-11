# loading and preparing the HBV pulseSILAC proteome data (C.Urban)
#
# Author: Alexey Stukalov
###############################################################################

project_id <- 'curban_hbv_psilac'
message('Project ID=', project_id)
data_version <- "20200705"
fit_version <- "20200705"
ms_folder <- 'mqoutput_20200705'
message('Dataset version is ', data_version)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

source(file.path(misc_scripts_path, 'fasta_utils.R'))
source(file.path(misc_scripts_path, 'matrix_utils.R'))
source(file.path(misc_scripts_path, "ggplot_ext.R"))

require(maxquantUtils)
require(msglm)
require(dplyr)
require(rjson)
require(stringr)
require(readr)
require(pheatmap)
require(tidyr)

msdata_path <- file.path(data_path, ms_folder)

data_info <- list(project_id = project_id,
                  data_ver = data_version, fit_ver = fit_version,
                  ms_folder = ms_folder,
                  instr_calib_protgroup_filename = "instr_protgroup_LFQ_calib_scaturro_qep5calib_20161110_borg.json",
                  instr_calib_pepmodstate_filename = "instr_QEP5_intensity_pepmodstate_calib_cov2_20200430.json",
                  quant_type = "LFQ", quant_col_prefix = "LFQ_Intensity",
                  pep_quant_type = "intensity")

message('Loading MS instrument calibration data from ', data_info$instr_calib_filename, '...')
instr_calib_protgroup <- fromJSON(file = file.path(data_path, data_info$instr_calib_protgroup_filename))$instr_calib
instr_calib_pepmodstate <- fromJSON(file = file.path(data_path, data_info$instr_calib_pepmodstate_filename))$instr_calib
instr_calib <- instr_calib_pepmodstate

bad_msruns <- c()

msruns.df <- read_tsv(file.path(msdata_path, "combined", "experimentalDesign.txt"),
                      col_names=TRUE, col_types = list(Fraction="i")) %>%
    dplyr::rename(raw_file=Name, msfraction=Fraction, msrun=Experiment, is_ptm=PTM) %>%
    extract(msrun, c("treatment", "timepoint", "channels", "replicate"),
            c("^(.+)_(\\d+)_(.+)_(\\d+)$"), remove = FALSE) %>%
    dplyr::mutate(replicate = parse_integer(replicate), timepoint_num = parse_double(timepoint)) %>%
    mutate(treatment = relevel(factor(treatment), "mock"),
           timepoint = factor(timepoint_num),
           condition = str_c(treatment, "_", timepoint, "h"),
           sample = str_c(condition, "_", replicate),
           is_used = !(msrun %in% bad_msruns))

fasta.dfs <- list(
    hbv = read_innate_uniprot_fasta(file.path(data_path, "fasta", "HBV_wt_and_delX_protzer.fasta")),
    human = read_innate_uniprot_fasta(file.path(data_path, "fasta", "uniprot-9606_proteome_human_reviewed_canonical_isoforms_191008.fasta"))
)

msdata.wide <- read.MaxQuant.ProteinGroups(file.path(msdata_path, 'combined/txt'),
                                           file_name = "proteinGroups.txt",
                                           import_data = c(data_info$quant_type, "ident_type", "ms2_count"))
mqevidence <- read.MaxQuant.Evidence(file.path(msdata_path, 'combined/txt'), 
                                     evidence.pepobj = "pepmodstate", correct_ratios = FALSE)
mqevidence$peptides <- read.MaxQuant.Peptides(file.path(msdata_path, 'combined/txt'),
                                              file_name = "peptides.txt",
                                              import_data = c("ident_type"))
mqevidence$peaks <- NULL # exclude big data frame

mqrdata_filepath <- file.path(scratch_path, str_c(project_id, '_mqdata_', ms_folder, '.RData'))
message('Saving imported MaxQuant MS data to ', mqrdata_filepath, '...')
save(data_info, msruns.df, fasta.dfs,
     msdata.wide, mqevidence, instr_calib_pepmodstate, instr_calib_protgroup,
     file = mqrdata_filepath, compress = "xz")
load(mqrdata_filepath)

msdata_colgroups <- attr(msdata.wide, "column_groups")

msdata_full <- list()

msdata_full$mschannels <- filter(msruns.df, is_used) %>%
    dplyr::inner_join(dplyr::filter(mqevidence$mschannels, quant_type == "SILAC") %>%
                      dplyr::select(msrun, mstag, mschannel)) %>%
    dplyr::semi_join(tibble(.name = msdata_colgroups$LFQ) %>%
                       extract(.name, c("mstag", "msrun"), remove=FALSE,
                               str_c("^", data_info$quant_col_prefix, "\\.(\\S+)\\s(\\S+)"))) %>%
    dplyr::select(-is_ptm, -channels, -msfraction) %>%
    dplyr::arrange(treatment, timepoint, replicate, mstag) %>%
    dplyr::mutate(mschannel = factor(mschannel, levels=mschannel),
                  msrun = factor(msrun, levels=unique(msrun)),
                  sample = factor(sample, levels=unique(sample)),
                  condition = factor(condition, levels=unique(condition)),
                  tagcondition = str_c(condition, "_", mstag),
                  tagcondition = factor(tagcondition, levels=unique(tagcondition)))

msdata_full <- append_protgroups_info(msdata_full, msdata.wide,
                                      proteins_info = dplyr::bind_rows(
                                          dplyr::mutate(fasta.dfs$hbv, is_viral = TRUE),
                                          dplyr::mutate(fasta.dfs$human, is_viral = FALSE)),
                                      import_columns = c("is_viral", "organism"))

msdata_full$proteins <- mutate(msdata_full$proteins,
                               protein_ac_noiso = str_remove(protein_ac, "-\\d+$"))

pacnoiso_stats.df <- left_join(msdata_full$protein2protgroup,
                               select(msdata_full$proteins, protein_ac, protein_ac_noiso)) %>%
    #filter(is_majority) %>%
    group_by(protein_ac_noiso) %>%
    summarise(n_noiso_pgs_have_razor = sum(npeptides_razor > 0),
              n_noiso_pgs = n_distinct(protgroup_id)) %>%
    ungroup()

msdata_full$proteins <- left_join(msdata_full$proteins, pacnoiso_stats.df)

pg_razor_stats.df <- filter(msdata_full$protein2protgroup, is_majority) %>%
    left_join(select(msdata_full$proteins, protein_ac, n_noiso_pgs, n_noiso_pgs_have_razor)) %>%
    group_by(protgroup_id) %>%
    summarise(nproteins_have_razor = sum(npeptides_razor > 0),
              nprotgroups_sharing_proteins = max(n_noiso_pgs)) %>%
    dplyr::ungroup()

strlist_label <- function(strs) {
  str_c(strs[[1]], if_else(n_distinct(strs) > 1, '...', ''))
}
strlist_label2 <- function(strs, delim=fixed(';')) {
  sapply(str_split(strs, delim), strlist_label)
}

msdata_full$protgroups <- left_join(msdata_full$protgroups, pg_razor_stats.df)
msdata_full$protgroups <- dplyr::mutate(msdata_full$protgroups,
                                        gene_label = strlist_label2(gene_names),
                                        protac_label = strlist_label2(protein_acs),
                                        protein_label = strlist_label2(protein_names),
                                        protgroup_label = case_when(is_viral ~ protein_label,
                                                                    !is.na(gene_label) ~ gene_label,
                                                                    !is.na(protac_label) ~ protac_label,
                                                                    TRUE ~ str_c('#', protgroup_id)))

msdata_full$peptides <- mqevidence$peptides %>%
    dplyr::left_join(select(msdata_full$proteins, lead_razor_protein_ac = protein_ac, is_viral)) %>%
    mutate(is_viral=replace_na(is_viral, FALSE))
msdata_full$pepmods <- mqevidence$pepmods %>%
    dplyr::left_join(select(msdata_full$peptides, peptide_id, is_viral)) %>%
    dplyr::mutate(is_used = TRUE,
                  pepmod_rank = case_when(!is_used ~ -1L,
                                          TRUE ~ 1L))
msdata_full$pepmodstates <- mqevidence$pepmodstates

# redefine protein groups (protregroups)
pepmods.df <- dplyr::select(msdata_full$pepmods, pepmod_id, protgroup_ids, protein_acs, lead_protein_acs, seq, modifs, charges, is_reverse, pepmod_rank)
proteins.df <- msdata_full$proteins
save(file = file.path(msdata_path, str_c(project_id, "_", ms_folder, '_', data_version, "_pepmods.RData")),
     pepmods.df, proteins.df)
# .. run protregroup.jl
msdata_full$protregroups <- read_tsv(file.path(data_path, ms_folder,
                                               str_c(project_id, "_", ms_folder, '_', data_version, "_protregroups_acs.txt")),
                                     col_types = list(protregroup_id = "i")) %>%
    dplyr::mutate(is_contaminant = str_detect(majority_protein_acs, "(;|^)CON__"),
                  is_reverse = str_detect(majority_protein_acs, "(;|^)REV__"))

msdata_full$protein2protregroup <- dplyr::select(msdata_full$protregroups, protregroup_id, protein_ac=majority_protein_acs) %>%
    separate_rows(protein_ac, sep=fixed(";"), convert=TRUE) %>%
    dplyr::mutate(is_majority = TRUE) %>%
    dplyr::group_by(protregroup_id) %>%
    dplyr::mutate(protein_ac_rank = row_number()) %>%
    dplyr::ungroup()

msdata_full$protregroup2pepmod <- bind_rows(
    select(msdata_full$protregroups, protregroup_id, pepmod_id=spec_pepmod_ids) %>%
        separate_rows(pepmod_id, sep=fixed(";"), convert=TRUE) %>%
        mutate(is_specific = TRUE),
    select(msdata_full$protregroups, protregroup_id, pepmod_id=pepmod_ids) %>%
        separate_rows(pepmod_id, sep=fixed(";"), convert=TRUE) %>%
        mutate(is_specific = FALSE)) %>%
    dplyr::group_by(protregroup_id, pepmod_id) %>%
    dplyr::summarise(is_specific = any(is_specific)) %>%
    dplyr::ungroup()
msdata_full$protregroup2pepmodstate <- dplyr::select(inner_join(msdata_full$protregroup2pepmod, msdata_full$pepmodstates),
                                                protregroup_id, pepmodstate_id, is_specific)

msdata_full$protregroups <- dplyr::inner_join(msdata_full$protregroups,
                                              dplyr::inner_join(msdata_full$protregroups, msdata_full$protein2protregroup) %>%
                                                  dplyr::inner_join(msdata_full$proteins) %>%
                                                  dplyr::arrange(protregroup_id, protein_ac_rank) %>%
                                                  dplyr::group_by(protregroup_id) %>%
                                                  dplyr::summarise(gene_names = str_c(gene_name, collapse=';'),
                                                                   gene_label = strlist_label(gene_name),
                                                                   protein_names = str_c(protein_name, collapse=';'),
                                                                   protein_label = strlist_label(protein_name),
                                                                   protein_descriptions = str_c(protein_description, collapse=';'),
                                                                   protein_description = strlist_label(protein_description),
                                                                   organism = str_c(unique(organism), collapse=';'),
                                                                   is_viral = any(coalesce(is_viral, FALSE))) %>%
                                                  dplyr::ungroup()) %>%
    dplyr::mutate(protac_label = sapply(str_split(majority_protein_acs, fixed(";")), strlist_label),
                  protregroup_label = case_when(is_viral ~ protein_label,
                                                !is.na(gene_label) ~ gene_label,
                                                !is.na(protac_label) ~ protac_label,
                                                TRUE ~ str_c('#', protregroup_id))) %>%
    dplyr::left_join(dplyr::inner_join(msdata_full$protregroup2pepmod, msdata_full$pepmods) %>%
                         dplyr::group_by(protregroup_id) %>%
                         dplyr::summarise(npeptides = n_distinct(peptide_id),
                                          npepmods = n_distinct(pepmod_id),
                                          npeptides_unique = n_distinct(peptide_id[is_specific]),
                                          npepmods_unique = n_distinct(pepmod_id[is_specific])) %>%
                         dplyr::ungroup() %>%
                         dplyr::mutate(npeptides_unique_razor = npeptides_unique,
                                       npeptides_razor = 0L,
                                       npepmods_unique_razor = npepmods_unique,
                                       npepmods_razor = 0L))
# prepare protgroup intensities (wider format: all mstags in one row)
intensity_prespec_df <- tibble(.name = msdata_colgroups$LFQ) %>%
    extract(.name, c("mstag", "msrun"), remove=FALSE,
            str_c("^", data_info$quant_col_prefix, "\\.(\\S+)\\s(\\S+)")) %>%
    mutate(.value = str_c("LFQ_Intensity.", mstag)) %>%
    dplyr::inner_join(dplyr::select(msdata_full$mschannels, msrun, mschannel, mstag, raw_file))

# prepare protgroup intensities (longer format: each mstag on ots own row)
protgroup_tagintensities_all.df <- tidyr::pivot_longer_spec(
    dplyr::select(msdata.wide, protgroup_id, !!msdata_colgroups$LFQ),
    mutate(intensity_prespec_df, .value = "intensity"))

protgroup_intensities_all.df <- tidyr::pivot_longer(
  dplyr::select(msdata.wide, protgroup_id, all_of(msdata_colgroups$LFQ)),
  msdata_colgroups$LFQ, names_to=c(".value", "msrun"),
  names_pattern=str_c("^", data_info$quant_col_prefix, "\\.(\\S+)\\s(\\S+)")) %>%
  dplyr::rename(intensity.L = L, intensity.H = H) %>%
  dplyr::mutate(log2_intensity.L = log2(intensity.L),
                log2_intensity.H = log2(intensity.H),
                log2_ratio.HL = log2_intensity.H - log2_intensity.L)

p <- ggplot(filter(protgroup_intensities_all.df, !is.na(log2_ratio)) %>%
         inner_join(msruns.df), aes(x=log2_ratio.HL, y=..count.., color=timepoint_num, group=msrun)) +
  geom_density(alpha=0.5) +
  scale_color_viridis_c(breaks = unique(msruns.df$timepoint_num)) +
  facet_grid(treatment ~ .) +
  theme_bw_ast(base_size=10, base_family = "")
ggsave(p, filename = file.path(analysis_path, "plots", str_c(ms_folder, '_', data_version),
                               paste0(project_id, "_", ms_folder, "_HL_log2ratio_", data_version, ".pdf")),
       width = 8, height = 8, device = cairo_pdf)

protgroup_stats.df <- dplyr::filter(protgroup_intensities_all.df, !is.na(log2_ratio.HL)) %>%
  dplyr::semi_join(dplyr::filter(msruns.df, timepoint_num > 0)) %>%
  dplyr::group_by(protgroup_id) %>%
  dplyr::summarise(n_ratios = n_distinct(msrun),
                   log2_ratio_min = min(log2_ratio.HL),
                   log2_ratio_max = max(log2_ratio.HL),
                   log2_ratio_median = median(log2_ratio.HL),
                   log2_ratio_25 = quantile(log2_ratio.HL, 0.25),
                   log2_ratio_75 = quantile(log2_ratio.HL, 0.75),
                   .groups = "drop") %>%
  dplyr::inner_join(dplyr::select(msdata_full$protgroups, protgroup_id, protgroup_label, is_contaminant, is_reverse, is_viral)) %>%
  dplyr::mutate(is_stable = !is_contaminant & !is_reverse & !is_viral & (log2_ratio_75 < 0) &
                            (log2_ratio_75 - log2_ratio_25 < 3) & (n_ratios > nrow(msruns.df) * 0.75))
  
protgroup_idents_all.df <- tidyr::pivot_longer(
    dplyr::select(msdata.wide, protgroup_id, !!msdata_colgroups$ident_type, !!msdata_colgroups$ms2_count),
    cols = c(msdata_colgroups$ident_type, msdata_colgroups$ms2_count),
    names_to = c(".value", "msrun"),
    names_pattern = "(ident_type|ms2_count)\\.(.*)") %>%
    dplyr::left_join(select(msruns.df, msrun) %>% distinct())
if (!has_name(protgroup_idents_all.df, "ident_type")) {
    # when matching disabled no ident_type is written
    protgroup_idents_all.df <- mutate(protgroup_idents_all.df,
                                      ident_type = factor(if_else(ms2_count > 0, "By MS/MS", NA_character_),
                                                          levels=c("By matching", "By MS/MS")))
}
msdata_full$protgroup_tagintensities <- protgroup_tagintensities_all.df %>%
    dplyr::semi_join(dplyr::select(msdata_full$mschannels, msrun, raw_file)) %>%
    dplyr::filter(!is.na(intensity)) %>%
    dplyr::select(-raw_file)
msdata_full$protgroup_idents <- protgroup_idents_all.df %>%
    dplyr::semi_join(dplyr::select(msdata_full$mschannels, msrun)) %>%
    dplyr::filter(!is.na(ident_type))

msdata_full$pepmodstate_intensities <- mqevidence$pepmodstate_intensities %>%
    dplyr::select(pepmod_id, pepmodstate_id, msrun, ident_type, starts_with("intensity.")) %>%
    dplyr::mutate(is_idented = ident_type %in% c("ISO-MSMS", "MULTI-MSMS", "MSMS"),
                  log_ratio.HL = log(intensity.H) - log(intensity.L))
msdata_full$pepmodstate_tagintensities <- tidyr::pivot_longer(msdata_full$pepmodstate_intensities,
  starts_with("intensity"), names_prefix = "intensity.", names_to = "mstag", values_to = "intensity") %>%
  dplyr::mutate(mstag = factor(mstag, levels=levels(msdata_full$mschannels$mstag))) %>%
  dplyr::left_join(dplyr::select(msdata_full$mschannels, msrun, mstag, mschannel))

msdata_full$protregroup_idents <- dplyr::inner_join(msdata_full$protregroup2pepmod, msdata_full$pepmodstate_tagintensities) %>%
    dplyr::group_by(msrun, protregroup_id) %>%
    dplyr::summarise(npepmods_quanted = sum(!is.na(intensity)),
                     nspecpepmods_quanted = sum(!is.na(intensity) & is_specific),
                     npepmods_idented = sum(is_idented),
                     nspecpepmods_idented = sum(is_idented & is_specific)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ident_type = factor(if_else(nspecpepmods_idented > 0, "By MS/MS", "By matching")))

# condition = treatment X timepoint
tagconditions.df <- dplyr::select(msdata_full$mschannels, tagcondition, condition, treatment, timepoint, timepoint_num, mstag) %>%
  mutate(mstag = relevel(mstag, "L")) %>%
  dplyr::distinct()
tagconditions.df <- dplyr::left_join(tagconditions.df,
    dplyr::select(tagconditions.df, treatment, mstag, timepoint0=timepoint_num)) %>%
    dplyr::mutate(is_modeled = (timepoint0 >= 0 & timepoint_num >= timepoint0) |
                               (timepoint0 == 1 & timepoint_num <= timepoint0),
                  timediff = str_c(if_else(timepoint_num >= timepoint0, "after", "before"), timepoint0, "h")) %>%
    tidyr::pivot_wider(all_of(c("tagcondition", "condition", "treatment", "timepoint", "timepoint_num", "mstag")),
                       names_from = timediff, values_from = is_modeled, values_fill=FALSE) %>%
    dplyr::arrange(mstag, treatment, timepoint_num)

msdata_full$mschannel_stats <- mschannel_statistics(msdata_full) %>%
  dplyr::mutate(na_ratio = n_missing / n)
  
msdata <- msdata_full[c('protgroup_tagintensities', 'pepmodstate_tagintensities', 'pepmods', 'pepmodstates',
                        'mschannels', 'protgroups', 'protregroups', 'protein2protregroup',
                        'protregroup2pepmod', 'protregroup2pepmodstate')]
#msdata$protregroup2pepmodstate <- dplyr::select(inner_join(msdata$protregroup2pepmod, msdata$pepmodstates),
#                                                protregroup_id, pepmodstate_id, is_specific)
#msdata$protregroups <- semi_join(msdata$protregroups, msdata$protregroup2pepmodstate)
#msdata$protgroup_tagintensities <- semi_join(msdata$protgroup_tagintensities, msdata$mschannel)
#msdata$pepmodstate_tagintensities <- semi_join(semi_join(msdata$pepmodstate_tagintensities, msdata$mschannels),
#                                               msdata$pepmodstates)

# setup experimental design matrices
tagconditionXeffect_orig.mtx <- model.matrix(
  formula(str_c("~ 1 + mstag:(", str_c(unique(str_subset(colnames(tagconditions.df), "^(after|before)\\d+h")), collapse =" + "), ")*treatment",
                " + mstag - treatment")),
  tagconditions.df
)

tagconditionXeffect.mtx <- tagconditionXeffect_orig.mtx[, str_detect(colnames(tagconditionXeffect_orig.mtx),
                                                               "^\\(Intercept\\)$|^mstagL:after0hTRUE$|^mstagL:before|^mstagH:after0hTRUE|^mstagH:after1hTRUE$", negate=TRUE) &
                                                        colnames(tagconditionXeffect_orig.mtx) != "mstagL" &
                                                        apply(abs(tagconditionXeffect_orig.mtx), 2, sum) != 0] # remove intercept and afterT0
dimnames(tagconditionXeffect.mtx) <- list(tagcondition = as.character(tagconditions.df$tagcondition),
                                          effect = colnames(tagconditionXeffect.mtx))

pheatmap(tagconditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, 
         filename = file.path(analysis_path, "plots", str_c(ms_folder, "_", fit_version),
                              paste0(project_id, "_exp_design_", ms_folder, "_", fit_version, ".pdf")),
         width = 10, height = 10)

effects.df <- tibble(effect=colnames(tagconditionXeffect.mtx)) %>%
  mutate(mstag = effect_factor(effect, 'mstag',
                               levels(tagconditions.df$mstag), NA),
         mstag = if_else(mstag == "H", mstag, factor(NA_character_, levels=levels(tagconditions.df$mstag))),
         treatment = effect_factor(effect, 'treatment',
                                   levels(tagconditions.df$treatment), NA),
         timepoint = case_when(str_detect(effect, "after") ~ effect_factor(effect, 'after', str_c(unique(tagconditions.df$timepoint_num), "hTRUE"), NA) %>%
                                   as.character() %>% str_remove_all("^after|hTRUE$"),
                               str_detect(effect, "before") ~ "0",
                               TRUE ~ NA_character_) %>%
           factor(levels = levels(tagconditions.df$timepoint)))

effects.df <- mutate(effects.df,
    effect_type = case_when(!is.na(mstag) & !is.na(timepoint) & !is.na(treatment) ~ "mstagXtreatmentXtimepoint",
                            !is.na(timepoint) & !is.na(treatment) ~ "treatmentXtimepoint",
                            !is.na(timepoint) & is.na(treatment) ~ "timepoint",
                            !is.na(treatment) ~ "treatment",
                            !is.na(mstag) ~ "mstag",
                            TRUE ~ NA_character_),
   effect_label = case_when(effect_type == "mstagXtreatmentXtimepoint" ~ str_c(mstag, "+", treatment, "+", timepoint, "h"),
                            effect_type == "treatmentXtimepoint" ~ str_c(treatment, "+", timepoint, "h"),
                            effect_type == "treatment" ~ as.character(treatment),
                            effect_type == "timepoint" ~ str_c(timepoint, "h"),
                            effect_type == "mstag" ~ as.character(mstag),
                            TRUE ~ NA_character_),
   prior_mean = 0,
   prior_tau = case_when(effect_type == "mstagXtreatmentXtimepoint" ~ 0.01,
                         effect_type == "treatmentXtimepoint" ~ 0.01,
                         effect_type == "timepoint" ~ 0.05,
                         effect_type == "mstag" ~ 5,
                         TRUE ~ NA_real_),
  prior_df2 = case_when(effect_type == "mstagXtreatmentXtimepoint" ~ 8.0,
                        effect_type == "treatmentXtimepoint" ~ 8.0,
                        effect_type == "timepoint" ~ 2.0,
                        effect_type == "timepoint_add" ~ 2.0,
                        TRUE ~ 1.0),
  is_positive = FALSE)

tagconditionXeffect.df <- conditionXeffect_frame(tagconditionXeffect.mtx, effects.df)

#TODO 0h?
#compound_metaconditions <- c("allTreatments_initialTimepoint")
compound_metaconditions <- c()
all_metaconditions <- c(levels(tagconditions.df$tagcondition), compound_metaconditions)
tagconditionXmetacondition.mtx <- false_matrix(tagcondition = levels(tagconditions.df$tagcondition),
                                               metacondition = all_metaconditions)
for (cname in levels(tagconditions.df$tagcondition)) {
  tagconditionXmetacondition.mtx[cname, cname] <- TRUE
}
#conditionXmetacondition.mtx[filter(conditions.df, timepoint_num == min(timepoint_num))$condition, "allTreatments_initialTimepoint"] <- TRUE
pheatmap(ifelse(tagconditionXmetacondition.mtx, 1.0, 0.0), cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(analysis_path, "plots", str_c(ms_folder, "_", fit_version),
                              paste0(project_id, "_tagmetaconditions_", ms_folder, "_", fit_version, ".pdf")),
         width = 8, height = 6)

tagconditionXmetacondition.df <- as_tibble(as.table(tagconditionXmetacondition.mtx)) %>%
  dplyr::filter(n != 0) %>% dplyr::select(-n)

all_contrasts <- character()

contrasts.df <- bind_rows(
  # all treatment pairs at each timepoint
  left_join(dplyr::mutate(dplyr::select(tagconditions.df, mstag, metacondition_lhs = tagcondition, treatment_lhs = treatment, timepoint_lhs = timepoint)),
            dplyr::mutate(dplyr::select(tagconditions.df, mstag, metacondition_rhs = tagcondition, treatment_rhs = treatment, timepoint_lhs = timepoint))) %>%
  filter(as.integer(treatment_lhs) > as.integer(treatment_rhs)) %>%
  mutate(timepoint_rhs = timepoint_lhs,
         contrast_kind = "treatment_vs_treatment"),
  # all timepoints of the same treatment
  left_join(dplyr::mutate(dplyr::select(tagconditions.df, mstag, metacondition_lhs = tagcondition, treatment_lhs = treatment, timepoint_lhs = timepoint)),
            dplyr::mutate(dplyr::select(tagconditions.df, mstag, metacondition_rhs = tagcondition, treatment_lhs = treatment, timepoint_rhs = timepoint))) %>%
  filter(as.integer(timepoint_lhs) > as.integer(timepoint_rhs)) %>%
  mutate(treatment_rhs = treatment_lhs,
         contrast_kind = "timepoint_vs_timepoint")
) %>%
  mutate(contrast = str_c(mstag, ":", treatment_lhs, "@", timepoint_lhs, "h_vs_", treatment_rhs, "@", timepoint_rhs, "h"),
         contrast_type = "comparison",
         offset = 0)
all_contrasts = contrasts.df$contrast

contrastXmetacondition.mtx <- zero_matrix(contrast = all_contrasts, metacondition = all_metaconditions)
for (i in 1:nrow(contrasts.df)) {
  contrastXmetacondition.mtx[contrasts.df$contrast[[i]],
                             c(contrasts.df$metacondition_lhs[[i]],
                               contrasts.df$metacondition_rhs[[i]])] <- c(1,-1)
}

contrastXmetacondition.df <- as_tibble(as.table(contrastXmetacondition.mtx)) %>% dplyr::filter(n != 0) %>%
  dplyr::rename(weight = n) %>%
  dplyr::mutate(contrast_type = 'comparison',
                condition_role = "signal")

contrastXtagcondition.df <- as_tibble(as.table(tagconditionXmetacondition.mtx)) %>% dplyr::filter(n != 0) %>%
  dplyr::select(-n) %>%
  dplyr::inner_join(contrastXmetacondition.df) %>%
  dplyr::arrange(contrast, contrast_type, metacondition, tagcondition)

pheatmap(contrastXmetacondition.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(analysis_path, "plots", str_c(ms_folder, "_", fit_version),
                              paste0(project_id, "_exp_design_contrasts_", ms_folder, "_", fit_version, ".pdf")),
         width = 8, height = 24)

##########
## normalization
msdata4norm.df <- msdata$pepmodstate_tagintensities %>% ungroup() %>%
  dplyr::filter(is_idented & !is.na(intensity)) %>% 
  dplyr::inner_join(dplyr::filter(msdata$protregroup2pepmod, is_specific)) %>%
  dplyr::semi_join(dplyr::filter(msdata$protregroups, !is_reverse & !is_contaminant & !is_viral)) %>%
  dplyr::semi_join(dplyr::inner_join(dplyr::filter(protgroup_stats.df, is_stable),
                                     mqevidence$protgroups2protgroup) %>%
                   dplyr::inner_join(dplyr::select(msdata_full$pepmods, protgroup_ids, pepmod_id)) %>%
                   dplyr::select(pepmod_id) %>% dplyr::distinct()) %>%
  dplyr::select(pepmodstate_id) %>% dplyr::distinct() %>%
  dplyr::inner_join(msdata$pepmodstate_tagintensities) %>%
  dplyr::inner_join(dplyr::select(msdata$mschannels, mschannel, msrun, mstag, tagcondition, treatment, timepoint)) %>%
  dplyr::filter(mstag == "L")

options(mc.cores=8)

msruns_hnorm <- multilevel_normalize_experiments(instr_calib_pepmodstate, msdata$mschannels,
                                                 msdata4norm.df,
                                                 quant_col = "intensity", obj_col = "pepmodstate_id", mschan_col = "mschannel",
                                                 mcmc.iter = 2000L,
                                                 #mcmc.chains = 6,
                                                 verbose=TRUE,
                                                 norm_levels = list(msrun = list(cond_col = "msrun", max_objs=1000L, missing_exp.ratio=0.1),
                                                                    condition = list(cond_col="tagcondition", max_objs=1000L, missing_exp.ratio=0.2),
                                                                    timepoint = list(cond_col="timepoint", max_objs=1000L, missing_exp.ratio=0.3)
                                                 ))
save(msruns_hnorm, file=file.path(scratch_path, str_c(project_id, "_", data_version, "_msruns_hnorm.RData")))

total_msrun_shifts.df <- dplyr::filter(msruns_hnorm$mschannel_shifts, str_detect(mschannel, "_L$")) %>%
  dplyr::select(-mschannel) %>%
  dplyr::rename(total_msrun_shift = total_mschannel_shift)
#total_msrun_shifts.df <- transmute(msdata$msruns, msrun, total_msrun_shift=0.0) # no normalization

## apply normalization
msdata$protgroup_tagintensities <- dplyr::left_join(msdata$protgroup_tagintensities,
                                                    dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift)) 
#dplyr::select(-total_msrun_shift)

# calculate the general trends for H/L ratios and the corresponding model effect priors
pepmodstate_ratio_stats.df <- dplyr::filter(msdata_full$pepmodstate_intensities, !is.na(log_ratio.HL) & intensity.L > 1E6) %>%
  dplyr::semi_join(dplyr::filter(msdata_full$pepmods, !is_contaminant & !is_reverse & !is_viral)) %>%
  dplyr::inner_join(dplyr::select(msruns.df, msrun, condition, treatment, timepoint)) %>%
  dplyr::group_by(treatment, condition, timepoint) %>%
  dplyr::summarise(npepmodstates = n_distinct(pepmodstate_id),
                   log_ratio_HL_sd = sd(log_ratio.HL),
                   log_ratio_HL_mean = mean(log_ratio.HL),
                   log_ratio_HL_median = median(log_ratio.HL),
                   .groups = "drop") %>%
  dplyr::mutate(log_ratio_HL_mock_1d_median = log_ratio_HL_median[treatment == "mock" & timepoint == 1]) %>%
  dplyr::group_by(timepoint) %>%
  dplyr::mutate(log_ratio_HL_mock_timepoint_delta = log_ratio_HL_median[treatment == "mock"] - log_ratio_HL_mock_1d_median,
                log_ratio_HL_treatmentXtimepoint_delta = log_ratio_HL_median - log_ratio_HL_mock_timepoint_delta - log_ratio_HL_mock_1d_median) %>%
  dplyr::ungroup()

HtagconditionXeffect.mtx <- tagconditionXeffect.mtx[str_detect(rownames(tagconditionXeffect.mtx), "_H"),
                                                    str_detect(colnames(tagconditionXeffect.mtx), "mstagH")]
Htagcondition_log_ratio_HL_median <- set_names(pepmodstate_ratio_stats.df$log_ratio_HL_median,
          str_c(pepmodstate_ratio_stats.df$condition, "_H"))

Htageffect_priors <- solve(HtagconditionXeffect.mtx, Htagcondition_log_ratio_HL_median[rownames(HtagconditionXeffect.mtx)])

effects.df <- dplyr::mutate(effects.df, new_prior_mean = Htageffect_priors[effect],
                            prior_mean = case_when(effect == "mstagH" ~ new_prior_mean,
                                                   !is.na(new_prior_mean) & effect_type == "timepoint" ~ pmin(new_prior_mean, 0),
                                                   !is.na(new_prior_mean) & effect_type == "mstagXtreatmentXtimepoint" & abs(new_prior_mean) > 0.05 ~ new_prior_mean,
                                                   !is.na(new_prior_mean) ~ 0.0,
                                                   TRUE ~ 0.0),
                            new_prior_mean = NULL)

global_protgroup_labu_shift <- 0.95*median(log(dplyr::filter(msdata$mschannels, mstag == "L") %>%
                                               dplyr::select(mschannel) %>% dplyr::distinct() %>%
                                               dplyr::inner_join(msdata$protgroup_tagintensities) %>% .$intensity), na.rm=TRUE)
global_pepmodstate_labu_shift <- 0.95*median(log(dplyr::filter(msdata$mschannels, mstag == "L") %>%
                                        dplyr::select(mschannel) %>% dplyr::distinct() %>%
                                        dplyr::inner_join(msdata$pepmodstate_tagintensities) %>% .$intensity), na.rm=TRUE)

pepmodstate_labu_min <- inner_join(msdata_full$pepmodstate_intensities, total_msrun_shifts.df) %>%
  mutate(intensity_norm.L = intensity.L * exp(-total_msrun_shift)) %>%
  .$intensity_norm.L %>% log() %>%
  quantile(0.001, na.rm=TRUE) - global_pepmodstate_labu_shift - 5

############################
# batch effects
# no batch effects so far
msrunXbatchEffect.mtx <- zero_matrix(msrun = unique(msdata$mschannels$msrun),
                                     batch_effect = c())

batch_effects.df <- tibble(batch_effect = character(0),
                           is_positive = logical(0),
                           prior_mean = double(0))

# no subbatch effects so far
msrunXsubbatchEffect.mtx <- zero_matrix(msrun = rownames(msrunXbatchEffect.mtx),
                                        subbatch_effect = c())

subbatch_effects.df <- tibble(subbatch_effect=character(0),
                              is_positive=logical(0))

rmsglmdata_filepath <- file.path(scratch_path, str_c(project_id, '_msglm_data_', ms_folder, '_', data_version, '.RData'))
message('Saving MS data for MSGLM to ', rmsglmdata_filepath, '...')
save(data_info, msdata,
     tagconditions.df, effects.df, contrasts.df,
     tagconditionXeffect.mtx, tagconditionXeffect.df,
     tagconditionXmetacondition.mtx, tagconditionXmetacondition.df,
     contrastXmetacondition.mtx, contrastXmetacondition.df, contrastXtagcondition.df,
     instr_calib_protgroup, instr_calib_pepmodstate,
     global_protgroup_labu_shift, global_pepmodstate_labu_shift, pepmodstate_labu_min,
     total_msrun_shifts.df, msruns_hnorm, 
     batch_effects.df, msrunXbatchEffect.mtx,
     subbatch_effects.df, msrunXsubbatchEffect.mtx,
     file = rmsglmdata_filepath)

rfulldata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_', ms_folder, '_', data_version, '.RData'))
message('Saving full MS data to ', rfulldata_filepath, '...')
save(data_info, msdata_full,
     file = rfulldata_filepath)

message('Done.')

set.seed(1232)
msdata_full$protgroup_tagintensities_all <- tidyr::expand(msdata_full$protgroup_tagintensities,
                                                          protgroup_id, mschannel) %>%
  left_join(dplyr::select(msdata_full$mschannels, mschannel, msrun, mstag)) %>%
  left_join(dplyr::select(msdata_full$protgroup_tagintensities, -any_of("total_msrun_shift"))) %>%
  impute_intensities(msdata_full$mschannel_stats) %>%
  #dplyr::mutate(total_msrun_shift=0) %>%
  dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift),
                intensity_imputed_norm = intensity_imputed*exp(-total_msrun_shift),
                log2_intensity_imputed_norm = log2(intensity_imputed_norm),
                is_imputed = is.na(intensity)) %>%
  dplyr::arrange(mstag, mschannel, protgroup_id)

protgroup_intensities4pca.df <- msdata_full$protgroup_tagintensities_all %>%
  dplyr::semi_join(dplyr::filter(msdata_full$protgroups, !is_reverse & !is_contaminant)) %>%
  dplyr::semi_join(dplyr::filter(msdata_full$mschannels, is_used)) %>%
  dplyr::arrange(mstag, mschannel, protgroup_id)

protgroup_intensities.mtx <- matrix(log2(protgroup_intensities4pca.df$intensity_norm),
                                    nrow = n_distinct(protgroup_intensities4pca.df$protgroup_id),
                                    dimnames = list(protgroup = unique(protgroup_intensities4pca.df$protgroup_id),
                                                    mschannel = unique(protgroup_intensities4pca.df$mschannel)))
protgroup_intensities_imp.mtx <- matrix(log2(protgroup_intensities4pca.df$intensity_imputed_norm),
                                        nrow = n_distinct(protgroup_intensities4pca.df$protgroup_id),
                                        dimnames = list(protgroup = unique(protgroup_intensities4pca.df$protgroup_id),
                                                        mschannel = unique(protgroup_intensities4pca.df$mschannel)))
require(Matrix)
msrun_intensities_pca <- stats::prcomp(protgroup_intensities_imp.mtx, scale. = TRUE)
msrun_intensities_pca.df <- as_tibble(msrun_intensities_pca$rotation,
                                      rownames="mschannel") %>%
  dplyr::inner_join(msdata$mschannels)

require(ggrepel)
p <- ggplot(filter(msrun_intensities_pca.df, mstag=="L" & timepoint_num > 0) %>% arrange(treatment, timepoint),
       aes(x=PC1, y=PC2, color=treatment)) +
    geom_point() +
    geom_text_repel(aes(label=mschannel), size=3, vjust=-1.1) +
    theme_bw_ast(base_family = "", base_size = 10) #+
ggsave(p, filename = file.path(analysis_path, "plots", str_c(ms_folder, "_", data_version),
                               str_c(project_id, "_", ms_folder, "_msruns_pca_norm_", data_version, ".pdf")),
       width = 12, height = 14)

require(pheatmap)
require(matrixStats)
mschannels_ordered <- unique(as.character(protgroup_intensities4pca.df$mschannel))
mstagL_mask <- str_detect(mschannels_ordered, "_L$")
protgroup_hclu = hclust(dist(protgroup_intensities_imp.mtx / rowMedians(protgroup_intensities_imp.mtx[, mstagL_mask])))
shown_protgroup_intensities.mtx <- protgroup_intensities.mtx[, mschannels_ordered]/rowMedians(protgroup_intensities_imp.mtx[, mstagL_mask])
pheatmap(shown_protgroup_intensities.mtx, cluster_cols=FALSE, cluster_rows=protgroup_hclu,
         file = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_", data_version, "_heatmap_intensity.pdf")), width=16, height=80)
pheatmap(log2(pmax(pmin(shown_protgroup_intensities.mtx,
                        quantile(shown_protgroup_intensities.mtx, 0.99, na.rm=TRUE)),
                   quantile(shown_protgroup_intensities.mtx, 0.01, na.rm=TRUE))),
         cluster_cols=FALSE, cluster_rows=protgroup_hclu,
         file = file.path(analysis_path, "plots", str_c(ms_folder, "_", fit_version),
                          paste0(project_id, "_", ms_folder, "_", data_version, "_heatmap_intensity.png")), width=35, height=60)

viral_mask <- unique(protgroup_intensities4pca.df$protgroup_id) %in% filter(msdata_full$protgroups, is_viral)$protgroup_id
viral_intensities.mtx <- protgroup_intensities.mtx[viral_mask, msruns_ordered]
rownames(bait_intensities.mtx) <- 
  tibble(protgroup_id = unique(protgroup_intensities4pca.df$protgroup_id)[viral_mask]) %>%
  inner_join(msdata_full$protgroups) %>% inner_join(msdata_full$protein2protgroup) %>%
  inner_join(mutate(msdata_full$proteins, full_code = str_c(str_replace(organism, "CoV-2", "CoV2"), '-', gene_name))) %>%
  group_by(protgroup_id) %>%
  summarise(full_code = str_c(full_code, collapse=';')) %>% ungroup() %>% .$full_code
pheatmap(bait_intensities.mtx, cluster_cols=FALSE, cluster_rows=FALSE,
         file = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_", data_version, "_viral_heatmap_intensity.pdf")), width=16, height=8)

require(uwot)
require(matrixStats)
protgroup_intensities_imp_pg.mtx <- protgroup_intensities_imp.mtx - rowMedians(protgroup_intensities_imp.mtx)
pg_umap2d = umap(protgroup_intensities_imp_norm.mtx, n_components=2,
                 n_neighbors=20, init="laplacian",
                 min_dist=0.2, metric = "euclidean")

objects_umap.df <- as_tibble(pg_umap2d, .name_repair="minimal") %>%
  set_names(c("x_2d", "y_2d")) %>%
  mutate(object_id = unique(protgroup_intensities4pca.df$protgroup_id)) %>%
  #dplyr::left_join(
  #  as_tibble(objects.umap3d, .name_repair="minimal") %>%
  #    set_names(c("x_3d", "y_3d", "z_3d")) %>%
  #    mutate(object_id = unique(object_intensities4pca.df$object_id))    
  #) %>%
  dplyr::inner_join(select(msdata$protgroups, object_id=protgroup_id, object_label=protgroup_label,
                           majority_protein_acs, gene_names, protein_descriptions,
                           is_viral, is_expected_stable, is_contaminant, is_reverse)) %>%
  mutate(category = case_when(replace_na(is_viral, FALSE) ~ "viral",
                              is_reverse ~ "reverse",
                              is_contaminant ~ "contaminant",
                              is.na(gene_names) ~ "NA",
                              is_expected_stable ~ "stable",
                              str_detect(gene_names, "(^|;)(?:HIST[123])?H[123][ABCDEF]") ~ "histone",
                              str_detect(gene_names, "(^|;)PPP[123456]") ~ "phosphatase",
                              str_detect(gene_names, "(^|;)RP[SL]\\d+") ~ "ribosome",
                              str_detect(gene_names, "(^|;)MRP[SL]\\d+") ~ "mito ribosome",
                              str_detect(gene_names, "(^|;)PSM[ABCD]\\d+") ~ "proteasome",
                              str_detect(gene_names, "(^|;)(MAVS|OAS\\d+|IFIT\\d+|PARP\\d+|DDX58|IFI\\d+|ISG\\d+|MX\\d+|TGFB\\d+|UNC93B1|ZC3HAV1)(;|$)") ~ "ISG",
                              TRUE ~ "default"
  ))# %>%
  #dplyr::left_join(filter(protregroup_effects.df, effect == "treatmentSC35M"))
  # FIXME duplicate point due to N-to-N protgroup to protregroup matching
  #dplyr::left_join(filter(object_effects_wide.df, std_type == "replicate") %>%
   #                  select(-gene_names, -majority_protein_acs))
require(plotly)
# plot umap + selected proteins (color = protein biological function)
umap <- ggplotly(
  ggplot(objects_umap.df,
         aes(x=x_2d, y=y_2d, color=category)) +
    geom_point(aes(text=object_label, size=!(category %in% c('default', 'NA')))) +
    #geom_text(data=filter(objects_umap.df, gene_names %in% c("RPL11", "RPL12")),
    #          aes(label=object_label), vjust=-1.1, color = "black", size = 3) +
    #geom_text_repel(data=filter(objects_umap.df, is_viral),
    #                aes(label=object_label)) +
    scale_size_manual(values = c("TRUE" = 1.0, "FALSE" = 0.5)) +
    scale_color_manual(values = c("stable" = "black", "viral" = "red", "ISG" = "blue", "hit" = "orange",
                                  "ribosome" = "darkgreen", "mito ribosome" = "darkred", "proteasome" = "khaki",
                                  "histone" = "red", "phosphatase" = "violet",
                                  "default" = "gray",
                                  "NA" = "pink", "reverse" = "pink", "contaminant" = "yellow")) +
    theme_bw_ast(base_family = "", base_size = 10),
  tooltip = "text"
)
umap

p <- ggplot(filter(msdata_full$protgroups, is_viral) %>% #str_detect(gene_names, "(^|;)(SERPIN)")) %>%
       inner_join(msdata_full$protgroup_intensities_all) %>%
       left_join(msdata_full$msruns),
       aes(x = timepoint_num, y = intensity_imputed_norm,
           color=treatment, fill=treatment)) +
  geom_smooth(alpha=0.25) +
  geom_point(aes(shape=is_imputed)) +
  scale_x_continuous(breaks = unique(msdata_full$msruns$timepoint_num)) +
  scale_shape_manual(values = c("FALSE"=16L, "TRUE"=1L)) +
  scale_color_manual(values=c("mock"="gray", "SARS_COV2"="red")) +
  scale_fill_manual(values=c("mock"="gray", "SARS_COV2"="red")) +
  facet_wrap(~ gene_names, ncol = 3, scales = "free_y") +
  theme_bw_ast()
p
ggsave(p, filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_SARS_CoV2_proteins_", data_version, ".pdf")),
       width = 12, height = 7, device = cairo_pdf)
