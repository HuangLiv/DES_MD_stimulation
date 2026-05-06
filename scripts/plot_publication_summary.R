#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- "--file="
script_path <- sub(script_arg, "", args[grep(script_arg, args)])
script_dir <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

std_dir <- file.path(root, "standard_md_analysis")
dyn_dir <- file.path(root, "dynamic_binding_analysis")
fig_dir <- file.path(root, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

reps <- paste0("rep", 1:3)
ligands <- c("CA1", "LYS1")
rep_cols <- c(rep1 = "#7F7F7F", rep2 = "#A6A6A6", rep3 = "#C7C7C7")
lig_cols <- c(CA1 = "#0072B2", LYS1 = "#D55E00")
metric_cols <- c(Protein = "#111111", CA1 = "#0072B2", LYS1 = "#D55E00")

read_xvg <- function(path, y_col = 2) {
  if (!file.exists(path)) return(data.frame())
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*[#@]", lines) & nzchar(trimws(lines))]
  if (!length(lines)) return(data.frame())
  d <- suppressWarnings(read.table(text = lines, fill = TRUE))
  d <- d[complete.cases(d), , drop = FALSE]
  if (ncol(d) < y_col) return(data.frame())
  time <- d[[1]]
  if (length(time) && max(time, na.rm = TRUE) > 1000) time <- time / 1000
  data.frame(time = time, value = d[[y_col]])
}

read_table_safe <- function(path) {
  if (!file.exists(path)) return(data.frame())
  read.delim(path, check.names = FALSE)
}

load_std_series <- function(metric, ligand = NULL, y_col = 2) {
  out <- list()
  for (rep in reps) {
    xvg_dir <- file.path(std_dir, rep, "xvg")
    fname <- switch(metric,
      protein_rmsd = "rmsd_protein.xvg",
      rg = "rg_protein.xvg",
      ligand_pose = sprintf("rmsd_%s_pose.xvg", ligand),
      ligand_self = sprintf("rmsd_%s_self.xvg", ligand),
      hbond = sprintf("hbond_protein_%s.xvg", ligand),
      stop("Unknown metric: ", metric)
    )
    d <- read_xvg(file.path(xvg_dir, fname), y_col = y_col)
    if (nrow(d)) {
      d$rep <- rep
      d$ligand <- ifelse(is.null(ligand), "Protein", ligand)
      d$metric <- metric
      out[[length(out) + 1]] <- d
    }
  }
  do.call(rbind, out)
}

load_dyn_series <- function(metric, ligand) {
  out <- list()
  for (rep in reps) {
    xvg_dir <- file.path(dyn_dir, rep, "xvg")
    fname <- switch(metric,
      mindist = sprintf("protein_%s_mindist.xvg", ligand),
      comdist = sprintf("protein_%s_comdist.xvg", ligand),
      stop("Unknown dynamic metric: ", metric)
    )
    d <- read_xvg(file.path(xvg_dir, fname))
    if (nrow(d)) {
      d$rep <- rep
      d$ligand <- ligand
      d$metric <- metric
      out[[length(out) + 1]] <- d
    }
  }
  do.call(rbind, out)
}

mean_sd <- function(d) {
  if (!nrow(d)) return(data.frame())
  m <- aggregate(value ~ time + ligand + metric, d, mean)
  s <- aggregate(value ~ time + ligand + metric, d, sd)
  names(m)[4] <- "mean"
  names(s)[4] <- "sd"
  out <- merge(m, s, by = c("time", "ligand", "metric"), all = TRUE)
  out$sd[is.na(out$sd)] <- 0
  out
}

save_plot <- function(p, name, width = 7.2, height = 4.6) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height, useDingbats = FALSE)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
}

theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_blank(),
      legend.key.size = unit(0.35, "cm"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank()
    )
}

draw_ts <- function(raw, summary, ylab, title = NULL, colors = metric_cols) {
  ggplot() +
    geom_line(data = raw, aes(time, value, group = interaction(rep, ligand), color = ligand),
              linewidth = 0.25, alpha = 0.35) +
    geom_ribbon(data = summary, aes(time, ymin = mean - sd, ymax = mean + sd, fill = ligand),
                alpha = 0.18, color = NA) +
    geom_line(data = summary, aes(time, mean, color = ligand), linewidth = 0.9) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = "Time (ns)", y = ylab, title = title) +
    theme_pub()
}

# Figure 1: protein stability.
protein_rmsd <- load_std_series("protein_rmsd")
protein_rg <- load_std_series("rg")
protein_rmsd$metric <- "Protein RMSD"
protein_rg$metric <- "Radius of gyration"
protein_stability <- rbind(protein_rmsd, protein_rg)
protein_stability$ligand <- protein_stability$metric
stab_sum <- mean_sd(protein_stability)
p1 <- ggplot() +
  geom_line(data = protein_stability, aes(time, value, group = interaction(rep, metric)),
            color = "#9E9E9E", linewidth = 0.25, alpha = 0.45) +
  geom_ribbon(data = stab_sum, aes(time, ymin = mean - sd, ymax = mean + sd),
              fill = "#4D4D4D", alpha = 0.18) +
  geom_line(data = stab_sum, aes(time, mean), color = "#111111", linewidth = 0.85) +
  facet_wrap(~metric, scales = "free_y", ncol = 1,
             labeller = as_labeller(c(`Protein RMSD` = "Protein RMSD (nm)",
                                      `Radius of gyration` = "Radius of gyration (nm)"))) +
  labs(x = "Time (ns)", y = NULL) +
  theme_pub()
save_plot(p1, "fig1_protein_stability_merged", 7.0, 5.8)

# Figure 2A: ligand RMSD.
lig_data <- list()
for (lig in ligands) {
  pose <- load_std_series("ligand_pose", lig); pose$panel <- "Pose RMSD (protein-fitted)"
  self <- load_std_series("ligand_self", lig); self$panel <- "Self RMSD"
  lig_data[[lig]] <- rbind(pose, self)
}
lig_raw <- do.call(rbind, lig_data)
lig_raw$metric <- lig_raw$panel
lig_sum <- mean_sd(lig_raw)
p2a <- ggplot() +
  geom_line(data = lig_raw, aes(time, value, group = interaction(rep, ligand), color = ligand),
            linewidth = 0.22, alpha = 0.32) +
  geom_ribbon(data = lig_sum, aes(time, ymin = pmax(mean - sd, 0), ymax = mean + sd, fill = ligand),
              alpha = 0.16, color = NA) +
  geom_line(data = lig_sum, aes(time, mean, color = ligand), linewidth = 0.85) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_color_manual(values = lig_cols) +
  scale_fill_manual(values = lig_cols) +
  labs(x = "Time (ns)", y = NULL) +
  theme_pub()
save_plot(p2a, "fig2a_ligand_rmsd_merged", 7.2, 5.2)

# Figure 2B: protein-ligand H-bonds.
hb_data <- list()
for (lig in ligands) {
  hb <- load_std_series("hbond", lig)
  hb$panel <- "Protein-ligand H-bonds"
  hb_data[[lig]] <- hb
}
hb_raw <- do.call(rbind, hb_data)
hb_raw$metric <- hb_raw$panel
hb_sum <- mean_sd(hb_raw)
p2b <- ggplot() +
  geom_line(data = hb_raw, aes(time, value, group = interaction(rep, ligand), color = ligand),
            linewidth = 0.12, alpha = 0.08) +
  geom_ribbon(data = hb_sum, aes(time, ymin = pmax(mean - sd, 0), ymax = mean + sd, fill = ligand),
              alpha = 0.12, color = NA) +
  geom_line(data = hb_sum, aes(time, mean, color = ligand), linewidth = 1.0) +
  scale_color_manual(values = lig_cols) +
  scale_fill_manual(values = lig_cols) +
  labs(x = "Time (ns)", y = "Hydrogen bonds") +
  theme_pub()
save_plot(p2b, "fig2b_protein_ligand_hbonds_merged", 7.2, 3.8)

hb_relaxed <- subset(hb_raw, time >= 20)
if (nrow(hb_relaxed)) {
  hb_rep <- aggregate(value ~ ligand + rep, hb_relaxed, mean)
  hb_mean <- aggregate(value ~ ligand, hb_rep, mean)
  hb_sd <- aggregate(value ~ ligand, hb_rep, sd)
  names(hb_mean)[2] <- "mean_hbond_20_100ns"
  names(hb_sd)[2] <- "sd_hbond_20_100ns"
  write.csv(merge(hb_mean, hb_sd, by = "ligand"),
            file.path(fig_dir, "hbond_20_100ns_summary.csv"), row.names = FALSE)
}

# Figure 3: dynamic contact distances.
dyn_list <- list()
for (lig in ligands) {
  mind <- load_dyn_series("mindist", lig); mind$panel <- "Minimum distance"
  com <- load_dyn_series("comdist", lig); com$panel <- "COM distance"
  dyn_list[[lig]] <- rbind(mind, com)
}
dyn_raw <- do.call(rbind, dyn_list)
dyn_raw$metric <- dyn_raw$panel
dyn_sum <- mean_sd(dyn_raw)
p3 <- ggplot() +
  geom_hline(data = data.frame(metric = "Minimum distance"), aes(yintercept = 0.35),
             linetype = "dashed", color = "grey40", linewidth = 0.35) +
  geom_line(data = dyn_raw, aes(time, value, group = interaction(rep, ligand), color = ligand),
            linewidth = 0.22, alpha = 0.32) +
  geom_ribbon(data = dyn_sum, aes(time, ymin = pmax(mean - sd, 0), ymax = mean + sd, fill = ligand),
              alpha = 0.16, color = NA) +
  geom_line(data = dyn_sum, aes(time, mean, color = ligand), linewidth = 0.85) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_color_manual(values = lig_cols) +
  scale_fill_manual(values = lig_cols) +
  labs(x = "Time (ns)", y = "Distance (nm)") +
  theme_pub()
save_plot(p3, "fig3_dynamic_contact_distances_merged", 7.2, 6.2)

# Figure 4: contact occupancy from dynamic summary.
dyn_summary <- read_table_safe(file.path(dyn_dir, "combined_summary.tsv"))
if (nrow(dyn_summary)) {
  occ <- dyn_summary[, c("rep", "ligand", "contact_lt_0.35", "contact_lt_0.45", "contact_lt_0.60")]
  names(occ) <- c("rep", "ligand", "0.35 nm", "0.45 nm", "0.60 nm")
  occ_long <- reshape(occ, varying = c("0.35 nm", "0.45 nm", "0.60 nm"),
                      v.names = "occupancy", timevar = "cutoff",
                      times = c("0.35 nm", "0.45 nm", "0.60 nm"),
                      direction = "long")
  occ_long$occupancy <- as.numeric(occ_long$occupancy) * 100
  occ_mean <- aggregate(occupancy ~ ligand + cutoff, occ_long, mean)
  occ_sd <- aggregate(occupancy ~ ligand + cutoff, occ_long, sd)
  names(occ_mean)[3] <- "mean"
  names(occ_sd)[3] <- "sd"
  occ_stat <- merge(occ_mean, occ_sd, by = c("ligand", "cutoff"))
  p4 <- ggplot(occ_stat, aes(cutoff, mean, fill = ligand)) +
    geom_col(position = position_dodge(width = 0.65), width = 0.58, color = "black", linewidth = 0.2) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 100)),
                  position = position_dodge(width = 0.65), width = 0.18, linewidth = 0.3) +
    geom_point(data = occ_long, aes(cutoff, occupancy, color = ligand),
               position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65),
               size = 1.4, alpha = 0.75, inherit.aes = FALSE) +
    scale_fill_manual(values = lig_cols) +
    scale_color_manual(values = lig_cols) +
    labs(x = "Contact cutoff", y = "Contact occupancy (%)") +
    ylim(0, 105) +
    theme_pub()
  save_plot(p4, "fig4_contact_occupancy_bar", 5.8, 4.0)
}

# Figure 5: protein RMSF using continuous residue numbering across the trimer.
rmsf_files <- list.files(std_dir, pattern = "rmsf_.*\\.csv$", recursive = TRUE, full.names = TRUE)
rmsf <- do.call(rbind, lapply(rmsf_files, read.csv))
if (!is.null(rmsf) && nrow(rmsf)) {
  rmsf$rep <- factor(rmsf$rep, levels = reps)
  protein_rmsf <- subset(rmsf, group == "Protein")
  if (nrow(protein_rmsf)) {
    protein_rmsf$segment_num <- suppressWarnings(as.integer(as.character(protein_rmsf$segment)))
    if (any(is.na(protein_rmsf$segment_num))) {
      protein_rmsf$segment_num <- as.integer(factor(protein_rmsf$segment, levels = unique(protein_rmsf$segment)))
    }
    seg_max <- aggregate(index ~ segment_num, protein_rmsf, max)
    seg_max <- seg_max[order(seg_max$segment_num), ]
    seg_max$offset <- c(0, head(cumsum(seg_max$index), -1))
    names(seg_max)[2] <- "segment_max_index"
    protein_rmsf <- merge(protein_rmsf, seg_max[, c("segment_num", "offset")], by = "segment_num")
    protein_rmsf$residue_id <- protein_rmsf$index + protein_rmsf$offset
    prot_mean <- aggregate(rmsf_nm ~ residue_id, protein_rmsf, mean)
    prot_sd <- aggregate(rmsf_nm ~ residue_id, protein_rmsf, sd)
    names(prot_mean)[2] <- "mean"
    names(prot_sd)[2] <- "sd"
    prot_stat <- merge(prot_mean, prot_sd, by = "residue_id")
    p5 <- ggplot() +
      geom_vline(xintercept = c(27.5, 54.5), linetype = "dashed", color = "grey70", linewidth = 0.3) +
      geom_line(data = protein_rmsf, aes(residue_id, rmsf_nm, group = rep),
                color = "#B0B0B0", linewidth = 0.25, alpha = 0.45) +
      geom_ribbon(data = prot_stat, aes(residue_id, ymin = pmax(mean - sd, 0), ymax = mean + sd),
                  fill = "#56B4E9", alpha = 0.20) +
      geom_line(data = prot_stat, aes(residue_id, mean), color = "#0072B2", linewidth = 0.8) +
      scale_x_continuous(breaks = seq(1, max(prot_stat$residue_id), by = 9)) +
      labs(x = "Continuous residue index across trimer", y = "Protein RMSF (nm)") +
      theme_pub()
    save_plot(p5, "fig5_protein_rmsf_continuous", 7.2, 3.6)
  }
}

# Summary tables.
std_summary <- read_table_safe(file.path(std_dir, "combined_summary.tsv"))
if (nrow(std_summary)) {
  std_num <- std_summary
  numeric_cols <- setdiff(names(std_num), c("rep", "metric", "group"))
  for (cc in numeric_cols) std_num[[cc]] <- suppressWarnings(as.numeric(std_num[[cc]]))
  m <- aggregate(mean ~ metric + group, std_num, mean, na.rm = TRUE)
  s <- aggregate(mean ~ metric + group, std_num, sd, na.rm = TRUE)
  names(m)[3] <- "mean_across_reps"
  names(s)[3] <- "sd_across_reps"
  write.csv(merge(m, s, by = c("metric", "group")), file.path(fig_dir, "standard_summary_across_replicates.csv"), row.names = FALSE)
}
if (nrow(dyn_summary)) {
  dyn_num <- dyn_summary
  for (cc in setdiff(names(dyn_num), c("rep", "ligand"))) dyn_num[[cc]] <- suppressWarnings(as.numeric(dyn_num[[cc]]))
  metrics <- c("com_mean_nm", "mindist_mean_nm", "contact_lt_0.35", "contact_lt_0.45", "contact_lt_0.60")
  out <- list()
  for (mm in metrics) {
    m <- aggregate(dyn_num[[mm]], list(ligand = dyn_num$ligand), mean, na.rm = TRUE)
    s <- aggregate(dyn_num[[mm]], list(ligand = dyn_num$ligand), sd, na.rm = TRUE)
    out[[mm]] <- data.frame(ligand = m$ligand, metric = mm, mean = m$x, sd = s$x)
  }
  write.csv(do.call(rbind, out), file.path(fig_dir, "dynamic_summary_across_replicates.csv"), row.names = FALSE)
}

cat("Publication figures written to:", fig_dir, "\n")
