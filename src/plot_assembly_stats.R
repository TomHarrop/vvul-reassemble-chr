#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(ggplot2)

###########
# GLOBALS #
###########

stats_files <- snakemake@input[["stats_files"]]
plot_file <- snakemake@output[["plot"]]

# dev
# stats_files <- list.files("output/060_stats",
#            recursive = TRUE,
#            pattern = "stats.tsv",
#            full.names = TRUE)
# plot_file <- "test/stats.pdf"

########
# MAIN #
########

# read stats
names(stats_files) <- basename(dirname(stats_files))
stats_list <- lapply(stats_files, fread)
stats <- rbindlist(stats_list, idcol = "meraculous_run")

# mung
stats[, c("chr",
          "version",
          "diplo",
          "k") :=
          tstrsplit(meraculous_run, ".", fixed = TRUE)]

# fill in missing values during melt. key order has to match order in CJ()
setkey(stats, chr, version, diplo, k)
pd <- melt(stats[CJ(unique(chr),
                    unique(version),
                    unique(diplo),
                    unique(k))],
           id.vars = c("chr", "version", "diplo", "k"),
           measure.vars = c("n_scaffolds", "scaf_bp", "scaf_N50", "scaf_L50"),
           fill = TRUE)

# set up labels
variable_order <- c("n_scaffolds" = "Scaffolds",
                    "scaf_bp" = "Assembled size (Kb)",
                    "scaf_N50" = "N50",
                    "scaf_L50" = "L50")

pd[, value := as.double(value)]
pd[variable == "scaf_bp", value := value / 1000]


pd[, k := as.numeric(gsub("[^[:digit:]]+", "", k))]
pd[, diplo := as.numeric(gsub("[^[:digit:]]+", "", diplo))]

pd[, variable := factor(plyr::revalue(variable, variable_order),
                        levels = variable_order)]

# set up colours
fill_colours = viridis::viridis_pal()(3)

# line for max per value
best_dt <- pd[, .(value = ifelse(
    variable %in% c("Scaffolds", "N50"),
    min(value, na.rm = TRUE),
    max(value, na.rm = TRUE))),
    by = variable]

# draw the plot
gp <- ggplot(pd, aes(x = as.factor(k),
                     y = value,
                     fill = as.factor(diplo))) +
    theme_grey(base_size = 8) +
    coord_flip() +
    theme(strip.placement = "outside",
          strip.background.x = element_blank(),
          strip.text.y = element_text(angle = 0)) +
    facet_grid(chr ~ variable, scales = "free", switch = "x") +
    xlab(expression(italic("k"))) + ylab(NULL) +
    scale_fill_manual(values = fill_colours[c(1,3)],
                      guide = guide_legend(title = "Diplo")) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.5) +
    geom_hline(mapping = aes(yintercept = value),
               data = best_dt,
               colour = fill_colours[2])

# write output
wo <- grid::convertUnit(grid::unit(483, "pt"), "mm", valueOnly = TRUE)
ho <- grid::convertUnit(grid::unit(664/3, "pt"), "mm", valueOnly = TRUE)
ggsave(filename = plot_file,
       plot = gp,
       device = cairo_pdf,
       width = wo,
       height = ho,
       units = "mm")

# write session info
sessionInfo()
