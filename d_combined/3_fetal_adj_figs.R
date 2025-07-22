###############################################################################
#                                                                             #
#          Plots comparing all MR estimates for main and fetal adjusted       #
#                                                                             #
###############################################################################

###############################################################################
#                                   Set up                                    #
###############################################################################

# Clear the work environment
rm(list=ls())

# Setting digits
options(digits = 10)

# Required libraries
## install.packages("paletteer")
x <- c("forestplot", "dplyr", "tidyr", "grid", "paletteer", "stringr")
lapply(x, require, character.only = TRUE)

# Set directories
home_dir <- paste0(Sys.getenv("AGE_AT_MENARCHE_DIR"), "/working/")
mr_results_dir <- paste0(home_dir, "/results/b_mr/")
mvmr_results_dir <- paste0(home_dir, "/results/c_mvmr/")
setwd(home_dir)
out_dir <- paste0(home_dir, "/results/d_combined/forest_plots/")
tables_out_dir <- paste0(home_dir, "/results/d_combined/tables/")

# Outcome names
outcome_names_df <- read.csv(paste0(home_dir, "/data/b_mr/outcome/outcomes.csv"))
outcomes <- c(outcome_names_df$outcome)
outcomes_not_in_fetal_adj <- c("Low birthweight", "Very preterm birth")

# Adapt long ones for plots:
outcome_names_df[1, 2] <- "Any hypertensive\ndisorders of pregnancy"
outcome_names_df[6, 2] <- "Small for gestational age "
outcome_names_df[8, 2] <- "Gestational diabetes\nmellitus"
outcome_names_df[10, 2] <- "Large for gestational age "

###############################################################################
#                          Read and clean results                             #
###############################################################################

### Raw files
univariate_mr_binary <- read.csv(paste0(mr_results_dir, "mr_results_binary.csv"))
univariate_mr_bw <- read.csv(paste0(mr_results_dir, "mr_results_bw.csv"))
multivariable_mr_binary <- read.csv(paste0(mvmr_results_dir, "mvmr_results_binary.csv"))
multivariable_mr_bw <- read.csv(paste0(mvmr_results_dir, "mvmr_results_bw.csv"))

### Univariate MR
mr_binary <- univariate_mr_binary %>%
  # update outcome names
  select(-outcome_clean) %>%
  left_join(outcome_names_df) %>%
  # restrict to IVW, main & fetal adjusted
  filter(method == "Inverse variance weighted" &
           (analysis == "Main" | analysis == "Adjusted for fetal genotype")) %>%
  # remove outcomes where fetal genotype not available
  filter(!(outcome_clean %in% outcomes_not_in_fetal_adj)) %>%
  mutate(
    mean = or,
    lower = or_lci95,
    upper = or_uci95
  )
mr_binary$`OR (95% CI)` <- mr_binary$`OR..95..CI.`
mr_binary <- mr_binary %>% 
  select(outcome_clean, analysis, mean, lower, upper, `OR (95% CI)`)

mr_bw <- univariate_mr_bw %>%
  # update outcome names
  select(-outcome_clean) %>%
  left_join(outcome_names_df) %>%
  # restrict to IVW, main & fetal adjusted
  filter(method == "Inverse variance weighted" &
  (analysis == "Main" | analysis == "Adjusted for fetal genotype")) %>%
  # remove outcomes where fetal genotype not available
  filter(!(outcome_clean %in% outcomes_not_in_fetal_adj)) %>%
  mutate(
    mean = b,
    lower = lo_ci,
    upper = up_ci,
    approach = "UVMR",
    # generate clean beta label to 3 dp
    `Beta (95% CI)` = paste0(format(round(b, 3), nsmall = 3), " (",
                             format(round(lo_ci, 3), nsmall = 3), " - ",
                             format(round(up_ci, 3), nsmall = 3), ")"),
  )
mr_bw <- mr_bw %>%
  select(outcome_clean, approach, analysis, mean, lower, upper, `Beta (95% CI)`)

### Multivariable MR
mvmr_binary <- multivariable_mr_binary %>%
  # update outcome names
  select(-outcome_clean) %>%
  left_join(outcome_names_df) %>%
  filter(id.exposure == "ieu-b-5136" &
    (analysis == "Main" | analysis == "Adjusted for fetal genotype")) %>%
  # remove outcomes where fetal genotype not available
  filter(!(outcome_clean %in% outcomes_not_in_fetal_adj)) %>%
  mutate(
      mean = or,
      lower = or_lci95,
      upper = or_uci95
    )
mvmr_binary$`OR (95% CI)` <- mvmr_binary$`OR..95..CI.`
mvmr_binary <- mvmr_binary %>%
  select(outcome_clean, analysis, mean, lower, upper, `OR (95% CI)`)

mvmr_bw <- multivariable_mr_bw %>%
  # update outcome names
  select(-outcome_clean) %>%
  left_join(outcome_names_df) %>%
  filter(id.exposure == "ieu-b-5136" &
    (analysis == "Main" | analysis == "Adjusted for fetal genotype")) %>%
  # remove outcomes where fetal genotype not available
  filter(!(outcome_clean %in% outcomes_not_in_fetal_adj)) %>%
  mutate(
      mean = b,
      lower = lo_ci,
      upper = up_ci,
      approach = "MVMR",
      # generate clean beta label
      `Beta (95% CI)` = paste0(format(round(b, 3), nsmall = 3), " (",
                               format(round(lo_ci, 3), nsmall = 3), " - ",
                               format(round(up_ci, 3), nsmall = 3), ")"),
  )
mvmr_bw <- mvmr_bw %>%
  select(outcome_clean, approach, analysis, mean, lower, upper, `Beta (95% CI)`)

###############################################################################
#           Birthweight results table (univariate & MVMR models)              #
###############################################################################

### Table summarising for continuous outcome (birthweight)
bw_table <- rbind(mr_bw, mvmr_bw) %>%
  mutate(approach = factor(approach, levels = c("UVMR",
    "MVMR")))

bw_table %>%
  select(approach, analysis, `Beta (95% CI)`) %>%
  arrange(approach, desc(analysis)) %>%
  write.csv(paste0(tables_out_dir, "bw_fetal_adj.csv"), row.names = FALSE)

### Quick plot
library(ggplot2)
plot <- bw_table %>%
  ggplot(aes(y = approach, x = mean, xmin = lower, xmax = upper, colour = analysis)) +
  geom_pointrange(position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Beta per 1 year increase in age at menarche") +
  theme_bw()
ggsave(plot, file = paste0(out_dir, "bw_fetal_adj_26032025.pdf"),
 width = 9, height = 5)

###############################################################################
#                             Plotting - MR IVW                               #
###############################################################################

### prepare data ########################################

binary_df <- mr_binary %>%
  mutate(`OR (95% CI)` = str_replace(`OR (95% CI)`, " - ", ", "))

# add empty rows and order by outcome, within that approach
empty_rows <- data.frame(outcome_clean = unique(binary_df$outcome_clean))
binary_df <- plyr::rbind.fill(empty_rows, binary_df) %>%
  arrange(match(outcome_clean, outcome_names_df$outcome_clean),
          match(analysis, c(NA, "Main", "Adjusted for fetal genotype"))) %>%
# replace outcome_clean with "" except for row 2 where we want the label
  mutate(
    outcome_clean = case_when(analysis == "Main" ~ outcome_clean,
    .default = ""))  %>%
  select(outcome_clean, `OR (95% CI)`, mean, lower, upper)
# headers
binary_df[1,] <- c("Outcome", "Odds ratio (95% CI)", NA, NA, NA)

# finalise text columns, add space for legend
alltext <- binary_df %>% select(outcome_clean, `OR (95% CI)`) %>%
  mutate(space = "                            ")
alltext <- alltext %>% replace(is.na(.), "")

# numeric results columns
resforalltext <- binary_df %>% select(mean, lower, upper)
resforalltext <- sapply(resforalltext, as.numeric)

### plot #############################################################

# UVMR IVW colour from main plot + orange for adjusted
# paletteer_d("rcartocolor::Temps")[1]
palette <- c("#2A788EFF", "#F18805")

styles <- fpShapesGp(
  box =
    rep(list(
      gpar(col = palette[1], fill = palette[1]), # for NA lines
      gpar(col = palette[1], fill = palette[1]),
      gpar(col = palette[2], fill = palette[2])
    ), 10), # number of outcomes
  lines = rep(list(
    gpar(col=palette[1]), # for NA lines
    gpar(col=palette[1]),
    gpar(col=palette[2])
  ), 10)
)

ticks <- c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4)

all <- resforalltext %>%
  forestplot(labeltext = alltext,
             txt_gp = fpTxtGp(
              label = gpar(cex = .8),
              ticks = gpar(cex = .8),
              xlab = gpar(cex = .9)),
             is.summary=c(TRUE, rep(FALSE,32)),
             align = c("l", "l", "r"), 
             graph.pos = 3,
             boxsize = 0.3,
             fn.ci_norm = fpDrawCircleCI,
             xlog = TRUE,
             xticks=ticks,
             xlab = "\nOdds ratio per 1 year increase in age at menarche",
             zero=1,
             col = fpColors(zero="grey50"),
             colgap=unit(0.01, "npc"),
             shapes_gp = styles) |>
  fp_decorate_graph(grid =
    structure(ticks[-c(1,length(ticks))],
    gp = gpar(lty = 2, col = "grey50")))
print(all)

tiff(paste0(out_dir, "supp_fig9_260325.tiff"), units="in", width=11, height=9, res=300)
plot.new()
all
legend("right",
       c("MR IVW",
         " ",
         "MR IVW adjusted\nfor fetal genotype",
         " "),
       title=as.expression(bquote(bold("Method"))),
       border="grey50", box.lwd=1,
       cex = 0.85,
       col=c(palette[1], "white", palette[2],"white"),
       pch=c(16))
dev.off()

###############################################################################
#                             Plotting - MVMR                                 #
###############################################################################

### prepare data ########################################

binary_df <- mvmr_binary %>%
  mutate(`OR (95% CI)` = str_replace(`OR (95% CI)`, " - ", ", "))

# add empty rows and order by outcome, within that approach
empty_rows <- data.frame(outcome_clean = unique(binary_df$outcome_clean))
binary_df <- plyr::rbind.fill(empty_rows, binary_df) %>%
  arrange(match(outcome_clean, outcome_names_df$outcome_clean),
          match(analysis, c(NA, "Main", "Adjusted for fetal genotype"))) %>%
# replace outcome_clean with "" except for row 2 where we want the label
  mutate(
    outcome_clean = case_when(analysis == "Main" ~ outcome_clean,
    .default = "")) %>%
  select(outcome_clean, `OR (95% CI)`, mean, lower, upper)
# headers
binary_df[1,] <- c("Outcome", "Odds ratio (95% CI)", NA, NA, NA)

# finalise text columns, add space for legend
alltext <- binary_df %>% select(outcome_clean, `OR (95% CI)`) %>%
  mutate(space = "                             ")
alltext <- alltext %>% replace(is.na(.), "")

# numeric results columns
resforalltext <- binary_df %>% select(mean, lower, upper)
resforalltext <- sapply(resforalltext, as.numeric)

### plot #############################################################

# MVMR IVW colour from main plot + orange for adjusted
palette <- c("#440154FF", "#F18805")

styles <- fpShapesGp(
  box =
    rep(list(
      gpar(col = palette[1], fill = palette[1]), # for NA lines
      gpar(col = palette[1], fill = palette[1]),
      gpar(col = palette[2], fill = palette[2])
    ), 10), # number of outcomes
  lines = rep(list(
    gpar(col=palette[1]), # for NA lines
    gpar(col=palette[1]),
    gpar(col=palette[2])
  ), 10)
)

ticks <- c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4)

all <- resforalltext %>%
  forestplot(labeltext = alltext,
             txt_gp = fpTxtGp(
              label = gpar(cex = .8),
              ticks = gpar(cex = .8),
              xlab = gpar(cex = .9)),
             is.summary=c(TRUE, rep(FALSE,29)),
             align = c("l", "l", "r"), 
             graph.pos = 3,
             boxsize = 0.3,
             fn.ci_norm = fpDrawCircleCI,
             xlog = TRUE,
             xticks=ticks,
             xlab = "\nOdds ratio per 1 year increase in age at menarche",
             zero=1,
             col = fpColors(zero="grey50"),
             colgap=unit(0.01, "npc"),
             shapes_gp = styles) |>
  fp_decorate_graph(grid =
    structure(ticks[-c(1,length(ticks))],
    gp = gpar(lty = 2, col = "grey50")))
print(all)

tiff(paste0(out_dir, "supp_fig10_260325.tiff"), units="in", width=11, height=9, res=300)
plot.new()
all
legend("right",
       c("MVMR",
         " ",
         "MVMR adjusted\nfor fetal genotype",
         " "),
       title=as.expression(bquote(bold("Method"))),
       border="grey50", box.lwd=1,
       cex = 0.85,
       col=c(palette[1], "white", palette[2],"white"),
       pch=c(16))
dev.off()

q('no')