###############################################################################
#                                                                             #
#                Create plots triangulating results across                    #
#                all regression models, MR & MVMR analyses                    #
#                                                                             #
###############################################################################

# fig 1 - HDP
# fig 2 - smaller size & preterm
# fig 3 - larger size & post term
# fig 4 - GDM
# supp fig 1 - continuous birthweight
# supp fig 2 - perinatal depression

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
reg_results_dir <- paste0(home_dir, "/results/a_regression/main/")
mr_results_dir <- paste0(home_dir, "/results/b_mr/")
mvmr_results_dir <- paste0(home_dir, "/results/c_mvmr/")
setwd(home_dir)
out_dir <- paste0(home_dir, "/results/d_combined/forest_plots/")

# Outcome names
outcome_names_df <- read.csv(paste0(home_dir, "/data/b_mr/outcome/outcomes.csv"))
outcomes <- c(outcome_names_df$outcome)

# Adapt long ones for plots:
outcome_names_df[1, 2] <- "Any hypertensive\ndisorders of pregnancy"
outcome_names_df[6, 2] <- "Small for gestational age "
outcome_names_df[8, 2] <- "Gestational diabetes\nmellitus"
outcome_names_df[10, 2] <- "Large for gestational age "

# Approach order
approach_order <- c(NA, "Regression M1", "Regression M2", "Regression M3",
                    "MR IVW", "MVMR, adiposity adjusted",
                    "MR Egger", "MR Weighted median", "MR Weighted mode")

###############################################################################
#                          Read and clean results                             #
###############################################################################

### Raw files
regression_binary <- read.csv(paste0(reg_results_dir, "linear_exp_binary_out.csv"))
regression_bw <- read.csv(paste0(reg_results_dir, "linear_exp_cont_out.csv"))
univariate_mr_binary <- read.csv(paste0(mr_results_dir, "mr_results_binary.csv"))
univariate_mr_bw <- read.csv(paste0(mr_results_dir, "mr_results_bw.csv"))
multivariable_mr_binary <- read.csv(paste0(mvmr_results_dir, "mvmr_results_binary.csv"))
multivariable_mr_bw <- read.csv(paste0(mvmr_results_dir, "mvmr_results_bw.csv"))

### regression
reg_binary <- regression_binary %>%
  # clean approach name
  mutate(
    approach = case_when(
      model_name == "raw" ~ "Regression M1",
      model_name == "adjusted" ~ "Regression M2",
      model_name == "adjusted_bmi" ~ "Regression M3"
    )) %>%
  # update outcome names
  select(-outcome_clean) %>%
  left_join(outcome_names_df) %>%
  # ORs
  mutate(
    or = beta_OR,
    or_lci95 = conf.low,
    or_uci95 = conf.high,
    `Cases / controls` = paste0(model_ncases, " / ", model_ncontrols),
    # generate clean OR label
    `OR (95% CI)` = paste0(format(round(or, 2), nsmall = 2), " (",
                           format(round(or_lci95, 2), nsmall = 2), " - ",
                           format(round(or_uci95, 2), nsmall = 2), ")")
  ) %>% 
  select(outcome_clean, approach, or, or_lci95, or_uci95, `OR (95% CI)`,
         `Cases / controls`)

reg_bw <- regression_bw %>%
  # clean approach name
  mutate(
    approach = case_when(
      model_name == "raw" ~ "Regression M1",
      model_name == "adjusted" ~ "Regression M2",
      model_name == "adjusted_bmi" ~ "Regression M3"
    )) %>%
  # update outcome names
  select(-outcome_clean) %>%
  left_join(outcome_names_df) %>%
  mutate(
    `Sample size` = paste0(model_nobs),
    b = estimate,
    lo_ci = conf.low,
    up_ci = conf.high,
    # generate clean beta label
    `Beta (95% CI)` = paste0(format(round(b, 3), nsmall = 3), " (",
                             format(round(lo_ci, 3), nsmall = 3), " - ",
                             format(round(up_ci, 3), nsmall = 3), ")"),
  ) %>% 
  # select columns needed
  select(outcome_clean, approach, b, lo_ci, up_ci, `Beta (95% CI)`)

### Univariate MR
mr_binary <- univariate_mr_binary %>%
  # update outcome names
  select(-outcome_clean) %>%
  left_join(outcome_names_df) %>%
  # restrict to IVW - all 4 methods
  filter(analysis == "Main") %>%
  mutate(approach = method) %>%
  # clean approach name
  mutate(
  approach = case_when(
    approach == "Inverse variance weighted" ~ "MR IVW",
    approach == "Weighted median" ~ "MR Weighted median",
    approach == "Weighted mode" ~ "MR Weighted mode",
    .default = approach
  ))
mr_binary$`OR (95% CI)` <- mr_binary$`OR..95..CI.`
mr_binary$`Cases / controls` <- mr_binary$`Cases...controls`
# select coluumns needed
mr_binary <- mr_binary %>% 
  select(outcome_clean, approach, or, or_lci95, or_uci95, `OR (95% CI)`,
         `Cases / controls`)

mr_bw <- univariate_mr_bw %>%
  # update outcome names
  select(-outcome_clean) %>%
  left_join(outcome_names_df) %>%
  # restrict to IVW - all 4 methods
  filter(analysis == "Main") %>%
  mutate(approach = method) %>%
           # clean approach name
           mutate(
             approach = case_when(
               approach == "Inverse variance weighted" ~ "MR IVW",
               approach == "Weighted median" ~ "MR Weighted median",
               approach == "Weighted mode" ~ "MR Weighted mode",
               .default = approach
             )) %>% 
          mutate(
         # generate clean beta label
         `Beta (95% CI)` = paste0(format(round(b, 3), nsmall = 3), " (",
                                  format(round(lo_ci, 3), nsmall = 3), " - ",
                                  format(round(up_ci, 3), nsmall = 3), ")"),
  ) %>%
  # select columns needed
  select(outcome_clean, approach, b, lo_ci, up_ci, `Beta (95% CI)`)

### Multivariable MR
mvmr_binary <- multivariable_mr_binary %>%
  # update outcome names
  select(-outcome_clean) %>%
  left_join(outcome_names_df) %>%
  filter(id.exposure == "ieu-b-5136" & analysis == "Main" & model == "model_1") %>%
  mutate(approach = "MVMR, adiposity adjusted")
mvmr_binary$`OR (95% CI)` <- mvmr_binary$`OR..95..CI.`
mvmr_binary$`Cases / controls` <- mvmr_binary$`Cases...control`
mvmr_binary <- mvmr_binary %>%
  # select columns needed
  select(outcome_clean, approach, or, or_lci95, or_uci95, `OR (95% CI)`,
         `Cases / controls`)

mvmr_bw <- multivariable_mr_bw %>%
  # update outcome names
  select(-outcome_clean) %>%
  left_join(outcome_names_df) %>%
  filter(id.exposure == "ieu-b-5136" & analysis == "Main" & model == "model_1") %>%
  mutate(approach = "MVMR, adiposity adjusted",
         # generate clean beta label
         `Beta (95% CI)` = paste0(format(round(b, 3), nsmall = 3), " (",
                                  format(round(lo_ci, 3), nsmall = 3), " - ",
                                  format(round(up_ci, 3), nsmall = 3), ")"),
  ) %>%
  # select columns needed
  select(outcome_clean, approach, b, lo_ci, up_ci, `Beta (95% CI)`)

###############################################################################
#                                  Palette                                    #
###############################################################################

# https://waldyrious.net/viridis-palette-generator/
# set styles
palette <- c(
  rep("#7AD151FF", 3), # Regression
  "#2A788EFF", # MR IVW
  "#440154FF", # MVMR, adiposity adjusted
  rep("#2A788EFF", 3)) # MR Egger / Weighted median / Weighted mode

###############################################################################
#                         Fig 1-4, S2 - Binary                                #
###############################################################################

### prepare data for plotting - binary ########################################

binary_df <- rbind(reg_binary, mr_binary, mvmr_binary) %>%
  mutate(`OR (95% CI)` = str_replace(`OR (95% CI)`, " - ", ", "))
binary_df$mean <- binary_df$or
binary_df$lower <- binary_df$or_lci95
binary_df$upper <- binary_df$or_uci95

# add empty rows and order by outcome, within that approach
empty_rows <- data.frame(outcome_clean = unique(binary_df$outcome_clean))
binary_df <- plyr::rbind.fill(empty_rows, binary_df) %>%
  arrange(match(outcome_clean, outcome_names_df$outcome_clean),
          match(approach, approach_order)) %>%
  # replace outcome_clean with "" except for row 4 where we want the label
  mutate(
    outcome_clean = case_when(approach == "MR IVW" ~ outcome_clean,
                              .default = ""))  %>%
  select(outcome_clean, approach, `OR (95% CI)`, mean, lower, upper)
# headers
binary_df[1,] <- c("Outcome", "Approach", "Odds ratio (95% CI)", NA, NA, NA)

# finalise text columns, add space for legend
alltext <- binary_df %>% select(outcome_clean, approach, `OR (95% CI)`)
alltext <- alltext %>% replace(is.na(.), "")

# numeric results columns
resforalltext <- binary_df %>% select(mean, lower, upper)
resforalltext <- sapply(resforalltext, as.numeric)

# 120 lines long!

# fig 1 - Any HDP, GH, PE #####################################################
ticks <- c(0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4)

styles <- fpShapesGp(
  box =
    rep(list(
      gpar(col = palette[1], fill = palette[1]), # for NA lines
      gpar(col = palette[1], fill = palette[1]),
      gpar(col = palette[2], fill = palette[2]),
      gpar(col = palette[3], fill = palette[3]),
      gpar(col = palette[4], fill = palette[4]),
      gpar(col = palette[5], fill = palette[5]),
      gpar(col = palette[6], fill = palette[6]),
      gpar(col = palette[7], fill = palette[7]),
      gpar(col = palette[8], fill = palette[8])
    ), 3), # number of outcomes
  lines = rep(list(
    gpar(col=palette[1]), # for NA lines
    gpar(col=palette[1]),
    gpar(col=palette[2]),
    gpar(col=palette[3]),
    gpar(col=palette[4]),
    gpar(col=palette[5]),
    gpar(col=palette[6]),
    gpar(col=palette[7]),
    gpar(col=palette[8])
  ), 3)
)

all <- resforalltext[1:27,] %>%
  forestplot(labeltext = alltext[1:27,],
             txt_gp = fpTxtGp(
               label = gpar(cex = .8),
               ticks = gpar(cex = .8),
               xlab = gpar(cex = .9)), #fontface = "bold")),
             is.summary=c(TRUE, rep(FALSE,71)),
             align = c("l", "l", "r"), 
             graph.pos = 3,
             boxsize = 0.3,
             #clip = c(-0.65, 1.3),
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

tiff(paste0(out_dir, "fig1_040725.tiff"), units="in", width=11, height=7, res=300)
plot.new()
all
dev.off()

# fig 3 - PTB, V PTB, SGA, LBW #################################################
#ticks <- c(0.4, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.8, 2.0) # wider than others
ticks <- c(0.5, 0.75, 1.00, 1.25, 1.5, 1.75, 2)

styles <- fpShapesGp(
  box =
    rep(list(
      gpar(col = palette[1], fill = palette[1]), # for NA lines
      gpar(col = palette[1], fill = palette[1]),
      gpar(col = palette[2], fill = palette[2]),
      gpar(col = palette[3], fill = palette[3]),
      gpar(col = palette[4], fill = palette[4]),
      gpar(col = palette[5], fill = palette[5]),
      gpar(col = palette[6], fill = palette[6]),
      gpar(col = palette[7], fill = palette[7]),
      gpar(col = palette[8], fill = palette[8])
    ), 4), # number of outcomes
  lines = rep(list(
    gpar(col=palette[1]), # for NA lines
    gpar(col=palette[1]),
    gpar(col=palette[2]),
    gpar(col=palette[3]),
    gpar(col=palette[4]),
    gpar(col=palette[5]),
    gpar(col=palette[6]),
    gpar(col=palette[7]),
    gpar(col=palette[8])
  ), 4)
)

all <- resforalltext[c(1,29:63),] %>%
  forestplot(labeltext = alltext[c(1,29:63),],
             txt_gp = fpTxtGp(
               label = gpar(cex = .8),
               ticks = gpar(cex = .8),
               xlab = gpar(cex = .9)), #fontface = "bold")),
             is.summary=c(TRUE, rep(FALSE,71)),
             align = c("l", "l", "r"), 
             graph.pos = 3,
             boxsize = 0.3,
             #clip = c(-0.65, 1.3),
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

tiff(paste0(out_dir, "fig3_040725.tiff"), units="in", width=11, height=8.5, res=300)
plot.new()
all
dev.off()

# fig 4 - post-term birth, LGA, HBW ###########################################
ticks <- c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4)
styles <- fpShapesGp(
  box =
    rep(list(
      gpar(col = palette[1], fill = palette[1]), # for NA lines
      gpar(col = palette[1], fill = palette[1]),
      gpar(col = palette[2], fill = palette[2]),
      gpar(col = palette[3], fill = palette[3]),
      gpar(col = palette[4], fill = palette[4]),
      gpar(col = palette[5], fill = palette[5]),
      gpar(col = palette[6], fill = palette[6]),
      gpar(col = palette[7], fill = palette[7]),
      gpar(col = palette[8], fill = palette[8])
    ), 3), # number of outcomes
  lines = rep(list(
    gpar(col=palette[1]), # for NA lines
    gpar(col=palette[1]),
    gpar(col=palette[2]),
    gpar(col=palette[3]),
    gpar(col=palette[4]),
    gpar(col=palette[5]),
    gpar(col=palette[6]),
    gpar(col=palette[7]),
    gpar(col=palette[8])
  ), 3)
)

all <- resforalltext[c(1,74:99),] %>%
  forestplot(labeltext = alltext[c(1,74:99),],
             txt_gp = fpTxtGp(
               label = gpar(cex = .8),
               ticks = gpar(cex = .8),
               xlab = gpar(cex = .9)), #fontface = "bold")),
             is.summary=c(TRUE, rep(FALSE,71)),
             align = c("l", "l", "r"), 
             graph.pos = 3,
             boxsize = 0.3,
             #clip = c(-0.65, 1.3),
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

tiff(paste0(out_dir, "fig4_040725.tiff"), units="in", width=11, height=7, res=300)
plot.new()
all
dev.off()

# fig 2 - GDM #################################################################
ticks <- c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4)
styles <- fpShapesGp(
  box =
    rep(list(
      gpar(col = palette[1], fill = palette[1]), # for NA lines
      gpar(col = palette[1], fill = palette[1]),
      gpar(col = palette[2], fill = palette[2]),
      gpar(col = palette[3], fill = palette[3]),
      gpar(col = palette[4], fill = palette[4]),
      gpar(col = palette[5], fill = palette[5]),
      gpar(col = palette[6], fill = palette[6]),
      gpar(col = palette[7], fill = palette[7]),
      gpar(col = palette[8], fill = palette[8])
    ), 1), # number of outcomes
  lines = rep(list(
    gpar(col=palette[1]), # for NA lines
    gpar(col=palette[1]),
    gpar(col=palette[2]),
    gpar(col=palette[3]),
    gpar(col=palette[4]),
    gpar(col=palette[5]),
    gpar(col=palette[6]),
    gpar(col=palette[7]),
    gpar(col=palette[8])
  ), 1)
)

all <- resforalltext[c(1,65:72),] %>%
  forestplot(labeltext = alltext[c(1,65:72),],
             txt_gp = fpTxtGp(
               label = gpar(cex = .8),
               ticks = gpar(cex = .8),
               xlab = gpar(cex = .9)), #fontface = "bold")),
             is.summary=c(TRUE, rep(FALSE,71)),
             align = c("l", "l", "r"), 
             graph.pos = 3,
             boxsize = 0.25,
             #clip = c(-0.65, 1.3),
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

tiff(paste0(out_dir, "fig2_040725.tiff"), units="in", width=11, height=3.5, res=300)
plot.new()
all
dev.off()

# fig S4 - Perinatal depr #################################################################
ticks <- c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)
styles <- fpShapesGp(
  box =
    rep(list(
      gpar(col = palette[1], fill = palette[1]), # for NA lines
      gpar(col = palette[1], fill = palette[1]),
      gpar(col = palette[2], fill = palette[2]),
      gpar(col = palette[3], fill = palette[3]),
      gpar(col = palette[4], fill = palette[4]),
      gpar(col = palette[5], fill = palette[5]),
      gpar(col = palette[6], fill = palette[6]),
      gpar(col = palette[7], fill = palette[7]),
      gpar(col = palette[8], fill = palette[8])
    ), 1), # number of outcomes
  lines = rep(list(
    gpar(col=palette[1]), # for NA lines
    gpar(col=palette[1]),
    gpar(col=palette[2]),
    gpar(col=palette[3]),
    gpar(col=palette[4]),
    gpar(col=palette[5]),
    gpar(col=palette[6]),
    gpar(col=palette[7]),
    gpar(col=palette[8])
  ), 1)
)

all <- resforalltext[c(1,101:108),] %>%
  forestplot(labeltext = alltext[c(1,101:108),],
             txt_gp = fpTxtGp(
               label = gpar(cex = .8),
               ticks = gpar(cex = .8),
               xlab = gpar(cex = .9)), #fontface = "bold")),
             is.summary=c(TRUE, rep(FALSE,71)),
             align = c("l", "l", "r"), 
             graph.pos = 3,
             boxsize = 0.25,
             #clip = c(-0.65, 1.3),
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

tiff(paste0(out_dir, "fig_s4_040725.tiff"), units="in", width=11, height=3.5, res=300)
plot.new()
all
dev.off()

###############################################################################
#                            Fig S3 - Birthweight                             #
###############################################################################

### prepare data for plotting - bw ############################################

bw_df <- rbind(reg_bw, mr_bw, mvmr_bw) %>%
  mutate(`Beta (95% CI)` = str_replace(`Beta (95% CI)`, " - ", ", "))
bw_df$mean <- bw_df$b
bw_df$lower <- bw_df$lo_ci
bw_df$upper <- bw_df$up_ci

# add empty rows and order by outcome, within that approach
empty_rows <- data.frame(outcome_clean = c(unique(bw_df$outcome_clean)))
bw_df <- plyr::rbind.fill(empty_rows, bw_df) %>%
  arrange(match(outcome_clean, outcome_names_df$outcome_clean),
          match(approach, approach_order)) %>%
  # replace outcome_clean with "" except for row 2 where we want the label
  mutate(
    outcome_clean = case_when(approach == "MR IVW" ~ outcome_clean,
                              .default = ""))  %>%
  select(outcome_clean, approach, `Beta (95% CI)`, mean, lower, upper)
# headers
bw_df[1,] <- c("Outcome", "Approach", "Beta (95% CI)", NA, NA, NA)
# add additional empty row
#bw_df <- rbind(bw_df[1,], c(rep(NA, 5)), bw_df[2:5,])

# finalise text columns, add space for legend
alltext <- bw_df %>% select(outcome_clean, approach, `Beta (95% CI)`)
alltext <- alltext %>% replace(is.na(.), "")

# numeric results columns
resforalltext <- bw_df %>% select(mean, lower, upper)
resforalltext <- sapply(resforalltext, as.numeric)

### plot - bw #################################################################

styles <- fpShapesGp(
  box =
    rep(list(
      gpar(col = palette[1], fill = palette[1]), # for NA lines
      gpar(col = palette[1], fill = palette[1]),
      gpar(col = palette[2], fill = palette[2]),
      gpar(col = palette[3], fill = palette[3]),
      gpar(col = palette[4], fill = palette[4]),
      gpar(col = palette[5], fill = palette[5]),
      gpar(col = palette[6], fill = palette[6]),
      gpar(col = palette[7], fill = palette[7]),
      gpar(col = palette[8], fill = palette[8])
    ), 1), # number of outcomes
  lines = rep(list(
    gpar(col=palette[1]), # for NA lines
    gpar(col=palette[1]),
    gpar(col=palette[2]),
    gpar(col=palette[3]),
    gpar(col=palette[4]),
    gpar(col=palette[5]),
    gpar(col=palette[6]),
    gpar(col=palette[7]),
    gpar(col=palette[8])
  ), 1)
)

ticks <- c(-0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07)

all <- resforalltext %>%
  forestplot(labeltext = alltext,
             txt_gp = fpTxtGp(
               label = gpar(cex = .7),
               ticks = gpar(cex = .7),
               xlab = gpar(cex = .8)),
             is.summary=c(TRUE, rep(FALSE, 9)),
             align = c("l", "l", "r"), 
             graph.pos = 3,
             boxsize = 0.25,
             fn.ci_norm = fpDrawCircleCI,
             xlog = FALSE,
             xticks=ticks,
             xlab = "\nBeta per 1 year increase in age at menarche",
             zero=0,
             col = fpColors(zero="grey50"),
             colgap=unit(0.01, "npc"),
             shapes_gp = styles) |>
  fp_decorate_graph(grid =
                      structure(ticks[-c(1,length(ticks))],
                                gp = gpar(lty = 2, col = "grey50")))
print(all)

tiff(paste0(out_dir, "fig_s3_040725.tiff"), units="in", width=11, height=3.5, res=300)
plot.new()
all
dev.off()

q('no')