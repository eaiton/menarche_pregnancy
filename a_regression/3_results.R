###############################################################################
#                                                                             #
#                        Generate tables and plots of                         #
#                      multivariable regression results                       #
#                                                                             #
###############################################################################

###############################################################################
#                                   Set up                                    #
###############################################################################

# Clear the work environment
rm(list = ls())

# Required libraries
x <- c("forestplot", "dplyr", "tidyr", "grid", "paletteer", "stringr")
lapply(x, require, character.only = TRUE)

# Set directories
home_dir <- paste0(Sys.getenv("AGE_AT_MENARCHE_DIR"), "/working/")
setwd(home_dir)
# where dataset is saved
data_dir <- paste0(home_dir, "/results/a_regression/main/")
# to write plots to
out_dir <- paste0(home_dir, "/results/a_regression/main/forest_plots")

# Outcomes
outcome_names_df <- read.csv(file.path("data/b_mr/outcome/outcomes.csv"))
all_outcomes <- outcome_names_df$outcome
binary_outcomes <- all_outcomes[all_outcomes != "zbw_all"]
continuous_outcomes <- c("zbw_all")

# Adapt long ones for plots:
outcome_names_df[1, 2] <- "Any hypertensive\ndisorders of pregnancy"
outcome_names_df[6, 2] <- "Small for gestational age "
outcome_names_df[8, 2] <- "Gestational diabetes\nmellitus"
outcome_names_df[10, 2] <- "Large for gestational age "

# Read in results
linear_exp_binary_out <- read.csv(paste0(data_dir, "/linear_exp_binary_out.csv"))
linear_exp_cont_out <- read.csv(paste0(data_dir, "/linear_exp_cont_out.csv"))
cat_exp_binary_out <- read.csv(paste0(data_dir, "/cat_exp_binary_out.csv"))
cat_exp_cont_out <- read.csv(paste0(data_dir, "/cat_exp_cont_out.csv"))

###############################################################################
#                           Supplementary Fig 5                               #
###############################################################################

# Supplementary Figure 5. Categoric exposure plots - OR in figure
# binary outcomes - to illustrate sensitivity analysis

### Binary outcome
cat_obs_binary <- cat_exp_binary_out %>%
  # update outcome names
  select(-outcome_clean) %>%
  left_join(outcome_names_df) %>%
  # clean method names
  mutate(
    approach = case_when(
      model_name == "raw" ~ "Regression M1",
      model_name == "adjusted" ~ "Regression M2",
      model_name == "adjusted_bmi" ~ "Regression M3"),
  # ORs
      mean = estimate,
      lower = conf.low,
      upper = conf.high,
  # generate clean OR label for figure
      `OR (95% CI)` = paste0(format(round(mean, 2), nsmall = 2), " (",
                           format(round(lower, 2), nsmall = 2), ", ",
                           format(round(upper, 2), nsmall = 2), ")"),
  # generate new outcome_clean names
    outcome_clean_2 = case_when(
      term == "age_at_menarche_cat_relevel1" ~ paste0(outcome_clean, " - early"),
      term == "age_at_menarche_cat_relevel3" ~ paste0(outcome_clean, " - late")
    ))

# add empty rows for outcome headings
empty_rows <- data.frame(
  outcome_clean =
    rep(unique(cat_obs_binary$outcome_clean),
    each = 2), # each outcome repeated twice in a row, for early / late
  outcome_clean_2 =
    unique(cat_obs_binary$outcome_clean_2),
  #term = unique(cat_obs_binary$term),
  approach = "Placeholder")
# order by: outcome early/late, within that approach
cat_obs_binary <- plyr::rbind.fill(empty_rows, cat_obs_binary) %>%
  arrange(match(outcome_clean, outcome_names_df$outcome_clean),
    outcome_clean_2, approach)

### format for plot
binary_df <- cat_obs_binary %>%
# replace outcome_clean with "" except for row 2 where we want the label
  mutate(
    outcome_clean_2 = case_when(approach == "Regression M2" ~ outcome_clean_2,
    .default = ""))  %>%
  select(outcome_clean_2, `OR (95% CI)`, mean, lower, upper)
# headers
binary_df[1,] <- c("Outcome\n", "Odds ratio (95% CI)\n", NA, NA, NA)
# add additional empty row
#binary_df <- rbind(binary_df[1,], c(rep(NA, 5)), binary_df[2:75,])

# finalise text columns, add space for legend
alltext <- binary_df %>% select(outcome_clean_2, `OR (95% CI)`) %>%
  mutate(space = "                    ")
alltext <- alltext %>% replace(is.na(.), "")

# numeric results columns
resforalltext <- binary_df %>% select(mean, lower, upper)
resforalltext <- sapply(resforalltext, as.numeric)

palette <-
  c("#fcae12", "#e65d2f", "#a92e5e")
styles <- fpShapesGp(
  box =
    #list(gpar(col = palette[1], fill = palette[1])),
    rep(list(
      gpar(col = palette[1], fill = palette[1]),  # for NA lines
      gpar(col = palette[1], fill = palette[1]),
      gpar(col = palette[2], fill = palette[2]),
      gpar(col = palette[3], fill = palette[3])
    ), 24), # number of outcomes
  lines =
   # gpar(col=palette[1]),
    rep(list(
      gpar(col=palette[1]), # for NA lines
      gpar(col=palette[1]),
      gpar(col=palette[2]),
      gpar(col=palette[3])
    ), 24)
)

ticks <- c(0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)

all <- resforalltext %>%
  forestplot(labeltext = alltext,
             txt_gp = fpTxtGp(
              label = gpar(cex = .8),
              ticks = gpar(cex = .8),
              xlab = gpar(cex = .9)),
             is.summary=c(TRUE, rep(FALSE,71)),
             align = c("l", "l", "r"), 
             graph.pos = 3,
             boxsize = 0.4,
             fn.ci_norm = fpDrawCircleCI,
             xlog = TRUE,
             xticks=ticks,
             xlab = "\nOdds ratio relative to normative age at menarche",
             zero=1,
             col = fpColors(zero="grey50"),
             colgap=unit(0.01, "npc"),
             shapes_gp = styles) |>
  fp_decorate_graph(grid =
    structure(ticks[-c(1,length(ticks))],
    gp = gpar(lty = 2, col = "grey50")))
print(all)

tiff(paste0(out_dir, "/supp_fig5_040725.tiff"), units="in", width=11, height=12, res=300)
plot.new()
all
legend("right",
       #inset = 0.03,
       c("Regression M1",
         " ",
         "Regression M2",
         " ",
         "Regression M3"),
       title=as.expression(bquote(bold("Approach"))),
       border="grey50", box.lwd=1,
       cex = 0.8,
       col=c(palette[1], "white", palette[2],"white", palette[3]),
       pch=c(16))
dev.off()

###############################################################################
#                           Supplementary Fig 6                               #
###############################################################################

# Supplementary Figure 5. Categoric exposure plots - OR in figure
# cont outcomes (bw)

### Continuous outcome
cat_obs_bw <- cat_exp_cont_out %>%
  # update outcome names
  select(-outcome_clean) %>%
  left_join(outcome_names_df) %>%
  # clean method names
  mutate(
    approach = case_when(
      model_name == "raw" ~ "Regression M1",
      model_name == "adjusted" ~ "Regression M2",
      model_name == "adjusted_bmi" ~ "Regression M3"),
  # Betas
      mean = estimate,
      lower = conf.low,
      upper = conf.high,
  # generate clean Beta label for figure
      `Beta (95% CI)` = paste0(format(round(mean, 3), nsmall = 3), " (",
                           format(round(lower, 3), nsmall = 3), ", ",
                           format(round(upper, 3), nsmall = 3), ")"),
  # generate new outcome_clean names
    outcome_clean_2 = case_when(
      term == "age_at_menarche_cat_relevel1" ~ paste0(outcome_clean, " - early "),
      term == "age_at_menarche_cat_relevel3" ~ paste0(outcome_clean, " - late ")
    ))

# add empty rows for outcome headings
empty_rows <- data.frame(
  outcome_clean =
    rep(unique(cat_obs_bw$outcome_clean),
    each = 2), # each outcome repeated twice in a row, for early / late
  outcome_clean_2 =
    unique(cat_obs_bw$outcome_clean_2),
  #term = unique(cat_obs_bw$term),
  approach = "Placeholder")
# order by: outcome early/late, within that approach
cat_obs_bw <- plyr::rbind.fill(empty_rows, cat_obs_bw) %>%
  arrange(match(outcome_clean, outcome_names_df$outcome_clean),
    outcome_clean_2, approach)

### format for plot
bw_df <- cat_obs_bw %>%
# replace outcome_clean with "" except for row 2 where we want the label
  mutate(
    outcome_clean_2 = case_when(approach == "Regression" ~ outcome_clean_2,
    .default = ""))  %>%
  select(outcome_clean_2, `Beta (95% CI)`, mean, lower, upper)
# headers
bw_df[1,] <- c("Outcome\n", "Beta (95% CI)\n", NA, NA, NA)
# add additional empty row
#binary_df <- rbind(binary_df[1,], c(rep(NA, 5)), binary_df[2:75,])

# finalise text columns, add space for legend
alltext <- bw_df %>% select(outcome_clean_2, `Beta (95% CI)`) %>%
  mutate(space = "                    ")
alltext <- alltext %>% replace(is.na(.), "")

# numeric results columns
resforalltext <- bw_df %>% select(mean, lower, upper)
resforalltext <- sapply(resforalltext, as.numeric)

palette <-
  c("#fcae12", "#e65d2f", "#a92e5e")
styles <- fpShapesGp(
  box =
    #list(gpar(col = palette[1], fill = palette[1])),
    rep(list(
      gpar(col = palette[1], fill = palette[1]),  # for NA lines
      gpar(col = palette[1], fill = palette[1]),
      gpar(col = palette[2], fill = palette[2]),
      gpar(col = palette[3], fill = palette[3])
    ), 2), # number of outcomes
  lines =
   # gpar(col=palette[1]),
    rep(list(
      gpar(col=palette[1]), # for NA lines
      gpar(col=palette[1]),
      gpar(col=palette[2]),
      gpar(col=palette[3])
    ), 2)
)

ticks <- c(-0.15, -0.10, -0.05, 0, 0.05, 0.10)

all <- resforalltext %>%
  forestplot(labeltext = alltext,
             txt_gp = fpTxtGp(
              label = gpar(cex = .7),
              ticks = gpar(cex = .7),
              xlab = gpar(cex = .8)),
             is.summary=c(TRUE, rep(FALSE, 10)),
             align = c("l", "l", "r"), 
             graph.pos = 3,
             boxsize = 0.15,
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

tiff(paste0(out_dir, "/supp_fig6_040725.tiff"), units="in", width=11, height=3.8, res=300)
plot.new()
all
legend("right",
       c("Regression M1",
         " ",
         "Regression M2",
         " ",
         "Regression M3"),
       title=as.expression(bquote(bold("Approach"))),
       border="grey50", box.lwd=1,
       cex = 0.8,
       col=c(palette[1], "white", palette[2],"white", palette[3]),
       pch=c(16))
dev.off()

q('no')