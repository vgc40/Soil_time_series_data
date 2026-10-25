# =============================================================================
# 21 Checks that explain the Figure 2 season tests (plain-language guide)
# =============================================================================
#
# These checks support Docs/Plain_language_guide_v4.docx. They do not change
# any result; they show why the season tests are done the way they are.
#
#   1. Sample counts per date in the original and corrected (v3) MicroResp
#      files: the January duplicate and the missing February plot.
#   2. The same whole-profile PERMANOVA run five ways (v3 method, date means,
#      all samples with whole dates shuffled, built-in vegan with whole dates
#      shuffled, and the mix-only version).
#   3. Made-up seasons: how often each method calls a meaningless grouping of
#      dates "significant".
#   4. Bray-Curtis distances: how big a uniform 34% increase is compared with
#      ordinary plot-to-plot differences.

rm(list = ls())
source(file.path("Scripts", "helpers.R"))
set.seed(20231023)

out_dir <- file.path("Results", "manuscript_v4", "Figure2", "explainers")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Sample counts per date: original vs corrected MicroResp release
# -----------------------------------------------------------------------------

count_samples <- function(path, release) {
  read_csv(path, col_types = cols(.default = col_character()), na = c("", "NA", "na")) %>%
    filter(Irrigation_trmt_pct == "100", Substrate != "Water") %>%
    distinct(Sample_cat, Collection_date, Plot) %>%
    mutate(Sampling_Date = mdy(Collection_date)) %>%
    group_by(Sampling_Date) %>%
    summarise(n_samples = n(), plots_listed = paste(sort(as.integer(Plot)), collapse = ", "),
              .groups = "drop") %>%
    mutate(release = release)
}

counts <- bind_rows(
  count_samples(file.path("DataHub March", "Soil_Phenotypes",
                          "SoilMicrobiomeSFA_prosser-time-series_MicroResp.csv"),
                "Original release"),
  count_samples(input_files$microresp, "Corrected release (v3)")
)
write_csv(counts, file.path(out_dir, "microresp_samples_per_date_by_release.csv"))

# -----------------------------------------------------------------------------
# 2. One whole-profile PERMANOVA, five ways
# -----------------------------------------------------------------------------

sir <- load_microresp_100() %>%
  select(Sample_cat, Sampling_Date, Plot_ID, Season, Substrate, SIR) %>%
  pivot_wider(names_from = Substrate, values_from = SIR) %>%
  arrange(Sampling_Date, Plot_ID)
community <- as.matrix(select(sir, -c(Sample_cat, Sampling_Date, Plot_ID, Season)))
dist_samples <- vegdist(community)

dates <- sort(unique(sir$Sampling_Date))
date_means <- rowsum(community, as.character(sir$Sampling_Date)) /
  as.vector(table(as.character(sir$Sampling_Date)))
date_means <- date_means[as.character(dates), ]
date_season <- assign_season(dates)

# (a) v3 method: samples as replicates, labels shuffled within each plot
v3 <- adonis2(dist_samples ~ Season, data = sir,
              permutations = how(nperm = 9999, blocks = sir$Plot_ID))

# (b) date means, ordinary adonis2 with dates shuffled
# Exact version: every choice of 5 dormant dates out of 13 (as in Figure 2).
date_mean_fit <- filter(date_level_season_test(community, sir$Sampling_Date), level == "Date means (primary)")

# (c) all samples, season labels shuffled among whole dates (every choice)
all_choices <- combn(length(dates), sum(date_season == "Plant-dormant"), simplify = FALSE)
f_for <- function(groups) suppressMessages(adonis2(dist_samples ~ groups, permutations = 0))$F[1]
real_f <- f_for(sir$Season)
fake_f <- map_dbl(all_choices, function(k) f_for(factor(sir$Sampling_Date %in% dates[k])))
whole_date_p <- mean(fake_f >= real_f - 1e-12)

# (d) built-in vegan whole-date shuffling; needs equal plots per date, so the
#     February 2024 date (5 plots) is dropped
balanced <- sir$Sampling_Date != as.Date("2024-02-12")
builtin <- adonis2(vegdist(community[balanced, ]) ~ Season, data = sir[balanced, ],
                   permutations = how(within = Within(type = "none"),
                                      plots = Plots(strata = factor(sir$Sampling_Date[balanced]),
                                                    type = "free"),
                                      nperm = 9999))

# (e) mix only: each sample scaled to its own total
mix <- date_level_season_test(community / rowSums(community), sir$Sampling_Date)

permanova_versions <- tribble(
  ~version, ~what_counts_as_a_replicate, ~n_rows, ~R2, ~p,
  "v3 method: plots shuffled one by one", "Each sample (plot-date)", nrow(sir), v3$R2[1], v3$`Pr(>F)`[1],
  "All samples, whole dates shuffled", "Each date", nrow(sir),
  suppressMessages(adonis2(dist_samples ~ Season, data = sir, permutations = 0))$R2[1], whole_date_p,
  "Date averages", "Each date", length(dates), date_mean_fit$R2, date_mean_fit$p_permanova,
  "Built-in vegan, whole dates shuffled, February dropped", "Each date", sum(balanced),
  builtin$R2[1], builtin$`Pr(>F)`[1],
  "Mix only (date averages, profiles scaled to their total)", "Each date", length(dates),
  mix$R2[mix$level == "Date means (primary)"], mix$p_permanova[mix$level == "Date means (primary)"]
)
write_csv(permanova_versions, file.path(out_dir, "permanova_five_ways.csv"))

# -----------------------------------------------------------------------------
# 3. Made-up seasons
# -----------------------------------------------------------------------------
# Pick 5 of the 13 dates at random and call them "dormant". The grouping is
# meaningless, so a trustworthy test should call it significant about 5% of
# the time.

date_mean_dist <- vegdist(date_means)
f_dates <- map_dbl(all_choices, function(k) {
  suppressMessages(adonis2(date_mean_dist ~ factor(seq_along(dates) %in% k), permutations = 0))$F[1]
})

n_fake <- 200
picks <- sample(length(all_choices), n_fake)
fake <- map_dfr(picks, function(pick) {
  k <- all_choices[[pick]]
  groups <- factor(sir$Sampling_Date %in% dates[k])
  tibble(
    fake_dormant_dates = paste(format(dates[k], "%b %Y"), collapse = "; "),
    p_v3_method = adonis2(dist_samples ~ groups,
                          permutations = how(nperm = 499, blocks = sir$Plot_ID))$`Pr(>F)`[1],
    p_date_method = mean(f_dates >= f_dates[pick] - 1e-12)
  )
})
fake_summary <- tibble(
  n_fake_seasons = n_fake,
  pct_significant_v3_method = 100 * mean(fake$p_v3_method < 0.05),
  pct_significant_date_method = 100 * mean(fake$p_date_method < 0.05),
  expected_pct_for_a_fair_test = 5
)
write_csv(fake, file.path(out_dir, "fake_seasons_each_test.csv"))
write_csv(fake_summary, file.path(out_dir, "fake_seasons_summary.csv"))

# -----------------------------------------------------------------------------
# 4. How big is a 34% increase on the NMDS distance scale?
# -----------------------------------------------------------------------------

level <- read_csv(file.path("Results", "manuscript_v4", "Figure2",
                            "Figure2A_overall_level_shift.csv"), show_col_types = FALSE)
multiplier <- 2^level$estimate
d <- as.matrix(dist_samples)
same_date <- outer(sir$Sampling_Date, sir$Sampling_Date, "==")
off_diagonal <- upper.tri(d)
distances <- tibble(
  comparison = c(sprintf("A sample vs itself with every substrate x%.2f", multiplier),
                 "Two plots sampled on the same date (median)",
                 "Two samples from different dates (median)"),
  bray_curtis = c(
    mean(map_dbl(seq_len(nrow(community)), function(r) {
      as.numeric(vegdist(rbind(community[r, ], community[r, ] * multiplier)))
    })),
    median(d[same_date & off_diagonal]),
    median(d[!same_date & off_diagonal])
  )
)
write_csv(distances, file.path(out_dir, "bray_curtis_scale_of_seasonal_shift.csv"))

print(counts %>% filter(Sampling_Date %in% as.Date(c("2024-01-02", "2024-02-12"))))
print(permanova_versions)
print(fake_summary)
print(distances)
