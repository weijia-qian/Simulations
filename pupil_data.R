# plot pupil light response curves
load(file = "pupil_data.rda")

sgrid1 <- unique(df$time)[-120] # remove the last frame

df_long <- df_post1_r %>%
  pivot_longer(
    cols = starts_with("pct"),
    names_to = "time",
    values_to = "pct_chg"
  ) %>%
  mutate(time = rep(sgrid1, times = nrow(df_post1_r)),
         use = factor(is_user, labels = c("Didn't smoke", "Smoked")))

# Sample 10 individuals
set.seed(916)  # for reproducibility
sample_ids <- df_post1_r %>%
  # group_by(is_user) %>%
  slice_sample(n = 10) %>%
  pull(ptid)

# Filter and reshape for plotting
df_long_sample <- df_long %>%
  filter(ptid %in% c("003-078", "003-114")) 

# Plot
p_pupil <- ggplot(df_long) +
  geom_line(aes(x = time, y = pct_chg, group = ptid), color = "lightgray") +
  geom_line(data = df_long_sample, aes(x = time, y = pct_chg, group = ptid, color = use), size = 0.8) +
  labs(
    x = "Time since light stimulus (seconds)",
    y = "Percent change in pupil diameter",
    color = ""
  ) +
  scale_color_manual(
    values = c("Didn't smoke" = "#FF7043", "Smoked" = "#2E7D32")
  )

# line plot for survival
plot_survival <- df_post1_r %>%
  filter(ptid %in% sample_ids) %>%
  mutate(
    id = as.integer(factor(ptid)),
    event_time = ifelse(time_since_use == 480, 75, time_since_use),
    line_type = ifelse(ptid %in% c("003-078", "003-114"), "dashed", "solid")
  )

# Plot with thicker lines, filled shapes (circle for death, triangle for censored)
p_survival <- ggplot(data = plot_survival, aes(x = event_time, y = id)) +
  geom_segment(
    aes(x = 0, xend = event_time, yend = id, linetype = line_type),
    color = "grey55",
    size = 1.2
  ) +  # Thicker lines
  geom_point(
    aes(color = factor(is_user), shape = factor(is_user)),
    size = 3
  ) +  # End markers
  scale_linetype_identity() +  # Use linetype directly from variable
  scale_color_manual(
    values = c("#FF7043", "#2E7D32"),
    labels = c("Didn't smoke", "Smoked")
  ) +
  scale_shape_manual(
    values = c(17, 16),
    labels = c("Didn't smoke", "Smoked")
  ) +
  labs(
    x = "Time since smoking (minutes)",
    y = "Individual",
    color = "",
    shape = ""
  ) +
  theme(
    axis.text.y = element_blank()
  )

p1 <- p_pupil + p_survival
# ggsave(file.path(fig_dir, "fig_curves_avg.jpeg"), p1, width = 8, height = 4, dpi = 500)