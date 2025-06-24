suppressPackageStartupMessages({
  library(mlr3proba)
  library(cowplot)
})

osipov_abbrev = "Osipov et al. (2024)"
wissel_abbrev = "Wissel et al. (2023)"
cao_abbrev    = "Cao et al. (2021)"

task_wissel2023 = readRDS(file.path("data/wissel2023/task_list.rds"))$clinical
task_osipov2024 = readRDS(file.path("data/osipov2024/task_list.rds"))$clinical
task_cao2021 = readRDS(file.path("data/cao2021/task_list.rds"))$clinical

p1 = autoplot(task_wissel2023, theme = theme_minimal(base_size = 14)) +
  ylim(c(0,1)) +
  labs(x = "Months", y = "Survival probability", title = wissel_abbrev) +
  theme(text = element_text(family = "Arial"))

p2 = autoplot(task_osipov2024, theme = theme_minimal(base_size = 14)) +
  ylim(c(0,1)) +
  labs(x = "Months", y = "Survival probability", title = osipov_abbrev) +
  theme(text = element_text(family = "Arial"))

p3 = autoplot(task_cao2021, theme = theme_minimal(base_size = 14)) +
  ylim(c(0,1)) +
  labs(x = "Months", y = "Survival probability", title = cao_abbrev) +
  theme(text = element_text(family = "Arial"))

plot_grid(p1, p2, p3, nrow = 1, labels = c("a", "b", "c"))
ggsave("data/img/km_curves.png", width = 10, height = 4, dpi = 600, bg = "white")
