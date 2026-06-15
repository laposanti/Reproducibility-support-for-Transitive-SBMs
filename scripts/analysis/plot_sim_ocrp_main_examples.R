#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

input_csv <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "output/simulation/raw/DemoOCRPvar_runs/DemoOCRPvar_run_20260609_131104/full_simulation_crossfit_final_DemoOCRPvar_run_20260609_131104.csv"
}

output_pdf <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  "output/paper/figures/current/sim_ocrp_main_examples.pdf"
}

sim <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)

osbm <- sim %>%
  filter(gen_model %in% c("WST", "SST"), hierch %in% c("hard", "easy")) %>%
  mutate(
    rep_id = as.integer(rep_id),
    fit_model = factor(fit_model, levels = c("WST", "SST", "DC-SBM")),
    gen_model = factor(gen_model, levels = c("WST", "SST")),
    hierch = factor(hierch, levels = c("hard", "easy"))
  )

dataset_wide <- osbm %>%
  select(rep_id, gen_model, hierch, fit_model, looic) %>%
  tidyr::pivot_wider(names_from = fit_model, values_from = looic) %>%
  mutate(
    best_ordered_model = ifelse(WST <= SST, "WST", "SST"),
    best_ordered_looic = pmin(WST, SST),
    dc_gain = `DC-SBM` - best_ordered_looic,
    scenario_code = paste0(substr(as.character(gen_model), 1, 1),
                           ifelse(hierch == "hard", "H", "E"),
                           rep_id)
  ) %>%
  arrange(gen_model, hierch, rep_id) %>%
  mutate(dataset_idx = seq_len(n()))

pick_example <- function(gen_name) {
  cand <- dataset_wide %>%
    filter(gen_model == gen_name, hierch == "easy", best_ordered_model == gen_name) %>%
    arrange(desc(dc_gain), rep_id)
  if (!nrow(cand)) {
    cand <- dataset_wide %>%
      filter(gen_model == gen_name, best_ordered_model == gen_name) %>%
      arrange(desc(dc_gain), rep_id)
  }
  cand[1, , drop = FALSE]
}

example_wst <- pick_example("WST")
example_sst <- pick_example("SST")
example_keys <- bind_rows(example_wst, example_sst) %>%
  mutate(
    example_tag = c("WST example", "SST example")
  )

left_df <- dataset_wide %>%
  left_join(
    example_keys %>% select(rep_id, gen_model, hierch, example_tag),
    by = c("rep_id", "gen_model", "hierch")
  ) %>%
  mutate(
    selected = !is.na(example_tag),
    point_col = ifelse(best_ordered_model == "WST", "WST", "SST")
  )

example_long <- osbm %>%
  inner_join(
    example_keys %>% select(rep_id, gen_model, hierch, example_tag),
    by = c("rep_id", "gen_model", "hierch")
  ) %>%
  mutate(
    example_tag = factor(example_tag, levels = c("WST example", "SST example")),
    fit_model = factor(fit_model, levels = c("WST", "SST", "DC-SBM"))
  )

cols <- c("WST" = "#1f77b4", "SST" = "#d62728", "DC-SBM" = "#4d4d4d")

p_left <- ggplot(left_df, aes(x = dataset_idx, y = dc_gain)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey60") +
  geom_point(aes(color = point_col), size = 2.4) +
  geom_point(
    data = subset(left_df, selected),
    shape = 21, stroke = 0.8, fill = "white", color = "black", size = 3.6
  ) +
  geom_text(
    data = subset(left_df, selected),
    aes(label = example_tag),
    nudge_y = 5, size = 3
  ) +
  scale_color_manual(values = cols, breaks = c("WST", "SST")) +
  scale_x_continuous(
    breaks = left_df$dataset_idx,
    labels = left_df$scenario_code
  ) +
  labs(
    x = NULL,
    y = expression(LOOIC[DC] - min(LOOIC[WST], LOOIC[SST])),
    color = "Best ordered"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

p_right <- ggplot(
  example_long,
  aes(x = fit_model, y = delta_looic_best, color = fit_model)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey60") +
  geom_pointrange(
    aes(
      ymin = pmax(delta_looic_best - delta_looic_best_se, 0),
      ymax = delta_looic_best + delta_looic_best_se
    ),
    linewidth = 0.5
  ) +
  facet_wrap(~ example_tag, nrow = 1) +
  scale_color_manual(values = cols) +
  labs(
    x = NULL,
    y = expression(Delta * LOOIC),
    color = "Model"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
pdf(output_pdf, width = 11.5, height = 4.8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1.5, 1), "null"))))
print(p_left, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_right, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

message("Wrote figure: ", output_pdf)
