# 07-02_Stan_model_diagnostics

# This script takes the output from the stan model runs in the /Stan/ folder and
# runs diagnostic checks

# First, need to load in all of the model runs and all of the packages.
source("R/07_model_analysis/07-01_load_stan_models.R")
source("R/07_model_analysis/alternative_spill_treatment/07-01_load_stan_models_transport.R")

##### Run diagnostic summaries: rhat and ess #####

# Create six-panel diagnostic plots for rhat, ess_bulk, and ess_tail

create_rhat_hist <- function(fit_summary, population){
  rhat_plot <- ggplot(fit_summary, aes(x = rhat)) +
    geom_histogram() + 
    # annotate("text", label = "Range:", x = max(fit_summary$rhat, na.rm = T), y = 600, hjust = 1) +
    # annotate("text", label = round(range(fit_summary$rhat, na.rm = T),3)[1], x = max(fit_summary$rhat, na.rm = T), y = 550,hjust = 1) +
    # annotate("text", label = paste0(" - ", round(range(fit_summary$rhat, na.rm = T),3)[2]), x = max(fit_summary$rhat, na.rm = T), y = 500, hjust = 1) +
    ylab("N parameters") +
    ggtitle(population)
  
  return(rhat_plot)
}

create_ess_bulk_hist <- function(fit_summary, population){
  ess_bulk_plot <- ggplot(fit_summary, aes(x = ess_bulk)) +
    geom_histogram() + 
    # annotate("text", label = "Range:", x = max(fit_summary$ess_bulk, na.rm = T), y = 600, hjust = 1) +
    # annotate("text", label = round(range(fit_summary$ess_bulk, na.rm = T),3)[1], x = max(fit_summary$ess_bulk, na.rm = T), y = 550,hjust = 1) +
    # annotate("text", label = paste0(" - ", round(range(fit_summary$ess_bulk, na.rm = T),3)[2]), x = max(fit_summary$ess_bulk, na.rm = T), y = 500, hjust = 1) +
    ylab("N parameters") +
    ggtitle(population)
  
  return(ess_bulk_plot)
}

create_ess_tail_hist <- function(fit_summary, population){
  ess_tail_plot <- ggplot(fit_summary, aes(x = ess_tail)) +
    geom_histogram() + 
    # annotate("text", label = "Range:", x = max(fit_summary$ess_tail, na.rm = T), y = 600, hjust = 1) +
    # annotate("text", label = round(range(fit_summary$ess_tail, na.rm = T),3)[1], x = max(fit_summary$ess_tail, na.rm = T), y = 550,hjust = 1) +
    # annotate("text", label = paste0(" - ", round(range(fit_summary$ess_tail, na.rm = T),3)[2]), x = max(fit_summary$ess_tail, na.rm = T), y = 500, hjust = 1) +
    ylab("N parameters") +
    ggtitle(population)
  
  return(ess_tail_plot)
}

# Run functions for all populations
UCW_rhat_plot <- create_rhat_hist(fit_summary = UCW_fit_summary, population = "Upper Columbia, Natural")
UCH_rhat_plot <- create_rhat_hist(fit_summary = UCH_fit_summary, population = "Upper Columbia, Hatchery")
MCW_rhat_plot <- create_rhat_hist(fit_summary = MCW_fit_summary, population = "Middle Columbia, Natural")
MCH_rhat_plot <- create_rhat_hist(fit_summary = MCH_fit_summary, population = "Middle Columbia, Hatchery")
SRW_T_rhat_plot <- create_rhat_hist(fit_summary = SRW_T_fit_summary, population = "Snake River, Natural, transported")
SRW_NT_rhat_plot <- create_rhat_hist(fit_summary = SRW_NT_fit_summary, population = "Snake River, Natural, in-river")
SRH_T_rhat_plot <- create_rhat_hist(fit_summary = SRH_T_fit_summary, population = "Snake River, Hatchery, transported")
SRH_NT_rhat_plot <- create_rhat_hist(fit_summary = SRH_NT_fit_summary, population = "Snake River, Hatchery, in-river")

UCW_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = UCW_fit_summary, population = "Upper Columbia, Natural")
UCH_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = UCH_fit_summary, population = "Upper Columbia, Hatchery")
MCW_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = MCW_fit_summary, population = "Middle Columbia, Natural")
MCH_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = MCH_fit_summary, population = "Middle Columbia, Hatchery")
SRW_T_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = SRW_T_fit_summary, population = "Snake River, Natural, transported")
SRW_NT_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = SRW_NT_fit_summary, population = "Snake River, Natural, in-river")
SRH_T_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = SRH_T_fit_summary, population = "Snake River, Hatchery, transported")
SRH_NT_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = SRH_NT_fit_summary, population = "Snake River, Hatchery, in-river")

UCW_ess_tail_plot <- create_ess_tail_hist(fit_summary = UCW_fit_summary, population = "Upper Columbia, Natural")
UCH_ess_tail_plot <- create_ess_tail_hist(fit_summary = UCH_fit_summary, population = "Upper Columbia, Hatchery")
MCW_ess_tail_plot <- create_ess_tail_hist(fit_summary = MCW_fit_summary, population = "Middle Columbia, Natural")
MCH_ess_tail_plot <- create_ess_tail_hist(fit_summary = MCH_fit_summary, population = "Middle Columbia, Hatchery")
SRW_T_ess_tail_plot <- create_ess_tail_hist(fit_summary = SRW_T_fit_summary, population = "Snake River, Natural, transported")
SRW_NT_ess_tail_plot <- create_ess_tail_hist(fit_summary = SRW_NT_fit_summary, population = "Snake River, Natural, in-river")
SRH_T_ess_tail_plot <- create_ess_tail_hist(fit_summary = SRH_T_fit_summary, population = "Snake River, Hatchery, transported")
SRH_NT_ess_tail_plot <- create_ess_tail_hist(fit_summary = SRH_NT_fit_summary, population = "Snake River, Hatchery, in-river")

# arrange them
rhat_comb_plots <- ggarrange(plotlist = list(UCW_rhat_plot, UCH_rhat_plot, MCW_rhat_plot, MCH_rhat_plot, 
                                             SRW_T_rhat_plot, SRW_NT_rhat_plot, SRH_T_rhat_plot, SRH_NT_rhat_plot),
                             ncol = 4, nrow = 2)

ggsave("figures/transport/diagnostic_plots/rhat_comb_plots.png", rhat_comb_plots, height = 8, width = 14)

ess_bulk_comb_plots <- ggarrange(plotlist = list(UCW_ess_bulk_plot, UCH_ess_bulk_plot, MCW_ess_bulk_plot, MCH_ess_bulk_plot,
                                                 SRW_T_ess_bulk_plot, SRW_NT_ess_bulk_plot, SRH_T_ess_bulk_plot, SRH_NT_ess_bulk_plot),
                             ncol = 4, nrow = 2)

ggsave("figures/transport/diagnostic_plots/ess_bulk_comb_plots.png", ess_bulk_comb_plots, height = 8, width = 14)

ess_tail_comb_plots <- ggarrange(plotlist = list(UCW_ess_tail_plot, UCH_ess_tail_plot, MCW_ess_tail_plot, MCH_ess_tail_plot,
                                                 SRW_T_ess_tail_plot, SRW_NT_ess_tail_plot, SRH_T_ess_tail_plot, SRH_NT_ess_tail_plot),
                             ncol = 4, nrow = 2)

ggsave("figures/transport/diagnostic_plots/ess_tail_comb_plots.png", ess_tail_comb_plots, height = 8, width = 14)




