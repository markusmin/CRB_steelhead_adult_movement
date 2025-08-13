### 07-09: Figure 8 Final fates and covariates

# This script takes the output from the final fates and covariates
# simulations and plots the output

# First, need to load in all of the model runs and all of the packages.
source("R/07_model_analysis/07-01_load_stan_models.R")
source("R/07_model_analysis/transport/07-01_load_stan_models_transport.R")

# Load all of the files

# Load the files
load("figures/final_fates_scenarios/simulation_runs/UMA_MCN_winterspill_homing.rda")
load("figures/final_fates_scenarios/simulation_runs/ENT_WEL_winterspill_homing.rda")
load("figures/final_fates_scenarios/simulation_runs/JDR_MCN_winterspill_homing.rda")
load("figures/final_fates_scenarios/simulation_runs/WAWA_ICH_winterspill_homing.rda")
load("figures/final_fates_scenarios/simulation_runs/WEN_RRE_winterspill_homing.rda")
load("figures/final_fates_scenarios/simulation_runs/YAK_PRA_winterspill_homing.rda")

# change names of these outputs
load("figures/transport/final_fates_scenarios/simulation_runs/TUC_LGR_not_transported_homing.rda")
TUC_not_transported_LGR_winterspill_homing <- TUC_LGR_winterspill_homing
load("figures/transport/final_fates_scenarios/simulation_runs/TUC_LGR_transported_homing.rda")
TUC_transported_LGR_winterspill_homing <- TUC_LGR_winterspill_homing

#### Figure 7: Final fates vs. covariates ####



plot_ff_cov <- function(data, dam_name, dam_spill_column){
  rear_colors <- c(hatchery = "#ff7f00", natural = "#33a02c")
  rear_shapes <- c(hatchery = 21, natural = 21)
  condition_colors <- c("coldest" = "#2c7bb6", "average" = "goldenrod2", "warmest" = "#d7191c")
  
  
  # fix the condition levels
  data %>% 
    mutate(conditions = ifelse(conditions == "cool", "coldest",
                               ifelse(conditions == "warm", "warmest", conditions))) -> data
  
  # fix the rear levels
  data %>% 
    mutate(rear_type = ifelse(rear_type == "wild", "natural", rear_type)) -> data
  
  # add a condition x rear
  data %>% 
    mutate(rear_condition = paste(rear_type, conditions)) %>% 
    mutate(rear_condition = factor(rear_condition, levels = c("natural coldest",
                                                              "natural average",
                                                              "natural warmest",
                                                              "hatchery coldest",
                                                              "hatchery average", 
                                                              "hatchery warmest"))) -> data
  
  rear_condition_fills <- c("natural coldest" =  "#2c7bb6",  "hatchery coldest" = "white",
                            "natural average" = "goldenrod2",  "hatchery average" = "white",
                            "natural warmest" = "#d7191c", "hatchery warmest" = "white")
  
  data$conditions <- factor(data$conditions, levels = c("coldest", "average", "warmest"))
  
  plot <- ggplot(data, aes(x = eval(parse(text = dam_spill_column)), y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = conditions,
                           fill = rear_condition, shape = rear_type)) +
    geom_linerange(aes(group = rear_condition), position=position_dodge(width=3)) +
    geom_point(aes(size = rear_type, group = rear_condition), position=position_dodge(width=3)) +
    ylab("Homing Probability") +
    xlab(paste0("Days of Winter Spill at ", dam_name)) +
    # ggtitle(" ") +
    # Create a scale common to all
    scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(lim = c(-3, 95), breaks = c(0, 30, 60, 90)) +
    scale_fill_manual(values = rear_condition_fills, guide = "none") +
    scale_size_manual(values = c("hatchery" = 4, "natural" = 3.5), guide = "none") +
    scale_shape_manual(values = rear_shapes) +
    scale_color_manual(values = condition_colors) +
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 15),
          panel.grid.major = element_line(color = "gray90"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
          legend.position = "none",
          plot.margin = unit(c(1.2, 0.2, 0.2, 0.2),"cm"))
  
  return(plot)
}


# John Day
JDR_FF_cov_plot <- plot_ff_cov(data = JDR_MCN_winterspill_homing, dam_name = "McNary Dam", dam_spill_column = "MCN_winterspill_actual")
# Umatilla
UMA_FF_cov_plot <- plot_ff_cov(data = UMA_MCN_winterspill_homing, dam_name = "McNary Dam", dam_spill_column = "MCN_winterspill_actual")
# Walla Walla
WAWA_FF_cov_plot <- plot_ff_cov(data = WAWA_ICH_winterspill_homing, dam_name = "Ice Harbor Dam", dam_spill_column = "ICH_winterspill_actual")
# Wenatchee
WEN_FF_cov_plot <- plot_ff_cov(data = WEN_RRE_winterspill_homing, dam_name = "Rocky Reach Dam", dam_spill_column = "RRE_winterspill_actual")
# Entiat
ENT_FF_cov_plot <- plot_ff_cov(data = ENT_WEL_winterspill_homing, dam_name = "Wells Dam", dam_spill_column = "WEL_winterspill_actual")
# Yakima
YAK_FF_cov_plot <- plot_ff_cov(data = YAK_PRA_winterspill_homing, dam_name = "Priest Rapids Dam", dam_spill_column = "PRA_winterspill_actual")

# Tucannon - transported vs. not transported
TUC_transported_FF_cov_plot <- plot_ff_cov(data = TUC_transported_LGR_winterspill_homing, dam_name = "Lower Granite Dam", dam_spill_column = "LGR_winterspill_actual")
TUC_not_transported_FF_cov_plot <- plot_ff_cov(data = TUC_not_transported_LGR_winterspill_homing, dam_name = "Lower Granite Dam", dam_spill_column = "LGR_winterspill_actual")


# Create the legend figure by itself
rear_colors <- c(hatchery = "#ff7f00", natural = "#33a02c")
rear_shapes <- c(hatchery = 21, natural = 21)
condition_colors <- c("coldest" = "#2c7bb6", "average" = "goldenrod2", "warmest" = "#d7191c")


# fix the condition levels
WAWA_ICH_winterspill_homing %>% 
  mutate(conditions = ifelse(conditions == "cool", "coldest",
                             ifelse(conditions == "warm", "warmest", conditions))) -> WAWA_ICH_winterspill_homing

# fix the rear levels
WAWA_ICH_winterspill_homing %>% 
  mutate(rear_type = ifelse(rear_type == "wild", "natural", rear_type)) -> WAWA_ICH_winterspill_homing

# add a condition x rear
WAWA_ICH_winterspill_homing %>% 
  mutate(rear_condition = paste(rear_type, conditions)) -> WAWA_ICH_winterspill_homing

rear_condition_fills <- c("natural coldest" =  "#2c7bb6",  "hatchery coldest" = "white",
                          "natural average" = "goldenrod2",  "hatchery average" = "white",
                          "natural warmest" = "#d7191c", "hatchery warmest" = "white")

WAWA_ICH_winterspill_homing$conditions <- factor(WAWA_ICH_winterspill_homing$conditions, levels = c("coldest", "average", "warmest"))

plot_for_legend <- ggplot(WAWA_ICH_winterspill_homing, aes(x = ICH_winterspill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = conditions,
                                                           fill = rear_condition, shape = rear_type)) +
  geom_point(size = 8, position=position_dodge(width=3)) +
  geom_linerange(position=position_dodge(width=3), show.legend = FALSE) +
  ylab("Homing Probability") +
  xlab(paste0("Days of Winter Spill at Ice Harbor Dam")) +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_fill_manual(values = rear_condition_fills, guide = "none") +
  scale_size_manual(values = c("hatchery" = 10, "natural" = 8), guide = "none") +
  scale_shape_manual(values = rear_shapes) +
  # scale_shape_manual(values = c(hatchery = 24, natural = 19)) +
  scale_color_manual(values = condition_colors) +
  guides(shape=guide_legend(title="Rearing Type", override.aes = list(fill = c("white", "black"),
                                                                      size = c(10, 10))),
         color = guide_legend(title = "Basin Conditions")) +
  # theme(plot.title = element_text(size = 12),
  #       # axis.text.y = element_text(color = rev(state_significance_colors)),
  #       axis.title = element_text(size = 12),
  #       axis.text.x = element_text(size = 12))
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(1.25, "cm"),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 20),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))

FF_cov_legend <- ggpubr::get_legend(plot_for_legend)
FF_cov_legend_gg <- as_ggplot(FF_cov_legend) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

FF_cov_combined_plot <- cowplot::ggdraw(ggarrange(JDR_FF_cov_plot, UMA_FF_cov_plot, WAWA_FF_cov_plot,
                                                  WEN_FF_cov_plot, TUC_transported_FF_cov_plot, TUC_not_transported_FF_cov_plot, ENT_FF_cov_plot,
                                                  YAK_FF_cov_plot, FF_cov_legend_gg, nrow = 3, ncol = 3,
                                                  # labels = c("(A) John Day River", " (B) Umatilla River", "(C) Walla Walla River", 
                                                  #            "(D) Wenatchee River", "(E) Tucannon River (transported)", "(F) Tucannon River (in-river)",
                                                  #            "(G) Entiat River", "(H) Yakima River"),
                                                  labels = c("(A) JDR", " (B) UMA", "(C) WAWA", 
                                                             "(D) WEN", "(E) TUC (transported)", "(F) TUC (in-river)",
                                                             "(G) ENT", "(H) YAK"),
                                                  label.x = 0.05, label.y = 0.95,  font.label = list(size = 18, face = "plain"),
                                                  hjust = 0, vjust = 0)) + theme(plot.background = element_rect(fill="white", color = NA))


ggsave(here::here("figures", "transport", "paper_figures", "Fig8_FF_cov_combined_plot_transport.png"), FF_cov_combined_plot, height = 16, width = 16)
