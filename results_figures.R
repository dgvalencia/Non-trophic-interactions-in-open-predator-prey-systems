#======================================================================#
#  Script to Reproduce Figures from the Manuscript                     #
#  "The complexity of a ‘simple' predator-prey system:                 #
#   Nontrophic positive interactions generate unsuspected dynamics"    #
#======================================================================#

# Load required libraries
library(tidyverse)
library(plyr)

#--------------------------------------------------------------------#
# 1. Load and Prepare Data
#--------------------------------------------------------------------#

# Read the simulation results file containing the last recorded values
last_time <- read.csv('last_time.csv')

# Set factor levels for interaction types to maintain consistency in plots
last_time$interaction_type <- factor(last_time$interaction_type, 
                                     levels = c("Only trophic", 
                                                "Recruitment facilitation",
                                                "Refuge provision", 
                                                "Both NTI"))

# Ensure prey abundance is non-negative
last_time$Prey[last_time$Prey < 0] <- 0

# Compute recruitment levels and observed predator mortality rates
last_time <- last_time %>% mutate(
  gL_P   = g * Lp,
  hL_V   = h * Lv,
  c_obs  = (g + (g_prime - g) * (Prey^a / (Prey^a + Vor^a))) * Lp,  # Observed predator recruitment
  m_obs  = (mmax - (mmax - mi) * (Prey^a / (Prey^a + Vos^a)))  # Observed predator mortality
)

#--------------------------------------------------------------------#
# 2. Subset Data for Plotting
#--------------------------------------------------------------------#

# Select relevant parameter combinations for visualization
absolute_change <- rbind(
  filter(last_time, interaction_type == 'Only trophic' & g == 0.1 & h == 0.5 & mmax == 0.255),
  filter(last_time, interaction_type == 'Recruitment facilitation' & g == 0.1 & g_prime %in% c(0.1, 0.5) & h == 0.5 & mmax == 0.255),
  filter(last_time, interaction_type == 'Refuge provision' & g == 0.1 & h == 0.5 & mi == 0.005 & mmax == 0.255),
  filter(last_time, interaction_type == 'Both NTI' & g == 0.1 & g_prime %in% c(0.1, 0.5) & h == 0.5 & mi == 0.005 & mmax == 0.255)
)

# Create a categorical variable to indicate predator presence or absence
absolute_change <- absolute_change %>% mutate(
  predator_pre_ab = ifelse(Lp == 0, "Absent", "Present")
)

#--------------------------------------------------------------------#
# 3. Generate Figures
#--------------------------------------------------------------------#

#----------------------#
# Figure 2a - Predator Recruitment
#----------------------#

# Function to generate recruitment plots for different predator larval availabilities

plot_recruitment <- function(Lp_value, file_name, y_limit, y_breaks) {
  ggplot(filter(absolute_change, a == 2 & d == 10 & Lp == Lp_value),
         aes(x = hL_V, y = c_obs)) +
    geom_line(aes(linetype = interaction_type, col = factor(Lp)), linewidth = 1.5) +
    scale_y_continuous(limits = y_limit, breaks = y_breaks) +
    scale_colour_manual(values = c("black", "black", "black")) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
    theme_bw() +
    labs(
      y = 'Predator recruitment (c)',
      x = expression('Prey arrival rate (' * italic(hL[V]) * ')'),
      linetype = 'Interaction type'
    ) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 18),
      strip.background = element_blank(),
      legend.title = element_text(size = 14),
      legend.key.size = unit(1, 'cm'),
      legend.position = 'right'
    )
  
  ggsave(file_name, width = 1334, height = 864, units = "px", dpi = 192)
}

# Generate recruitment plots for Lp = 2, 20, 200
plot_recruitment(2, "fig_2a_cobs_LP2.png", c(0, 1), seq(0, 1, by = 0.2))
plot_recruitment(20, "fig_2a_cobs_LP20.png", c(0, 10), seq(0, 10, by = 2))
plot_recruitment(200, "fig_2a_cobs_LP200.png", c(0, 100), seq(0, 100, by = 20))

#-------------------------------------#
# Figure 2b - Predator Mortality
#-------------------------------------#

# Function to generate mortality plots for different predator larval availabilities
plot_mortality <- function(Lp_value, file_name, y_limit, y_breaks) {
  ggplot(filter(absolute_change, a == 2 & d == 10 & Lp == Lp_value),
         aes(x = hL_V, y = m_obs)) +
    geom_line(aes(linetype = interaction_type, col = factor(Lp)), linewidth = 1.5) +
    scale_y_continuous(limits = y_limit, breaks = y_breaks) +
    scale_colour_manual(values = c("black", "black", "black")) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
    theme_bw() +
    labs(
      y = expression('Predator mortality ('*italic(m["P"])*')'),
      x = expression('Prey arrival rate (' * italic(hL[V]) * ')'),
      linetype = 'Interaction type'
    ) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 18),
      strip.background = element_blank(),
      legend.title = element_text(size = 14),
      legend.key.size = unit(1, 'cm'),
      legend.position = 'right'
    )
  
  ggsave(file_name, width = 1334, height = 864, units = "px", dpi = 192)
}

# Generate mortality plots for Lp = 2, 20, 200
plot_mortality(2, "fig_2b_mobs_Lp2.png", c(0.09, 0.26), seq(0.1, 0.25, by = 0.05))
plot_mortality(20, "fig_2b_mobs_Lp20.png", c(0.09, 0.26), seq(0.1, 0.25, by = 0.05))
plot_mortality(200, "fig_2b_mobs_Lp200.png", c(0.09, 0.26), seq(0.1, 0.25, by = 0.05))

#-------------------------------------#
# Figure 2c - Predator Abundance
#-------------------------------------#

# Function to generate predator abundance plots for different predator larval availabilities
plot_predator_abundance <- function(Lp_value, file_name, y_limit, y_breaks) {
  ggplot(filter(absolute_change, a == 2 & d == 10 & Lp == Lp_value),
         aes(x = hL_V, y = Predator)) +
    geom_line(aes(linetype = interaction_type, col = factor(Lp)), linewidth = 1.5) +
    scale_y_continuous(limits = y_limit, breaks = y_breaks) +
    scale_colour_manual(values = c("black", "black", "black")) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
    theme_bw() +
    labs(
      y = expression('Predator abundance'),
      x = expression('Prey arrival rate (' * italic(hL[V]) * ')'),
      linetype = 'Interaction type'
    ) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 18),
      strip.background = element_blank(),
      legend.title = element_text(size = 14),
      legend.key.size = unit(1, 'cm'),
      legend.position = 'right'
    )
  
  ggsave(file_name, width = 1334, height = 864, units = "px", dpi = 192)
}

# Generate predator abundance plots for Lp = 2, 20, 200
plot_predator_abundance(2, "fig_2c_predatorab_lp2.png", c(0, 40), seq(0, 40, by = 10))
plot_predator_abundance(20, "fig_2c_predatorab_lp20.png", c(0, 40), seq(0, 40, by = 10))
plot_predator_abundance(200, "fig_2c_predatorab_lp200.png", c(0, 40), seq(0, 40, by = 10))

#-------------------------------------#
# Figure 2d - Prey Abundance
#-------------------------------------#

# Function to generate prey abundance plots for different predator larval availabilities
plot_prey_abundance <- function(Lp_value, file_name) {
  ggplot(filter(absolute_change, a == 2 & d == 10 & Lp %in% c(0, Lp_value)),
         aes(x = hL_V, y = Prey)) +
    geom_line(aes(linetype = interaction_type, col = predator_pre_ab), linewidth = 1.5) +
    scale_colour_manual(values = c("gray", "black")) +
    theme_bw() +
    geom_hline(aes(yintercept = 3000), color = 'gray', linetype = "dashed") +
    geom_hline(aes(yintercept = 7000), color = 'gray', linetype = "dashed") +
    scale_linetype_manual(values = c("solid", "dashed", "12", "dotdash")) +
    labs(
      y = expression('Prey abundance'),
      x = expression('Prey arrival rate (' * italic(hL[V]) * ')'),
      linetype = "Interaction type"
    ) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 18),
      strip.background = element_blank(),
      legend.title = element_text(size = 14),
      legend.key.size = unit(1, 'cm'),
      legend.position = 'right'
    )
  
  ggsave(file_name, width = 1334, height = 864, units = "px", dpi = 192)
}

# Generate prey abundance plots for Lp = 2, 20, 200
plot_prey_abundance(2, "fig_2d_preyab_lp2.png")
plot_prey_abundance(20, "fig_2d_preyab_lp20.png")
plot_prey_abundance(200, "fig_2d_preyab_lp200.png")


#---------------------------------------------------------------------------------#
# Figure 3a - Predator vs Prey Abundance: ony trophic interaction
#---------------------------------------------------------------------------------#

ggplot(filter(last_time, interaction_type %in% "Only trophic" &
                Lp %in% c(2,10, 20,30) &
                mmax %in% c(0.255)),
       aes(x=Prey, y=Predator, shape = factor(Lp)))+
  geom_line(col='black')+
  geom_hline(aes(yintercept=45.45), linetype='dashed',col="gray")+
  geom_vline(aes(xintercept=10000), linetype='dashed',col="gray")+
  geom_point(aes(fill=hL_V, shape=factor(Lp)), size=5, col='black')+
  scale_x_continuous(breaks = seq(from = 0, to =10000, by=2500),limits = c(0,10000))+
  facet_wrap(interaction_type~.)+
  scale_fill_viridis_c(breaks=seq(0,2500, by=500))+
  scale_shape_manual(values=c(21,22,23,24))+
  labs(x='Prey abundance', y=bquote('Predator abundance'),
       fill=bquote('Prey local arrival rate'~'('~italic(hL[V])*')'),
       shape=bquote('Predator local larval availability'~'('~italic(L[p])*')'))+
  theme_bw()+
  theme(axis.title = element_text(size=18),
        axis.text  = element_text(size=16),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.key.height = unit(1, 'cm'),
        legend.key.size = unit(1, 'cm'),
        panel.grid = element_blank(),
        strip.background = element_rect(fill='white'),
        strip.text = element_text(size=18),
        legend.position = 'none')

ggsave("fig_3a_trophic.png", 
       width = 1250, height = 1250, units = "px", dpi=300)

#---------------------------------------------------------------------------------#
# Figure 3b - Predator vs Prey Abundance: predation and recruitment facilitation
#---------------------------------------------------------------------------------#


ggplot(filter(last_time, interaction_type %in% "Recruitment facilitation" &
                g_prime == 0.5 &
                Lp %in% c(2, 10, 20,30) &
                mmax %in% c(0.255)),
       aes(x=Prey, y=Predator, shape = factor(Lp)))+
  geom_line(col='black')+
  geom_point(aes(fill=hL_V, shape=factor(Lp)), size=5, col='black')+
  scale_x_continuous(breaks = seq(from = 0, to =10000, by=2500),limits = c(0,10000))+
  geom_hline(aes(yintercept=45.45), linetype='dashed',col="gray")+
  geom_vline(aes(xintercept=10000), linetype='dashed',col="gray")+
  facet_wrap(interaction_type~.)+
  scale_fill_viridis_c(breaks=seq(0,2500, by=500))+
  scale_shape_manual(values=c(21,22,23,24))+
  labs(x='Prey abundance', y=bquote('Predator abundance'),
       fill=bquote('Prey local arrival rate'~'('~italic(hL[V])*')'),
       shape=bquote('Predator local larval availability'~'('~italic(L[p])*')'))+
  theme_bw()+
  theme(axis.title = element_text(size=18),
        axis.text  = element_text(size=16),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.key.height = unit(1, 'cm'),
        legend.key.size = unit(1, 'cm'),
        panel.grid = element_blank(),
        strip.background = element_rect(fill='white'),
        strip.text = element_text(size=18),
        legend.position = 'none')

ggsave("fig_3b_recruitment_facilitation.png", 
       width = 1250, height = 1250, units = "px", dpi=300)

#---------------------------------------------------------------------------------#
# Figure 3c - Predator vs Prey Abundance: predation and refuge provision
#---------------------------------------------------------------------------------#

ggplot(filter(last_time, interaction_type %in% "Refuge provision" &
                Lp %in% c(2,10, 20,30) &
                mmax %in% c(0.255),
              mi %in% c(0.005)),
       aes(x=Prey, y=Predator, shape = factor(Lp)))+
  geom_line(col='black')+
  geom_point(aes(fill=hL_V, shape=factor(Lp)), size=5, col='black')+
  scale_x_continuous(breaks = seq(from = 0, to =10000, by=2500),limits = c(0,10000))+
  geom_hline(aes(yintercept=45.45), linetype='dashed',col="gray")+
  geom_vline(aes(xintercept=10000), linetype='dashed',col="gray")+
  facet_wrap(interaction_type~.)+
  scale_fill_viridis_c(breaks=seq(0,2500, by=500))+
  scale_shape_manual(values=c(21,22,23,24))+
  labs(x='Prey abundance', y=bquote('Predator abundance'),
       fill=bquote('Prey local arrival rate'~'('~italic(hL[V])*')'),
       shape=bquote('Predator local larval availability'~'('~italic(L[p])*')'))+
  theme_bw()+
  theme(axis.title = element_text(size=18),
        axis.text  = element_text(size=16),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.key.height = unit(1, 'cm'),
        legend.key.size = unit(1, 'cm'),
        panel.grid = element_blank(),
        strip.background = element_rect(fill='white'),
        strip.text = element_text(size=18),
        legend.position = 'none')

ggsave("fig_3c_refuge_provision.png", 
       width = 1250, height = 1250, units = "px", dpi=300)

#---------------------------------------------------------------------------------#
# Figure 3d - Predator vs Prey Abundance: all interactions simultaneously
#---------------------------------------------------------------------------------#

ggplot(filter(last_time, interaction_type %in% "Both NTI" &
                g_prime == 0.5 &
                Lp %in% c(2,10, 20,30) &
                mi %in% c(0.005),
              mmax %in% c(0.255)),
       aes(x=Prey, y=Predator, shape = factor(Lp)))+
  geom_line(col='black')+
  geom_point(aes(fill=hL_V, shape=factor(Lp)), size=5, col='black')+
  scale_x_continuous(breaks = seq(from = 0, to =10000, by=2500),limits = c(-1,10000))+
  geom_hline(aes(yintercept=45.45), linetype='dashed',col="gray")+
  geom_vline(aes(xintercept=10000), linetype='dashed',col="gray")+
  facet_wrap(interaction_type~.)+
  scale_fill_viridis_c(breaks=seq(0,2500, by=500))+
  scale_shape_manual(values=c(21,22,23,24))+
  labs(x='Prey abundance', y=bquote('Predator abundance'),
       fill=bquote('Prey local arrival rate'~'('~italic(hL[V])*')'),
       shape=bquote('Predator local larval availability'~'('~italic(L[p])*')'))+
  theme_bw()+
  theme(axis.title = element_text(size=18),
        axis.text  = element_text(size=16),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.key.height = unit(1, 'cm'),
        legend.key.size = unit(1, 'cm'),
        panel.grid = element_blank(),
        strip.background = element_rect(fill='white'),
        strip.text = element_text(size=18),
        legend.position = 'none')

ggsave("fig_3d_bothNTI.png", 
       width = 1250, height = 1250, units = "px", dpi=300)

#--------------------------------------------------------------------#
# 4. Compute Predator-Prey Coupling Strength (Derivative Analysis)
#--------------------------------------------------------------------#

# Define conditions for filtering based on interaction type

interaction_filters <- list(
  "Only trophic" = quote(mmax == 0.255),
  "Recruitment facilitation" = quote(mmax == 0.255 & g_prime == 0.5),
  "Refuge provision" = quote(mmax == 0.255 & mi == 0.005),
  "Both NTI" = quote(mmax == 0.255 & mi == 0.005 & g_prime == 0.5)
)

# Function to filter data based on conditions
filter_interaction_data <- function(interaction_type, condition) {
  last_time %>%
    filter(interaction_type == !!interaction_type) %>%
    filter(!!condition)
}

# Apply filter function dynamically
regresion_data <- map_dfr(names(interaction_filters), 
                          ~filter_interaction_data(.x, interaction_filters[[.x]]))

# Keep only cases where Lp <= 30
regresion_data <- regresion_data %>% filter(Lp <= 30)

# Order data for derivative calculation
regresion_data <- regresion_data %>%
  arrange(interaction_type, gL_P, Prey)

#----------------------------------------------#
# Function to Compute Predator-Prey Coupling Strength
#----------------------------------------------#

compute_derivatives <- function(df, variable) {
  # Calculate numerical derivative (ΔPredator / Δvariable)
  derivatives <- diff(df$Predator) / diff(df[[variable]])
  hL_V <- head(df$hL_V, -1)
  data.frame(hL_V = hL_V, dpdv = derivatives)
}

# Apply derivative function for different variables
derivatives_prey <- plyr::ddply(regresion_data, ~ interaction_type + gL_P, 
                                compute_derivatives, variable = "Prey")

derivatives_Lv <- plyr::ddply(regresion_data, ~ interaction_type + Lp, 
                              compute_derivatives, variable = "Lv")

#----------------------------------------------#
# Plot Predator-Prey Coupling Strength
#----------------------------------------------#

ggplot(filter(derivatives_prey, gL_P != 0 & hL_V != 0), 
       aes(x = hL_V, y = gL_P, z = dpdv)) +
  geom_contour_filled(bins = 20) +
  geom_contour(breaks = 0.000495, linewidth = 1, col = 'white', linetype = 'dashed') +
  facet_wrap(. ~ interaction_type) +
  labs(
    x = bquote('Prey recruitment rate (' * italic(s) * ')'),
    y = bquote('Predator recruitment rate (' * italic(c) * ')'),
    fill = 'Derivative'
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 14),
    legend.title.align = 0.5,
    strip.background = element_blank(),
    strip.text = element_text(size = 16),
    legend.position = 'right'
  )

ggsave("fig_4.png", width = 3276, height = 2303, units = "px", dpi = 330)