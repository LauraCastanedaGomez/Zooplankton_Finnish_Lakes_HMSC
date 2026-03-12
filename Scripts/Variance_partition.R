
####VP######

####VP TRAITS####

###Plot by variables###

library(scales)
library(tidyverse)
library(ggplot2)

# Extract the values for traits
vp_traits <- VP[["R2T"]][["Beta"]]

vp_traits

# Convert to data.frame and remove the intercept
vp_traits_df <- tibble(
  Trait = names(vp_traits),
  Variance = as.numeric(vp_traits)
) %>%
  filter(Trait != "(Intercept)") %>%   # <-- remove intercept
  arrange(Variance)

vp_traits_df

vp_traits_df <- vp_traits_df %>%
  mutate(Trait = trimws(Trait))

unique(vp_traits_df$Trait)

# Assign groups
vp_traits_df <- vp_traits_df %>%
  mutate(
    VarGroup = case_when(
      Trait %in% c("conductivity","pH","water_temp","water_depth",
                   "colour","chla","lake_area") ~ "Physicochemical",
      Trait %in% c("Tavg","precip") ~ "Climate",
      Trait %in% c("human","forest") ~ "Land-use"
    )
  )

# Plot
ggplot(vp_traits_df, aes(x = reorder(Trait, Variance), y = Variance, fill=VarGroup)) +
  geom_col(color = "white", width = 0.7) + # bars per trait
  coord_flip()  +
  scale_fill_manual(values = c(
    "Physicochemical" = "#5E81AC",
    "Climate" = "#A3BE8C",
    "Land-use" = "#BF616A"
  )) +
  scale_x_discrete(labels = c(
    "conductivity" = "Conductivity",
    "pH" = "pH",
    "water_temp" = "Water temp",
    "water_depth" = "Water depth",
    "colour" = "Colour",
    "chla" = "Chla",
    "lake_area" = "Lake area",
    "Tavg" = "Temp",
    "precip" = "Precip",
    "human" = "Human impact",
    "forest" = "Forest"
  )) +
  geom_text(aes(label = paste0(round(Variance*100, 1), "%")), 
            hjust = -0.1, size = 3.5) +  # labels outside bars
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Explained variance by functional traits",
    x = "Environmental variable",
    y = "Explained variance (R²)",
    fill = "Variable type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 11, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
  )


###Plot by variable group###

library(dplyr)
library(tibble)
library(tidyr)

# 1️⃣ Extract trait values
vp_traits <- VP[["R2T"]][["Beta"]]

# 2️⃣ Create data.frame
vp_traits_df <- tibble(
  Trait = names(vp_traits),
  Variance = as.numeric(vp_traits)
) %>%
  filter(Trait != "(Intercept)")

# 3️⃣ Create variable groups
var_groups <- tibble(
  Trait = c("conductivity","pH","water_temp","water_depth",
            "colour","chla","lake_area","Tavg","precip",
            "human","forest"),
  VarGroup = c("Physicochemical","Physicochemical","Physicochemical","Physicochemical",
               "Physicochemical","Physicochemical","Physicochemical",
               "Climate","Climate",
               "Land-use","Land-use")
)

# 4️⃣ Join groups and calculate mean variance per group
vp_grouped <- vp_traits_df %>%
  left_join(var_groups, by = "Trait") %>%
  group_by(VarGroup) %>%
  summarise(Variance_mean = mean(Variance, na.rm = TRUE), .groups = "drop") %>%
  arrange(Variance_mean)

# 5️⃣ Plot vertical bars
ggplot(vp_grouped, aes(x = reorder(VarGroup, Variance_mean), y = Variance_mean)) +
  geom_col(fill = "#5E81AC", color = "white", width = 0.7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Explained variance by functional traits",
    x = "Variable type",
    y = "Mean explained variance (R²)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),  # rotate labels
    axis.text.y = element_text(size = 11, color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
  )

vp_grouped <- vp_grouped %>%
  mutate(Variance_mean_pct = Variance_mean * 100)

vp_grouped

vp_grouped <- vp_grouped %>%
  mutate(VarGroup = factor(VarGroup, levels = c("Climate", "Land-use", "Physicochemical")))

ggplot(vp_grouped, aes(x = VarGroup, y = Variance_mean, fill = VarGroup)) +
  geom_col(color = "white", width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(Variance_mean*100,1),"%")), 
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c(
    "Climate" = "#A3BE8C",
    "Land-use" = "#BF616A",
    "Physicochemical" = "#5E81AC"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0,0.15))) +
  labs(
    title = "Explained variance by functional traits",
    x = "Variable type",
    y = "Mean explained variance (R²)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "black"),
    axis.title = element_text(size = 11, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

####VP TAXONOMIC####

library(scales)
library(tidyverse)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)

# Compute variance partitioning
VP = computeVariancePartitioning(m)

# Transform to long format
vp_long <- VP$vals %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%  # each row is an environmental variable
  pivot_longer(
    cols = -Variable,
    names_to = "Species",
    values_to = "Variance"
  )

# Clean species names
species_groups <- read_excel("//ad.helsinki.fi/home/l/lauracas/Desktop/NoCapture_BUENO/data/Groups.xlsx") %>%
  mutate(Sp = trimws(Sp))

vp_long <- vp_long %>%
  left_join(species_groups, by = c("Species" = "Sp")) %>%
  filter(!is.na(Group))  # remove species without group

# Compute mean variance by variable and taxon group
vp_grouped <- vp_long %>%
  group_by(Variable, Group) %>%
  summarise(Variance_mean = mean(Variance, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = scales::percent(Variance_mean, accuracy = 1))

# Optional: rename groups
vp_grouped <- vp_grouped %>%
  mutate(Group = recode(Group,
                        "Rotifer" = "Rotifera",
                        "Cladocera" = "Cladocera"))

# Horizontal plot with side-by-side bars
ggplot(vp_grouped, aes(x = reorder(Variable, Variance_mean), y = Variance_mean, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), color = "white", width = 0.7) +
  geom_text(aes(label = label),
            position = position_dodge(width = 0.8), # labels above each bar
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_x_discrete(labels = c(
    "conductivity" = "Conductivity",
    "pH" = "pH",
    "water_temp" = "Water temp",
    "water_depth" = "Water depth",
    "colour" = "Colour",
    "chla" = "Chla",
    "lake_area" = "Lake area",
    "Tavg" = "Temp",
    "precip" = "Precip",
    "human" = "Human impact",
    "forest" = "Forest",
    "Random: Index_random_geo"="Random"
  )) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c(
    "Cladocera" = "#D08770",
    "Rotifera" = "#88C0D0"
  )) +
  labs(
    title = "Explained variance by environmental variables and taxon group",
    x = "Environmental variable",
    y = "Explained variance (R²)",
    fill = "Taxon group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

