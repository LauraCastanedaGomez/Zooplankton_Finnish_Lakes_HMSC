
###############################################################
### LIBRARIES
###############################################################
library(readxl)
library(cluster)   
library(dplyr)
library(ggplot2)
library(tidyr)


###############################################################
### LOAD TRAITS
###############################################################

traits_df <- read_excel("//ad.helsinki.fi/home/l/lauracas/Desktop/NoCapture_BUENO/data/Todos_traits.xlsx")

# First column = species names
species <- traits_df[[1]]
traits_df <- traits_df[,-1]
rownames(traits_df) <- species

# Convert binary/categorical traits to factors
traits_df[] <- lapply(traits_df, function(x){
  if(is.numeric(x) && all(x %in% c(0,1))) as.factor(x) else x
})

# Keep body size as continuous
traits_df$Body_size <- as.numeric(traits_df$Body_size)

rownames(traits_df) <- species


###############################################################
### LOAD RESIDUAL CO-OCCURRENCE MATRIX (Ω)
###############################################################
omega_df <- read_excel("//ad.helsinki.fi/home/l/lauracas/Desktop/NoCapture_BUENO/data/matrix_coocurrence.xlsx")
omega_mat <- as.matrix(omega_df[,-1])
rownames(omega_mat) <- colnames(omega_mat) <- omega_df[[1]]


###############################################################
### FUNCTION: PERMUTATION TEST
###############################################################

permutation_test <- function(dist_trait, omega_mat, nperm = 9999){
  
  obs_trait  <- dist_trait[upper.tri(dist_trait)]
  obs_omega  <- omega_mat[upper.tri(omega_mat)]
  
  # observed correlation
  obs_rho <- cor(obs_trait, obs_omega, method = "spearman")
  
  null_rho <- numeric(nperm)
  species <- rownames(dist_trait)
  
  for(i in 1:nperm){
    perm <- sample(species)
    dist_perm <- dist_trait[perm, perm]
    perm_trait_vals <- dist_perm[upper.tri(dist_perm)]
    
    null_rho[i] <- cor(perm_trait_vals, obs_omega, method = "spearman")
  }
  
  p_val <- (sum(abs(null_rho) >= abs(obs_rho)) + 1) / (nperm + 1)
  
  list(
    observed_rho = obs_rho,
    p_value = p_val,
    null_distribution = null_rho
  )
}


###############################################################
### RUN THE TEST FOR ALL TRAITS
###############################################################

run_permutation_for_all_traits <- function(traits_df, omega_mat, nperm = 9999){
  
  results <- data.frame(
    Trait = character(),
    Observed_rho = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE)
  
  for(trait in colnames(traits_df)){
    
    trait_vec <- data.frame(T = traits_df[[trait]])
    
    dist_trait <- as.matrix(daisy(trait_vec, metric="gower"))
    rownames(dist_trait) <- rownames(traits_df)
    colnames(dist_trait) <- rownames(traits_df)
    
    res <- permutation_test(dist_trait, omega_mat, nperm)
    
    results <- rbind(
      results,
      data.frame(
        Trait = trait,
        Observed_rho = res$observed_rho,
        P_value = res$p_value
      )
    )
    
    assign(paste0("null_", trait), res$null_distribution, envir = .GlobalEnv)
  }
  
  results
}

results <- run_permutation_for_all_traits(traits_df, omega_mat, nperm = 9999)


###############################################################
### FIGURE 1 — FOREST PLOT
###############################################################

results_plot <- results %>%
  mutate(Trait = factor(Trait, levels = Trait[order(Observed_rho)]))

ggplot(results_plot, aes(x = Observed_rho, y = Trait)) +
  geom_point(aes(color = P_value < 0.05), size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values=c("black", "red")) +
  labs(title = "Trait effects on residual co-occurrence (Ω)",
       x = "Spearman rho",
       y = "Trait",
       color = "Significant") +
  theme_minimal(base_size = 14)


###############################################################
### FIGURE 2 — HEATMAP
###############################################################

results_heat <- results_plot %>%
  mutate(Sign = ifelse(P_value < 0.05, "*", ""))

ggplot(results_heat, aes(x = Trait, y = 1, fill = Observed_rho)) +
  geom_tile() +
  geom_text(aes(label = Sign), size = 6) +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  theme_minimal(base_size = 14) +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = element_text(angle=45, hjust=1)) +
  labs(title = "Correlation between traits and Ω",
       fill = "rho")

