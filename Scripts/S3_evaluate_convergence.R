
##################################################################################################
# MCMC CONVERGENCE ANALYSIS
##################################################################################################

set.seed(1)

showBeta = TRUE
showGamma = TRUE
showOmega = TRUE
maxOmega = 100
showRho = TRUE
showAlpha = TRUE

localDir = "."
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)

library(Hmsc)
library(colorspace)
library(vioplot)

samples_list <- c(5, 250, 250, 250, 500)
thin_list <- c(1, 1, 10, 100, 500)
nst <- length(thin_list)
nChains <- 4

text.file <- file.path(resultDir, "MCMC_convergence.txt")
cat("MCMC Convergence statistics\n\n", file = text.file, sep = "")

ma.beta <- list()
ma.gamma <- list()
ma.omega <- list()

na.beta <- c()
na.gamma <- c()
na.omega <- c()

Lst <- 1
while (Lst <= nst) {
  thin <- thin_list[Lst]
  samples <- samples_list[Lst]
  
  filename <- file.path(modelDir, paste0("models_thin_", thin,
                                         "_samples_", samples,
                                         "_chains_", nChains, ".Rdata"))
  if (file.exists(filename)) {
    load(filename)
    cat("\n", filename, "\n\n", file = text.file, append = TRUE)
    nm <- length(models)
    for (j in 1:nm) {
      mpost <- convertToCodaObject(models[[j]], spNamesNumbers = c(TRUE, FALSE),
                                   covNamesNumbers = c(TRUE, FALSE))
      nr <- models[[j]]$nr
      cat("\n", names(models)[j], "\n\n", file = text.file, append = TRUE)
      
      if (showBeta && !is.null(mpost$Beta)) {
        psrf <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
        if (!is.null(psrf)) {
          ma.beta[[length(ma.beta) + 1]] <- psrf[, 1]
          na.beta <- c(na.beta, paste0(thin, ",", samples))
        }
      }
      
      if (showGamma && !is.null(mpost$Gamma)) {
        psrf <- gelman.diag(mpost$Gamma, multivariate = FALSE)$psrf
        if (!is.null(psrf)) {
          ma.gamma[[length(ma.gamma) + 1]] <- psrf[, 1]
          na.gamma <- c(na.gamma, paste0(thin, ",", samples))
        }
      }
      
      if (showOmega && nr > 0) {
        for (k in 1:nr) {
          tmp <- mpost$Omega[[k]]
          if (!is.null(tmp)) {
            z <- dim(tmp[[1]])[2]
            if (z > maxOmega) {
              sel <- sample(1:z, size = maxOmega)
              for (i in 1:length(tmp)) tmp[[i]] <- tmp[[i]][, sel]
            }
            psrf <- gelman.diag(tmp, multivariate = FALSE)$psrf
            if (!is.null(psrf)) {
              ma.omega[[length(ma.omega) + 1]] <- psrf[, 1]
              na.omega <- c(na.omega, paste0(thin, ",", samples))
            }
          }
        }
      }
    }
  }
  Lst <- Lst + 1
}

# UNIFORM COLOR
col_single <- "lightblue"

pdf(file = file.path(resultDir, "MCMC_convergence.pdf"))

if (showBeta && length(ma.beta) > 0) {
  par(mfrow = c(2, 1))
  vioplot(ma.beta, col = col_single, names = na.beta,
          ylim = c(0, max(unlist(ma.beta), na.rm = TRUE)),
          main = "psrf(beta)")
  vioplot(ma.beta, col = col_single, names = na.beta,
          ylim = c(0.90, 1.10),
          main = "psrf(beta)")
}

if (showGamma && length(ma.gamma) > 0) {
  par(mfrow = c(2, 1))
  vioplot(ma.gamma, col = col_single, names = na.gamma,
          ylim = c(0, max(unlist(ma.gamma), na.rm = TRUE)),
          main = "psrf(gamma)")
  vioplot(ma.gamma, col = col_single, names = na.gamma,
          ylim = c(0.90, 1.10),
          main = "psrf(gamma)")
}

if (showOmega && length(ma.omega) > 0) {
  par(mfrow = c(2, 1))
  vioplot(ma.omega, col = col_single, names = na.omega,
          ylim = c(0, max(unlist(ma.omega), na.rm = TRUE)),
          main = "psrf(omega)")
  vioplot(ma.omega, col = col_single, names = na.omega,
          ylim = c(0.90, 1.10),
          main = "psrf(omega)")
}

dev.off()

library(coda)

# Example with Beta
ess_beta <- effectiveSize(mpost$Beta)

# Example with Gamma
ess_gamma <- effectiveSize(mpost$Gamma)

# Example with Omega (per block)
ess_omega <- lapply(mpost$Omega, effectiveSize)


summary(ess_beta)   # Min, 1Q, Median, Mean, 3Q, Max

summary(ess_gamma)

summary(unlist(ess_omega))  # combined summary



library(coda)

# Assuming mpost is the object from convertToCodaObject() already loaded
# ESS summaries
ess_beta  <- effectiveSize(mpost$Beta)
ess_gamma <- effectiveSize(mpost$Gamma)

# Omega: it is a list (possibly by levels), combine
ess_omega_list <- lapply(mpost$Omega, effectiveSize)
ess_omega <- unlist(ess_omega_list)

# Basic summaries
cat("ESS Beta:\n"); print(summary(ess_beta))
cat("\nESS Gamma:\n"); print(summary(ess_gamma))
cat("\nESS Omega (combined):\n"); print(summary(ess_omega))

# Proportion relative to total saved samples (example)
# total_samples_saved = samples * nChains  # adjust 'samples' and 'nChains' as needed
total_saved <- length(as.mcmc(mpost$Beta[[1]])) * length(mpost$Beta) # example for list
# better to calculate total saved per parameter:
niter_per_chain <- nrow(mpost$Beta[[1]])    # iterations saved per chain for Beta
nchains <- length(mpost$Beta)               # number of chains (if Beta is list per chain)
total_saved_per_param <- niter_per_chain * nchains
cat("\nExample total saved per parameter:", total_saved_per_param, "\n")

# Identify the worst (minimum) ESS and check their traces
worst_ess_omega <- sort(ess_omega)[1:10]   # 10 parameters with lowest ESS
print(worst_ess_omega)

# Optional: plot traceplots of worst parameters
# Need to locate the index/position in the original mpost$Omega; example:
bad_names <- names(worst_ess_omega)
par(mfrow = c(3,1))
for (nm in bad_names[1:min(6, length(bad_names))]) {
  # locate in Omega object (adjust according to structure)
  # suppose unif_omega_mat is the matrix with columns named by parameters
  # plot(as.mcmc(unif_omega_mat[, nm]), main = paste("Trace:", nm))
  NULL
}

###################################################

# Filter only the final run
final_index_beta <- which(na.beta == "500,500")
final_index_gamma <- which(na.gamma == "500,500")
final_index_omega <- which(na.omega == "500,500")

ma.beta.final <- ma.beta[final_index_beta]
ma.gamma.final <- ma.gamma[final_index_gamma]
ma.omega.final <- ma.omega[final_index_omega]

library(vioplot)

# Open pdf to save
pdf("MCMC_convergence_final.pdf", width=10, height=8)

par(mfrow=c(3,2))  # 3 rows (Beta, Gamma, Omega) x 2 columns (zoom out / zoom in)

# --- BETA ---
vioplot(ma.beta.final, col="lightblue", names=final_index_beta,
        ylim=c(0, max(unlist(ma.beta.final), na.rm=TRUE)),
        main="PSRF Beta (zoomed out)")
vioplot(ma.beta.final, col="lightblue", names=final_index_beta,
        ylim=c(0.9, 1.1),
        main="PSRF Beta (zoomed in)")

# --- GAMMA ---
vioplot(ma.gamma.final, col="lightgreen", names=final_index_gamma,
        ylim=c(0, max(unlist(ma.gamma.final), na.rm=TRUE)),
        main="PSRF Gamma (zoomed out)")
vioplot(ma.gamma.final, col="lightgreen", names=final_index_gamma,
        ylim=c(0.9, 1.1),
        main="PSRF Gamma (zoomed in)")

# --- OMEGA ---
vioplot(ma.omega.final, col="lightpink", names=final_index_omega,
        ylim=c(0, max(unlist(ma.omega.final), na.rm=TRUE)),
        main="PSRF Omega (zoomed out)")
vioplot(ma.omega.final, col="lightpink", names=final_index_omega,
        ylim=c(0.9, 1.1),
        main="PSRF Omega (zoomed in)")

dev.off()

