# Red panda spatial occupancy analysis for supplementary material.
# Run from this folder. Input contains red panda camera detection histories
# (including non-detections) and only the covariates needed by these models.

library(spOccupancy)

input <- readRDS("inputs/red_panda_spatial_model_ready_data.rds")
model_data <- input$data
y <- model_data$y
occ <- model_data$occ.covs
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE)

# All models use the same exponential NNGP process and priors.
phi_primary <- c(3 / 100, 3 / 5)  # effective spatial range: 5–100 km
sigma_sq_prior <- c(0, 5)
mcmc <- list(n.batch = 500, batch.length = 200, n.burn = 20000,
             n.thin = 10, n.chains = 3, n.neighbors = 15,
             n.omp.threads = 8)

fit_model <- function(occ_formula, det_formula, phi_bounds, seed) {
  set.seed(seed)
  spPGOcc(
    occ.formula = occ_formula, det.formula = det_formula,
    data = model_data,
    inits = list(phi = mean(phi_bounds), sigma.sq = 1),
    priors = list(
      beta.normal = list(mean = 0, var = 2.72),
      alpha.normal = list(mean = 0, var = 2.72),
      phi.unif = phi_bounds,
      sigma.sq.unif = sigma_sq_prior
    ),
    cov.model = "exponential", NNGP = TRUE,
    n.neighbors = mcmc$n.neighbors, search.type = "cb",
    n.batch = mcmc$n.batch, batch.length = mcmc$batch.length,
    n.burn = mcmc$n.burn, n.thin = mcmc$n.thin,
    n.chains = mcmc$n.chains, n.omp.threads = mcmc$n.omp.threads,
    verbose = TRUE, n.report = 25
  )
}

# The Rmd selected the lowest-WAIC model only among models meeting its
# Rhat <= 1.01 and ESS >= 400 screen. Retain this selection rule.
selection_pass <- function(fit) {
  parts <- c("beta", "alpha", "theta")
  rh <- unlist(fit$rhat[parts], use.names = FALSE)
  es <- unlist(fit$ESS[parts], use.names = FALSE)
  length(rh) > 0 && length(es) > 0 &&
    all(is.finite(rh)) && all(rh <= 1.01) &&
    all(is.finite(es)) && all(es >= 400)
}

run_stage <- function(formulas, stage, fixed_formula, seed_start) {
  dir.create(file.path(out_dir, stage), showWarnings = FALSE)
  rows <- vector("list", length(formulas))
  files <- setNames(character(length(formulas)), names(formulas))
  for (i in seq_along(formulas)) {
    name <- names(formulas)[i]
    occ_formula <- if (stage == "detection") fixed_formula else formulas[[i]]
    det_formula <- if (stage == "detection") formulas[[i]] else fixed_formula
    fit <- fit_model(occ_formula, det_formula, phi_primary, seed_start + i)
    file <- file.path(out_dir, stage, paste0(name, ".rds"))
    saveRDS(fit, file)
    files[name] <- file
    w <- waicOcc(fit)
    rows[[i]] <- data.frame(
      model = name, WAIC = as.numeric(w["WAIC"]),
      elpd = as.numeric(w["elpd"]), pD = as.numeric(w["pD"]),
      selection_pass = selection_pass(fit)
    )
    rm(fit)
    gc()
  }
  table <- do.call(rbind, rows)
  table <- table[order(table$WAIC), ]
  table$delta_WAIC <- table$WAIC - min(table$WAIC)
  table$weight <- exp(-table$delta_WAIC / 2)
  table$weight <- table$weight / sum(table$weight)
  write.csv(table, file.path(out_dir, paste0(stage, "_model_comparison.csv")),
            row.names = FALSE)
  eligible <- table[table$selection_pass, ]
  if (nrow(eligible) == 0) stop("No ", stage, " model passed the MCMC screen.")
  selected <- eligible$model[which.min(eligible$WAIC)]
  list(selected = selected, file = files[[selected]], table = table, files = files)
}

# Stage 1: vary detection while holding the full quadratic occupancy formula.
full_occ <- ~ elevation + elevation_sq + slope + northness + eastness +
  forest_cover + road_density + stream_density + FLII +
  understory_PAVD + canopy_height

detection_formulas <- list(
  null = ~ 1,
  season = ~ season,
  season_northness = ~ season + det_northness,
  terrain = ~ det_slope + det_northness + det_eastness,
  full = ~ season + det_slope + det_northness + det_eastness
)

detection <- run_stage(detection_formulas, "detection", full_occ, 5000)
best_det <- detection_formulas[[detection$selected]]

# Stage 2: vary occupancy while holding the selected detection formula.
occupancy_formulas <- list(
  full_quadratic = full_occ,
  full_linear = ~ elevation + slope + northness + eastness + forest_cover +
    road_density + stream_density + FLII + understory_PAVD + canopy_height,
  no_canopy_height = ~ elevation + elevation_sq + slope + northness + eastness +
    forest_cover + road_density + stream_density + FLII + understory_PAVD,
  no_understory = ~ elevation + elevation_sq + slope + northness + eastness +
    forest_cover + road_density + stream_density + FLII + canopy_height,
  no_forest_cover = ~ elevation + elevation_sq + slope + northness + eastness +
    road_density + stream_density + FLII + understory_PAVD + canopy_height,
  no_FLII = ~ elevation + elevation_sq + slope + northness + eastness +
    forest_cover + road_density + stream_density + understory_PAVD +
    canopy_height,
  no_stream_density = ~ elevation + elevation_sq + slope + northness + eastness +
    forest_cover + road_density + FLII + understory_PAVD + canopy_height,
  no_road = ~ elevation + elevation_sq + slope + northness + eastness +
    forest_cover + stream_density + FLII + understory_PAVD + canopy_height,
  no_slope = ~ elevation + elevation_sq + northness + eastness + forest_cover +
    road_density + stream_density + FLII + understory_PAVD + canopy_height,
  no_aspect = ~ elevation + elevation_sq + slope + forest_cover + road_density +
    stream_density + FLII + understory_PAVD + canopy_height,
  no_elevation = ~ slope + northness + eastness + forest_cover + road_density +
    stream_density + FLII + understory_PAVD + canopy_height,
  occupancy_core = ~ elevation + elevation_sq + forest_cover + FLII +
    stream_density + understory_PAVD,
  topography_only = ~ elevation + elevation_sq + slope + northness + eastness,
  vegetation_only = ~ forest_cover + understory_PAVD + canopy_height,
  landscape_context = ~ FLII + stream_density + road_density,
  null = ~ 1
)

occupancy <- run_stage(occupancy_formulas, "occupancy", best_det, 8000)
best_occ <- occupancy_formulas[[occupancy$selected]]
final_fit <- readRDS(occupancy$file)
saveRDS(final_fit, file.path(out_dir, "final_spatial_model.rds"))

# Prior sensitivity: refit the selected structures under neighboring
# spatial range priors. These fits do not change the selected formulas.
range_priors <- list(
  local_3_30km = c(3 / 30, 3 / 3),
  primary_5_100km = phi_primary,
  broad_10_200km = c(3 / 200, 3 / 10)
)
prior_rows <- list()
dir.create(file.path(out_dir, "prior_sensitivity"), showWarnings = FALSE)
for (i in seq_along(range_priors)) {
  name <- names(range_priors)[i]
  fit <- if (name == "primary_5_100km") final_fit else
    fit_model(best_occ, best_det, range_priors[[i]], 12000 + i)
  saveRDS(fit, file.path(out_dir, "prior_sensitivity", paste0(name, ".rds")))
  w <- waicOcc(fit)
  prior_rows[[i]] <- data.frame(prior = name, WAIC = as.numeric(w["WAIC"]),
                                elpd = as.numeric(w["elpd"]),
                                pD = as.numeric(w["pD"]))
  if (name != "primary_5_100km") rm(fit)
  gc()
}
write.csv(do.call(rbind, prior_rows),
          file.path(out_dir, "prior_sensitivity_WAIC.csv"), row.names = FALSE)

# Final spatial-model parameters and camera occupancy probabilities.
posterior_summary <- function(samples, component) {
  s <- as.matrix(samples)
  data.frame(
    component = component, parameter = colnames(s),
    mean = colMeans(s), median = apply(s, 2, median),
    lower_95 = apply(s, 2, quantile, 0.025),
    upper_95 = apply(s, 2, quantile, 0.975),
    row.names = NULL
  )
}
parameters <- rbind(
  posterior_summary(final_fit$beta.samples, "occupancy"),
  posterior_summary(final_fit$alpha.samples, "detection"),
  posterior_summary(final_fit$theta.samples, "spatial covariance")
)
write.csv(parameters, file.path(out_dir, "final_parameter_summary.csv"),
          row.names = FALSE)

psi <- as.matrix(final_fit$psi.samples)
site_occupancy <- data.frame(
  camstn = rownames(y), mean = colMeans(psi),
  lower_95 = apply(psi, 2, quantile, 0.025),
  upper_95 = apply(psi, 2, quantile, 0.975)
)
write.csv(site_occupancy, file.path(out_dir, "site_occupancy.csv"),
          row.names = FALSE)
overall_psi <- rowMeans(psi)
write.csv(data.frame(mean = mean(overall_psi),
                     lower_95 = quantile(overall_psi, 0.025),
                     upper_95 = quantile(overall_psi, 0.975)),
          file.path(out_dir, "overall_occupancy.csv"), row.names = FALSE)

# Detection probability is averaged only over sampled camera occasions.
p <- fitted(final_fit)$p.samples
sampled <- !is.na(y)
mean_p_by_draw <- apply(p, 1, function(draw) mean(draw[sampled]))
write.csv(data.frame(mean = mean(mean_p_by_draw),
                     lower_95 = quantile(mean_p_by_draw, 0.025),
                     upper_95 = quantile(mean_p_by_draw, 0.975)),
          file.path(out_dir, "overall_detection.csv"), row.names = FALSE)

# Freeman–Tukey posterior predictive checks, grouped by camera site and
# by seven-day sampling occasion, as in the source Rmd.
ppc_site <- ppcOcc(object = final_fit, fit.stat = "freeman-tukey", group = 1)
ppc_occasion <- ppcOcc(object = final_fit, fit.stat = "freeman-tukey", group = 2)
saveRDS(list(by_site = ppc_site, by_occasion = ppc_occasion),
        file.path(out_dir, "final_spatial_PPC.rds"))
writeLines(c("By camera site:", capture.output(summary(ppc_site)),
             "", "By occasion:", capture.output(summary(ppc_occasion))),
           file.path(out_dir, "final_spatial_PPC_summary.txt"))

# The full-quadratic candidate isolates the shape of the elevation response.
quad <- readRDS(occupancy$files[["full_quadratic"]])
b <- as.matrix(quad$beta.samples)
peak_z <- -b[, "elevation"] / (2 * b[, "elevation_sq"])
hump <- b[, "elevation_sq"] < 0 & is.finite(peak_z)
elevation_m <- read.csv("inputs/site_elevation_m.csv")
peak_m <- mean(elevation_m$elevation_m) +
  sd(elevation_m$elevation_m) * peak_z[hump]
write.csv(data.frame(
  probability_negative_quadratic = mean(b[, "elevation_sq"] < 0),
  median_peak_z_given_hump = if (any(hump)) median(peak_z[hump]) else NA_real_,
  median_peak_m_given_hump = if (any(hump)) median(peak_m) else NA_real_
), file.path(out_dir, "quadratic_elevation_summary.csv"), row.names = FALSE)

saveRDS(list(input = "inputs/red_panda_spatial_model_ready_data.rds",
             selected_detection = detection$selected,
             selected_occupancy = occupancy$selected,
             mcmc = mcmc, primary_phi_prior = phi_primary),
        file.path(out_dir, "analysis_manifest.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
