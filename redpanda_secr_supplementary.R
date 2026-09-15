# =============================================================================
# Red panda density in Bhutan: spatially explicit capture-recapture
# Supplementary analysis script
#
# Reproduces the saved habitat-mask analysis: density, abundance, model
# selection, sex mixture, and threshold sensitivity.
#
# Requires:  R >= 4.3, secr >= 5.4, sf, terra
# Inputs  :  inputs/capthist.rds            capture history with sex covariates
#            inputs/secr_captures.txt       session, ID, occasion, detector
#            inputs/secr_traps.txt          detector, x, y, effort
#            inputs/individual_covariates.csv
#            Layers/Bhutan_DEM.tif          DEM clipped to the national boundary
#            Layers/forest_only.shp         forest cover
# Outputs :  outputs/ - fitted models (.rds) and result tables (.csv)
# Runtime :  the six-model set takes several hours; the threshold sensitivity
#            adds two further fits. Completed fits are cached in outputs/ and
#            reused on a re-run.
#
# Detector coordinates are UTM zone 46N (EPSG:32646). The forest layer is in
# DRUKREF 03 / Bhutan National Grid and is reprojected on the fly.
# =============================================================================

library(secr)
library(sf)
library(terra)

## ---- settings ---------------------------------------------------------------
in_dir      <- "inputs"
out_dir     <- "outputs"
layers_dir  <- "Layers"
trap_epsg   <- 32646
far_thr     <- 10000     # within-identity spatial spread above which the rule applies (m)
buffer_m    <- 20000     # integration buffer around detectors
spacing_m   <- 2000      # mask spacing
elev_range  <- c(2400, 3700)   # suitable elevation band, m (Dorji et al. 2011)
detfn       <- "HN"      # halfnormal; parameter g0
ncores      <- 4
reuse       <- FALSE     # package has no supplied fits; set TRUE for repeat runs
thresholds  <- c(5000, 20000)  # additional thresholds for the sensitivity analysis

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260908)
inp <- function(...) file.path(in_dir,  ...)
out <- function(...) file.path(out_dir, ...)

## ---- 1. capture history -----------------------------------------------------
CH0 <- readRDS(inp("capthist.rds"))
stopifnot(sum(CH0) == 345, nrow(CH0) == 252, sum(usage(traps(CH0))) == 4834)

cap0 <- read.delim(inp("secr_captures.txt"), header = FALSE, stringsAsFactors = FALSE)
names(cap0) <- c("Session", "ID", "Occasion", "Detector")
trap_xy <- as.data.frame(traps(CH0)); trap_xy$Detector <- rownames(traps(CH0))
stopifnot(nrow(cap0) == 345, !anyNA(match(cap0$Detector, trap_xy$Detector)))

## ---- 2. far-recapture exclusion rule ----------------------------------------
# Where the detections assigned to one identity span more than `threshold`,
# retain only the largest spatially connected group. Components come from
# single-linkage clustering. Ties on record count are resolved on the number of
# distinct cells, then the lowest component label, so the rule is deterministic.
# On these data the tie-break is immaterial: tied components hold one record
# each, so retained individuals, detections and pooled spatial variance are
# unchanged under random tie-breaking.
apply_far_rule <- function(cap, threshold) {
  keep <- rep(TRUE, nrow(cap))
  for (id in unique(cap$ID)) {
    ii <- which(cap$ID == id)
    xy <- trap_xy[match(cap$Detector[ii], trap_xy$Detector), , drop = FALSE]
    if (nrow(xy) < 2) next
    d <- dist(xy[, c("x", "y")])
    if (all(as.numeric(d) <= threshold)) next
    labs <- cutree(hclust(d, method = "single"), h = threshold)
    tab  <- sort(table(labs), decreasing = TRUE)
    top  <- names(tab)[tab == tab[1]]
    if (length(top) > 1) {
      ncell <- vapply(top, function(L)
        nrow(unique(xy[labs == as.integer(L), c("x", "y"), drop = FALSE])), integer(1))
      top <- top[order(-ncell, as.integer(top))]
    }
    keep[ii[labs != as.integer(top[1])]] <- FALSE
  }
  keep
}

build_CH <- function(threshold, tag) {
  keep <- apply_far_rule(cap0, threshold)
  f <- out(sprintf("capture_history_%s.txt", tag))
  write.table(cap0[keep, ], f, sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  ch <- read.capthist(f, inp("secr_traps.txt"), detector = "count",
                      fmt = "trapID", noccasions = 1, binary.usage = FALSE)
  covariates(ch) <- covariates(CH0)[match(rownames(ch), rownames(covariates(CH0))), ,
                                    drop = FALSE]
  ch
}

## ---- 3. covariates ----------------------------------------------------------
# Effort (transect points per cell) already enters through usage(), which makes
# the expected count proportional to it. c_logeff is the centred log of the same
# quantity, so a coefficient of zero means detections scale exactly with effort
# and a non-zero coefficient measures the departure.
# Unknown sex must be NA, not a third level, for hcov to integrate over it.
add_covariates <- function(ch) {
  e <- as.numeric(usage(traps(ch))[, 1])
  covariates(traps(ch)) <- data.frame(
    effort = e, c_logeff = as.numeric(scale(log(e), center = TRUE, scale = FALSE)))
  s <- as.character(covariates(ch)$Sex)
  s[s %in% c("U", "u", "")] <- NA_character_
  covariates(ch)$SexHC <- factor(s, levels = c("F", "M"))
  ch
}

CH <- add_covariates(build_CH(far_thr, "primary"))
stopifnot(!verify(CH, report = 0)$errors)

## ---- 4. state space ---------------------------------------------------------
# Full buffer, then clipped to forest within the suitable elevation band.
# subset() is used so the mask keeps its spacing and area attributes, which
# maskarea() and region.N() depend on.
mask_full <- make.mask(traps(CH), buffer = buffer_m, spacing = spacing_m,
                       type = "trapbuffer")

pts  <- terra::vect(cbind(mask_full$x, mask_full$y), type = "points",
                    crs = paste0("EPSG:", trap_epsg))
dem  <- terra::rast(file.path(layers_dir, "Bhutan_DEM.tif"))
if (!terra::same.crs(pts, dem)) pts <- terra::project(pts, terra::crs(dem))
elev <- as.numeric(terra::extract(dem, pts)[[2]])
# The DEM is clipped to the national boundary, so mask points outside Bhutan
# return NA and are excluded: the state space is habitat within Bhutan.

forest  <- st_make_valid(st_read(file.path(layers_dir, "forest_only.shp"), quiet = TRUE))
mask_sf <- st_as_sf(data.frame(x = mask_full$x, y = mask_full$y),
                    coords = c("x", "y"), crs = st_crs(trap_epsg))
if (st_crs(forest) != st_crs(trap_epsg)) mask_sf <- st_transform(mask_sf, st_crs(forest))
in_forest <- lengths(st_intersects(mask_sf, forest)) > 0

keep_hab <- !is.na(elev) & elev >= elev_range[1] & elev <= elev_range[2] & in_forest
mask_hab <- subset(mask_full, subset = keep_hab)

cat(sprintf("State space: %d points, %.0f km2 (%.0f%% of the %.0f km2 buffer)\n",
            nrow(mask_hab), maskarea(mask_hab) / 100,
            100 * nrow(mask_hab) / nrow(mask_full), maskarea(mask_full) / 100))

## ---- 5. model set -----------------------------------------------------------
# Density is constant throughout; the models differ in how the baseline
# detection probability varies. h2 is the two-class sex mixture from hcov, so
# every model estimates pmix. effort_free sets ignoreusage = TRUE, so effort
# acts only through the covariate; comparing it with `effort` tests whether the
# proportional-usage assumption is adequate.
specs <- list(
  null                = list(model = list(D ~ 1, g0 ~ 1,             sigma ~ 1),  iu = FALSE),
  sex                 = list(model = list(D ~ 1, g0 ~ h2,            sigma ~ 1),  iu = FALSE),
  effort              = list(model = list(D ~ 1, g0 ~ c_logeff,      sigma ~ 1),  iu = FALSE),
  sex_effort          = list(model = list(D ~ 1, g0 ~ h2 + c_logeff, sigma ~ 1),  iu = FALSE),
  sex_effort_sexsigma = list(model = list(D ~ 1, g0 ~ h2 + c_logeff, sigma ~ h2), iu = FALSE),
  effort_free         = list(model = list(D ~ 1, g0 ~ c_logeff,      sigma ~ 1),  iu = TRUE))

fit_model <- function(nm, ch, mask, tag, start) {
  f <- out(sprintf("%s_%s.rds", tag, nm))
  if (reuse && file.exists(f)) { cat("reusing", basename(f), "\n"); return(readRDS(f)) }
  cat("fitting", basename(f), "...\n")
  fit <- secr.fit(ch, mask = mask, detectfn = detfn, CL = FALSE, hcov = "SexHC",
                  binomN = 0, model = specs[[nm]]$model,
                  details = list(ignoreusage = specs[[nm]]$iu),
                  start = start, method = "BFGS", ncores = ncores, trace = FALSE)
  saveRDS(fit, f); fit
}

ini   <- autoini(CH, mask = mask_hab, detectfn = "HN")
start <- list(D = ini$D, g0 = ini$g0, sigma = ini$sigma)

fits <- list()
for (nm in names(specs))
  fits[[nm]] <- fit_model(nm, CH, mask_hab, "habitat",
                          if (length(fits)) fits[[1]] else start)

## ---- 6. model selection -----------------------------------------------------
model_table <- do.call(rbind, lapply(names(fits), function(nm) {
  a <- AIC(fits[[nm]])
  data.frame(model = nm,
             usage = if (specs[[nm]]$iu) "ignored" else "proportional",
             npar = a$npar[1], logLik = a$logLik[1], AIC = a$AIC[1])
}))
model_table <- model_table[order(model_table$AIC), ]
model_table$dAIC   <- model_table$AIC - min(model_table$AIC)
model_table$weight <- exp(-model_table$dAIC / 2) / sum(exp(-model_table$dAIC / 2))
print(model_table, row.names = FALSE, digits = 4)
write.csv(model_table, out("model_selection.csv"), row.names = FALSE)

# Likelihood-ratio tests for the hypotheses the model set was built to address.
LL <- setNames(model_table$logLik, model_table$model)
NP <- setNames(model_table$npar,   model_table$model)
lrt <- function(a, b, q) {
  chi <- 2 * (LL[[b]] - LL[[a]]); df <- NP[[b]] - NP[[a]]
  data.frame(question = q, chisq = chi, df = df,
             p = pchisq(chi, df, lower.tail = FALSE))
}
lrt_table <- rbind(
  lrt("null", "sex",        "g0 differs by sex"),
  lrt("null", "effort",     "detection departs from proportionality with effort"),
  lrt("sex_effort", "sex_effort_sexsigma", "sigma differs by sex"))
print(lrt_table, row.names = FALSE, digits = 3)
write.csv(lrt_table, out("likelihood_ratio_tests.csv"), row.names = FALSE)

# The minimum-AIC model discards the measured effort weighting. Among models
# within 2 AIC of the best, we report the most parsimonious one that retains it.
competitive <- model_table[model_table$dAIC < 2 & model_table$usage == "proportional", ]
primary_nm  <- competitive$model[order(competitive$npar, competitive$AIC)][1]
primary_fit <- fits[[primary_nm]]
cat("Reported model:", primary_nm,
    "| minimum-AIC model:", model_table$model[1], "\n")

## ---- 7. estimates -----------------------------------------------------------
estimates <- function(fit, mask, label) {
  pr <- predict(fit)
  pr <- if (is.data.frame(pr)) pr else pr[[1]]
  D <- pr["D", ]
  sigma <- pr["sigma", ]
  rn <- region.N(fit, mask = mask)
  data.frame(label = label,
             D_per_km2 = 100 * D[["estimate"]],
             LCL = 100 * D[["lcl"]], UCL = 100 * D[["ucl"]],
             CV_pct = 100 * D[["SE.estimate"]] / D[["estimate"]],
             sigma_m = sigma[["estimate"]],
             HR95_km2 = 18.8 * sigma[["estimate"]]^2 / 1e6,
             area_km2 = maskarea(mask) / 100,
             N = rn["E.N", "estimate"],
             N_lcl = rn["E.N", "lcl"], N_ucl = rn["E.N", "ucl"])
}

primary <- estimates(primary_fit, mask_hab, primary_nm)
print(primary, row.names = FALSE, digits = 4)
write.csv(primary, out("density_reported.csv"), row.names = FALSE)
print(predict(primary_fit), digits = 4)          # includes pmix, the sex ratio

# 95% home range from sigma assumes a bivariate-normal activity centre:
# area = pi * (sigma * sqrt(2 * log(20)))^2 = 18.8 * sigma^2.
cat(sprintf("sigma = %.0f m implies a 95%% home range of %.0f km2; published red panda\n",
            primary$sigma_m, primary$HR95_km2))
cat(sprintf("home ranges of 1-10 km2 correspond to sigma of %.0f-%.0f m.\n",
            sqrt(1e6 / 18.8), sqrt(1e7 / 18.8)))
cat(sprintf("buffer / sigma = %.1f (conventional guideline: 4)\n",
            buffer_m / primary$sigma_m))

## ---- 8. sensitivity to the exclusion threshold ------------------------------
# The reported model refitted at each threshold. This is the dominant source of
# uncertainty in the analysis and is not represented in the confidence interval
# above.
sens <- list(cbind(threshold_km = far_thr / 1000, detections = sum(CH),
                   estimates(primary_fit, mask_hab, primary_nm)))
for (thr in thresholds) {
  ch  <- add_covariates(build_CH(thr, paste0("thr", thr)))
  fit <- fit_model(primary_nm, ch, mask_hab, paste0("thr", thr), primary_fit)
  sens[[length(sens) + 1L]] <- cbind(threshold_km = thr / 1000, detections = sum(ch),
                                     estimates(fit, mask_hab, primary_nm))
}
sens <- do.call(rbind, sens)
sens <- sens[order(sens$threshold_km), ]
sens$sigma_over_threshold <- sens$sigma_m / (sens$threshold_km * 1000)
print(sens, row.names = FALSE, digits = 4)
write.csv(sens, out("threshold_sensitivity.csv"), row.names = FALSE)

cat(sprintf(paste0("\nAcross thresholds of %s km: density %.3f-%.3f km-2 (%.1f-fold),\n",
                   "sigma %.0f-%.0f m, and sigma/threshold %.2f-%.2f - a near-constant ratio,\n",
                   "indicating the detection scale is set by the exclusion rule rather than\n",
                   "by the animals. The three analyses differ by %d of %d detection records.\n"),
            paste(sens$threshold_km, collapse = ", "),
            min(sens$D_per_km2), max(sens$D_per_km2),
            max(sens$D_per_km2) / min(sens$D_per_km2),
            min(sens$sigma_m), max(sens$sigma_m),
            min(sens$sigma_over_threshold), max(sens$sigma_over_threshold),
            max(sens$detections) - min(sens$detections), nrow(cap0)))

## ---- 9. session -------------------------------------------------------------
sessionInfo()
