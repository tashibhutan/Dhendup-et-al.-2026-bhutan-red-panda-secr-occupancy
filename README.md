# Red panda SECR and occupancy analyses

This repository contains reproducible R scripts and prepared input data for analyses of red pandas in Bhutan.

The SECR analysis uses a prepared count-detector capture history, detector coordinates and effort, individual sex covariates, and a genetic ID-to-detector audit. Habitat masking uses a clipped Bhutan DEM and forest-cover shapefile. The script rebuilds capture histories, applies the specified far-recapture rule, fits half-normal SECR models with sex and effort alternatives, estimates density and abundance, and evaluates sensitivity to the exclusion threshold.

The occupancy analysis uses a red panda-only, camera-by-occasion dataset containing detection and non-detection histories, seven-day sampling effort, season, detection covariates, occupancy covariates, and projected camera coordinates. A small elevation table supports interpretation of the quadratic elevation response. The script fits spatial single-species occupancy models with spOccupancy, using an exponential NNGP process. It compares detection and occupancy structures with WAIC, evaluates alternative spatial-range priors, summarizes occupancy and detection probabilities, assesses the elevation response, and performs Freeman–Tukey posterior predictive checks.

Both scripts are simplified versions of the analyses used for the manuscript and are intended to be run from their respective folders with the required R packages installed. No fitted model outputs are included in this repository.
