## Southern Highbush Blueberry GxE with Factor Analytic Model

Multi-environment trial (MET) analysis for yield and bloom-timing traits
across AL, FL, NC breeding-program materials, plus a bloom-date forecasting
service. 


### Layout

```
R/pipeline/         targets pipeline stages (stage00 ... stage09)
R/R/                helper functions imported by stages
python/service/     FastAPI deployment service for bloom forecasting
data/raw/           raw souce files
data/processed/     output, processed data
.github/workflows/  CI: renv restore, pipeline test
```

### Running the R pipeline

```r
# one-time setup
install.packages("renv")

####
file.exists('renv/activate.R')
source('renv/activate.R')
renv::project()
####

renv::restore()

# run the full pipeline
targets::tar_make()

# check the pipeline graph
targets::tar_visnetwork()
```

Each stage is defined in `R/pipeline/stageNN_*.R` and wired together in
`R/pipeline/_targets.R`. Stages are numbered in dependency order.

### Stages

```
- stage00_dta_ingestion: takes al, nc, & fl phenotypic and metadta
- stage01_pedigree_process: takes pedigree from multiple sources for processing
- stage02_marker_ingestion: ingest marker dta
- stage03_cleaning: clean pheno dta with metadata file - also clean names.
- stage04_bloom_harmonization: gets the dates for bloom5, bloom50, bloom100, bloom_duration, and auc
- stage05_per_trial_qc: plot basis h2, CV, and reliability with asreml and sommer
- stage06_chill_heat: chill hrs and GDD
- stage07_fa_model: fa model for yield with asreml and sommer, includes diagnostic layer with post fa FAST (Smith & Cullis, 2018), produce plots.
- stage07b_fa_bloom_model: fa model and fast for bloom traits
- stage07c_fa_auc_model: fa model and fast for bloom auc.
- stage08_bloom_rxn_norm: chill/heat accummulation-based reaction norm.
- stage09_validation: leave-one-year-out validation for the forecasting service, mimics deployment scenario.
```


### Status
```
- 10/11/2026: FA models completed with pedigree data (Stages 00 - 07 completed).
- Still to complete: FA models with markers, move FW, AMMI, GGE into pipeline, reaction norm with chilling hrs, and validate forcasting before finishing with service deployment.
```

```r
# 1. from repo root, install deps once
renv::restore()

# 2. see the graph before running anything (good sanity check first)
targets::tar_visnetwork()

# 3. run everything
targets::tar_make()

# 4. read out a specific target's value once it's built
targets::tar_read(prep_yield) # prepared clean dta that goes into fa model
targets::tar_read(al_bundle)  # the raw loaded list, if you want to poke at MetaData/OffID directly
targets::tar_read(al_metadata_raw) # metadata

# 5. quick sanity checks specific to AL dta
bm <- targets::tar_read(al_bloom_long_clean)
nrow(bm)
range(bm$obs_date)
table(bm$location, bm$season)

# 6. see what's stale / what failed, without re-running
targets::tar_progress()
targets::tar_meta(fields = error) |> dplyr::filter(!is.na(error))
```
