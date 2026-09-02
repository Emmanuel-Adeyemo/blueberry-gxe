# Southern Highbush Blueberry GxE Modeling

Multi-environment trial (MET) analysis for yield and bloom-timing traits
across AL, FL, NC breeding-program materials, plus a bloom-date forecasting
service.

**Start here:** [`DESIGN.md`](./DESIGN.md) is the authoritative design
specification — every modeling decision, its rationale, and its literature
grounding live there. This README covers repo layout and how to run things;
DESIGN.md covers *why* the pipeline is shaped this way.

## Layout

```
R/pipeline/        targets-orchestrated pipeline stages (stage00 ... stage09)
R/R/                plain, unit-tested helper functions imported by stages
tests/testthat/     unit tests for R/R/
python/service/     FastAPI deployment service for Model 2 (bloom forecasting)
data/raw/           untouched source files, one subfolder per source
data/processed/     targets-generated intermediate/output data (gitignored)
artifacts/model2/   versioned, dated Model 2 fit objects served by the API
.github/workflows/  CI: renv restore, unit tests, pipeline smoke test
```

## Running the R pipeline

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

# inspect the pipeline graph
targets::tar_visnetwork()
```

Each stage is defined in `R/pipeline/stageNN_*.R` and wired together in
`R/pipeline/_targets.R`. Stages are numbered in dependency order — see
DESIGN.md §5 for what each stage does and why it's sequenced that way.

## Running the deployment service

```bash
cd python
pip install -e .
uvicorn service.main:app --reload
```

The service loads a versioned Model 2 artifact from `artifacts/model2/` and
serves bloom-date predictions given current-year weather input. It does not
retrain — see DESIGN.md §6.2 / §8 for the refresh workflow.

## Status

Pre-implementation scaffold. Stage files contain function signatures and
`TODO` markers referencing the relevant DESIGN.md section, not working
statistical code yet 



```r
# 1. from repo root, install deps once
renv::restore()

# 2. see the graph before running anything (good sanity check first)
targets::tar_visnetwork()

# 3. run everything
targets::tar_make()

# 4. read out a specific target's value once it's built
targets::tar_read(al_bloom_long)
targets::tar_read(al_bundle)          # the raw loaded list, if you want to poke at MetaData/OffID directly
targets::tar_read(al_metadata_raw)

# 5. quick sanity checks specific to AL
b <- targets::tar_read(al_bloom_long)
nrow(b)
range(b$obs_date)
table(b$location, b$season)

# 6. see what's stale / what failed, without re-running
targets::tar_progress()
targets::tar_meta(fields = error) |> dplyr::filter(!is.na(error))
```
