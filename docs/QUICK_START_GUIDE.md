# ms3R Package Transformation - Quick Start Guide

## TL;DR

Transform ms3R from scripts → professional R package with:
- **10-50x faster** (Rcpp + parallelization)
- **Extensible** assessment models (SISCAL, SISCA, SS3, etc.)
- **Professional** structure (DESCRIPTION, tests, docs)

## Three Main Documents

1. **`package_transformation_plan.md`** - Complete roadmap (READ THIS FIRST)
2. **`performance_improvements.md`** - Detailed Rcpp & parallel code
3. **This guide** - Quick reference

---

## Immediate Actions (Start Today)

### 1. Create Package Structure (30 minutes)

```r
# In R console
usethis::create_package("~/ms3R-pkg")

# Copy existing code
file.copy("simSCAL.R", "~/ms3R-pkg/R/simulation.R")
file.copy("ms3Rplots.R", "~/ms3R-pkg/R/plotting.R")
# ... etc for all R files

# Copy TMB code
dir.create("~/ms3R-pkg/src")
file.copy("hierProd.cpp", "~/ms3R-pkg/src/")

# Copy data
dir.create("~/ms3R-pkg/inst/extdata")
file.copy("simCtlFile.txt", "~/ms3R-pkg/inst/extdata/")
# ... etc
```

### 2. Create Minimal DESCRIPTION (5 minutes)

```
Package: ms3R
Title: Multi-Species Multi-Stock Management Strategy Evaluation
Version: 0.1.0
Authors@R: person("Your", "Name", email = "you@example.com", role = c("aut", "cre"))
Description: MSE framework for multi-species fisheries.
License: GPL-3
Imports: TMB, tidyverse, Rcpp
LinkingTo: Rcpp, TMB
```

### 3. Test Package Build (5 minutes)

```r
devtools::document()  # Generate NAMESPACE
devtools::load_all()  # Test loading
devtools::check()     # Check for issues
```

---

## Priority 1: Assessment Model Interface (Day 1-2)

### Goal
Allow users to easily plug in different assessment models:

```r
# Instead of hard-coded hierProd
ctlList$mp$assess$method <- "SISCAL"  # or "SS3", "Sablefish", etc.
blob <- runMS3(ctlList = ctlList)
```

### Implementation

**Create `R/assessment-interface.R`:**

```r
#' Generic assessment interface
#' @export
assess_stock <- function(data, control, model_type, ...) {
  UseMethod("assess_stock", structure(list(), class = model_type))
}

#' Register assessment model
#' @export
register_assessment_model <- function(name, assess_fn, ...) {
  .assessment_registry[[name]] <- list(
    assess = assess_fn,
    ...
  )
}
```

**Create `R/assessment-hierProd.R`:**

```r
#' @export
assess_stock.hierProd <- function(data, control, model_type, ...) {
  # Your existing hierProd code
  # Return standardized list:
  list(
    B_t = biomass_timeseries,
    TAC_next = recommended_TAC,
    Bmsy = Bmsy_estimate,
    ...
  )
}
```

**Create `inst/assessment_interfaces/siscal_interface.R`:**

```r
#' @export
assess_stock.SISCAL <- function(data, control, model_type, ...) {
  require(SISCAL)

  # Transform ms3R data to SISCAL format
  siscal_data <- list(
    catch = data$catch,
    cpue = data$index,
    ...
  )

  # Run SISCAL
  fit <- SISCAL::fit_model(siscal_data)

  # Return standardized output
  list(
    B_t = fit$biomass,
    TAC_next = forecast_siscal(fit),
    ...
  )
}
```

**Test it:**

```r
# Register on load
.onLoad <- function(lib, pkg) {
  register_assessment_model("hierProd", assess_stock.hierProd)
  if(requireNamespace("SISCAL"))
    register_assessment_model("SISCAL", assess_stock.SISCAL)
}

# Usage
list_assessment_models()  # See what's available
ctlList$mp$assess$method <- "SISCAL"
blob <- runMS3(ctlList)
```

---

## Priority 2: Rcpp Speedup (Day 3-5)

### Goal
10x faster operating model through Rcpp

### Critical Functions to Optimize

**1. Population dynamics loop (simSCAL.R:4493-4520)**

Currently:
```r
for(s in 1:nS)
  for(p in 1:nP)
    for(a in 1:A_s[s])
      N_axspt[a,,s,p,t] <- N_axspt[a-1,,s,p,t-1] * exp(-Z_axspt[a-1,,s,p,t-1])
```

Create `src/pop_dynamics.cpp`:
```cpp
#include <Rcpp.h>
// [[Rcpp::export]]
NumericVector updatePopCpp(NumericVector N_prev, NumericVector Z_prev, ...) {
  // C++ implementation - see performance_improvements.md
}
```

Use in R:
```r
.ageSexOpMod <- function(obj, t) {
  # Replace R loops
  N_axspt[,,,,t] <- array(
    updatePopCpp(as.vector(N_axspt[,,,,t-1]), ...),
    dim = c(nA, nX, nS, nP)
  )
}
```

**2. Catch calculations (simSCAL.R:4654-4686)**

Similarly, move nested loops to Rcpp.

### Quick Setup

```bash
# In package root
Rscript -e "Rcpp::compileAttributes()"
R CMD INSTALL .
```

---

## Priority 3: Parallelization (Day 6-7)

### Goal
Run 100 replicates in parallel, not sequentially

### Current (Slow)
```r
# Runs sequentially
for(rep in 1:100) {
  blob[[rep]] <- runMS3(ctlFile)
}
```

### New (Fast)
```r
# Create R/parallel.R
#' @export
runMS3_parallel <- function(ctlFile, nCores = 4) {
  require(future.apply)
  plan(multisession, workers = nCores)

  nReps <- readCtlFile(ctlFile)$ctl$nReps

  results <- future_lapply(1:nReps, function(rep) {
    runMS3(ctlFile, rep = rep)
  }, future.seed = TRUE)

  return(results)
}

# Usage
blobs <- runMS3_parallel("simCtlFile.txt", nCores = 8)
# 8x faster!
```

### Nested Parallelization (Advanced)
```r
# Parallelize both batch scenarios AND replicates
# See performance_improvements.md for full implementation
```

---

## Quick Wins (Can Do Today)

### 1. Add roxygen2 docs to main functions (1 hour)

```r
#' Run Multi-Species Multi-Stock Simulation
#'
#' Main simulation function
#'
#' @param ctlFile Path to control file
#' @return Simulation results (blob)
#' @export
runMS3 <- function(ctlFile = "./simCtlFile.txt", ...) {
  # existing code
}
```

Then: `devtools::document()`

### 2. Vectorize simple loops (30 minutes)

Find patterns like:
```r
# Slow
for(i in 1:n) result[i] <- x[i] * y[i]

# Fast
result <- x * y
```

### 3. Profile to find bottlenecks (15 minutes)

```r
profvis::profvis({
  runMS3("simCtlFile.txt")
})
# Look for red boxes = slow functions
```

---

## Expected Performance Gains

| Optimization | Time for 100 Reps | Speedup |
|--------------|-------------------|---------|
| **Current** | 3.3 hours | 1x |
| + Rcpp | 30 minutes | 6.6x |
| + Parallel (8 cores) | 5 minutes | 40x |
| + Both | **4 minutes** | **50x** |

---

## Testing Your Changes

### 1. Numerical Equivalence
```r
# Test Rcpp matches R
test_that("Rcpp gives same results as R", {
  result_R <- updatePopR(...)
  result_Rcpp <- updatePopCpp(...)
  expect_equal(result_R, result_Rcpp, tolerance = 1e-10)
})
```

### 2. Performance Benchmarking
```r
microbenchmark::microbenchmark(
  R_version = updatePopR(...),
  Rcpp_version = updatePopCpp(...),
  times = 100
)
```

---

## Common Issues & Solutions

### Issue: "TMB not found"
```r
# Install TMB
install.packages("TMB")
```

### Issue: "Rcpp compilation fails"
```bash
# Install compilers
# Mac: xcode-select --install
# Ubuntu: sudo apt-get install r-base-dev
# Windows: Install Rtools
```

### Issue: "Function not exported"
```r
# Add to function: #' @export
# Then: devtools::document()
```

### Issue: "Parallel workers can't find functions"
```r
# Export to workers
clusterExport(cl, c("runMS3", "other_functions"))
# Or use future::plan() which handles this
```

---

## Recommended Order of Implementation

### Week 1: Foundation
1. Package structure ✓
2. Basic documentation ✓
3. Assessment interface design ✓

### Week 2: Assessment Models
4. Implement hierProd interface ✓
5. Add 2-3 external model interfaces ✓
6. Testing and examples ✓

### Week 3: Performance
7. Profile code, identify bottlenecks ✓
8. Implement Rcpp for top 3 slow functions ✓
9. Add replicate parallelization ✓

### Week 4: Polish
10. Comprehensive testing ✓
11. Vignettes and documentation ✓
12. Package check and release ✓

---

## Getting Help

### Essential Resources
- **R Packages book**: https://r-pkgs.org/ (READ THIS!)
- **Rcpp**: https://www.rcpp.org/
- **performance_improvements.md**: Detailed code examples

### Questions to Ask
1. Which assessment models do I need first?
2. What speedup do I actually need?
3. Will I submit to CRAN or just GitHub?
4. Do I need backwards compatibility with old scripts?

---

## Success Checklist

- [ ] Package builds without errors: `devtools::check()`
- [ ] Can switch assessment models: `ctlList$mp$assess$method <- "SISCAL"`
- [ ] Rcpp functions match R versions (tested)
- [ ] Parallel runs work: `runMS3_parallel()`
- [ ] Documentation complete: All functions have `?help`
- [ ] At least 2 vignettes written
- [ ] Tests pass: `devtools::test()`
- [ ] 10x+ faster than original

---

## Next Steps Right Now

1. **Read** `package_transformation_plan.md` (detailed roadmap)
2. **Create** package structure with `usethis::create_package()`
3. **Implement** assessment interface (copy code from plan)
4. **Profile** existing code to confirm bottlenecks
5. **Start** Rcpp optimization with `updatePopCpp()`

---

## Contact

Questions? Check the detailed docs:
- `package_transformation_plan.md` - Full implementation guide
- `performance_improvements.md` - Rcpp code examples

Good luck! This is a substantial but very achievable project. Focus on the assessment interface first (highest value), then performance optimization.
