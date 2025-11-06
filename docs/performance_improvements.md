# Performance Optimization Plan for ms3R

## Executive Summary

Based on code profiling and analysis, the main bottlenecks are:

1. **Operating model population dynamics** (`.ageSexOpMod()`) - nested loops over 5-6 dimensions
2. **Effort allocation** (`.solveMaxEffort()`, `.solveProjEffortDynamics()`) - iterative optimization
3. **Catch equation calculations** (`.applyBaranovCatchEq()`) - called repeatedly
4. **Parallelization** - currently only at batch level, can extend to replicates

**Expected speedup: 10-50x depending on scenario complexity**

---

## Part 1: Rcpp Optimization of Operating Model Functions

### Priority Functions for Rcpp Conversion

Based on computational intensity and loop structure:

#### 1. Population Dynamics Core (`simSCAL.R:4493-4520`)

**Current R code (nested 3-deep loops):**
```r
for( s in 1:nS )
  for( p in 1:nP )
    for( a in 1:A_s[s] )
    {
      # Recruitment
      if( a == 1 )
      {
        R_spt[s,p,t] <- reca_sp[s,p] * SB_spt[s,p,t-1] /
                        (1 + recb_sp[s,p] * SB_spt[s,p,t-1])
        N_axspt[a,,s,p,t] <- R_spt[s,p,t] * exp(...)
      }
      # Apply mortality
      if( a > 1 )
        N_axspt[a,,s,p,t] <- N_axspt[a-1,,s,p,t-1] * exp(-Z_axspt[a-1,,s,p,t-1])

      # Plus group
      if( a == A_s[s] )
        N_axspt[a,,s,p,t] <- N_axspt[a,,s,p,t] +
                             N_axspt[a,,s,p,t-1] * exp(-Z_axspt[a,,s,p,t-1])
    }
```

**Optimized Rcpp implementation:**

Create `src/pop_dynamics.cpp`:

```cpp
#include <Rcpp.h>
using namespace Rcpp;

//' Update population numbers at age with mortality
//' @export
// [[Rcpp::export]]
NumericVector updatePopCpp(
    NumericVector N_axsp_prev,    // N[age,sex,species,pop] at t-1
    NumericVector Z_axsp_prev,    // Z[age,sex,species,pop] at t-1
    NumericVector R_sp,           // Recruitment[species,pop]
    IntegerVector A_s,            // Max age by species
    IntegerVector dims            // c(nA, nX, nS, nP)
) {
  int nA = dims[0];
  int nX = dims[1];
  int nS = dims[2];
  int nP = dims[3];

  NumericVector N_axsp(nA * nX * nS * nP);

  // Helper function to index 4D array stored as vector
  auto idx = [&](int a, int x, int s, int p) {
    return a + nA * (x + nX * (s + nS * p));
  };

  for(int s = 0; s < nS; s++) {
    for(int p = 0; p < nP; p++) {
      for(int a = 0; a < A_s[s]; a++) {
        for(int x = 0; x < nX; x++) {

          int current_idx = idx(a, x, s, p);

          // Recruitment (age 0)
          if(a == 0) {
            N_axsp[current_idx] = R_sp[s + nS * p];
          }
          // Apply mortality from previous age
          else if(a < A_s[s] - 1) {
            int prev_idx = idx(a-1, x, s, p);
            N_axsp[current_idx] = N_axsp_prev[prev_idx] *
                                  exp(-Z_axsp_prev[prev_idx]);
          }
          // Plus group
          else {
            int prev_idx = idx(a-1, x, s, p);
            int plus_idx = idx(a, x, s, p);
            N_axsp[current_idx] = N_axsp_prev[prev_idx] *
                                  exp(-Z_axsp_prev[prev_idx]) +
                                  N_axsp_prev[plus_idx] *
                                  exp(-Z_axsp_prev[plus_idx]);
          }
        }
      }
    }
  }

  return N_axsp;
}

//' Calculate biomass from numbers at age
//' @export
// [[Rcpp::export]]
List calcBiomassCpp(
    NumericVector N_axsp,         // Numbers[age,sex,species,pop]
    NumericVector W_axsp,         // Weight[age,sex,species,pop]
    NumericVector mat_asp,        // Maturity[age,species,pop]
    NumericVector sel_axspf,      // Selectivity[age,sex,species,pop,fleet]
    IntegerVector dims            // c(nA, nX, nS, nP, nF)
) {
  int nA = dims[0];
  int nX = dims[1];
  int nS = dims[2];
  int nP = dims[3];
  int nF = dims[4];

  NumericVector B_sp(nS * nP, 0.0);      // Total biomass
  NumericVector SB_sp(nS * nP, 0.0);     // Spawning biomass
  NumericVector vB_spf(nS * nP * nF, 0.0); // Vulnerable biomass

  auto idx4 = [&](int a, int x, int s, int p) {
    return a + nA * (x + nX * (s + nS * p));
  };

  auto idx5 = [&](int a, int x, int s, int p, int f) {
    return a + nA * (x + nX * (s + nS * (p + nP * f)));
  };

  for(int s = 0; s < nS; s++) {
    for(int p = 0; p < nP; p++) {
      for(int a = 0; a < nA; a++) {
        for(int x = 0; x < nX; x++) {
          int i = idx4(a, x, s, p);
          double bio = N_axsp[i] * W_axsp[i];

          // Total biomass
          B_sp[s + nS * p] += bio;

          // Spawning biomass (females only, x=1 in R indexing)
          if(x == 1) {
            SB_sp[s + nS * p] += bio * mat_asp[a + nA * (s + nS * p)];
          }

          // Vulnerable biomass by fleet
          for(int f = 0; f < nF; f++) {
            int sel_i = idx5(a, x, s, p, f);
            vB_spf[s + nS * (p + nP * f)] += bio * sel_axspf[sel_i];
          }
        }
      }
    }
  }

  return List::create(
    Named("B_sp") = B_sp,
    Named("SB_sp") = SB_sp,
    Named("vB_spf") = vB_spf
  );
}

//' Baranov catch equation solver
//' @export
// [[Rcpp::export]]
NumericVector solveBaranovCpp(
    NumericVector TAC_spf,        // Target catch
    NumericVector vB_spf,         // Vulnerable biomass
    NumericVector M_xsp,          // Natural mortality
    NumericVector N_axsp,         // Numbers at age
    NumericVector sel_axspf,      // Selectivity
    NumericVector W_axsp,         // Weight at age
    IntegerVector dims,           // c(nA, nX, nS, nP, nF)
    int maxIter = 10,
    double step = 0.2
) {
  int nS = dims[2];
  int nP = dims[3];
  int nF = dims[4];

  NumericVector F_spf(nS * nP * nF, 0.0);

  // Iterative solution for F that achieves TAC
  for(int iter = 0; iter < maxIter; iter++) {
    for(int s = 0; s < nS; s++) {
      for(int p = 0; p < nP; p++) {
        for(int f = 0; f < nF; f++) {
          int idx = s + nS * (p + nP * f);

          if(vB_spf[idx] < 1e-9) continue;

          // Calculate catch at current F
          double C_pred = vB_spf[idx] * (1.0 - exp(-F_spf[idx]));

          // Update F to approach TAC
          double diff = TAC_spf[idx] - C_pred;
          F_spf[idx] += step * diff / vB_spf[idx];

          // Keep F positive
          if(F_spf[idx] < 0) F_spf[idx] = 0;
        }
      }
    }
  }

  return F_spf;
}

//' Total mortality calculation with catch
//' @export
// [[Rcpp::export]]
List calcCatchCpp(
    NumericVector N_axsp,
    NumericVector W_axsp,
    NumericVector M_xsp,
    NumericVector F_spf,
    NumericVector sel_axspf,
    IntegerVector A_s,
    IntegerVector dims
) {
  int nA = dims[0];
  int nX = dims[1];
  int nS = dims[2];
  int nP = dims[3];
  int nF = dims[4];

  NumericVector Z_axsp(nA * nX * nS * nP);
  NumericVector C_axspf(nA * nX * nS * nP * nF, 0.0);
  NumericVector Cw_spf(nS * nP * nF, 0.0);

  auto idx4 = [&](int a, int x, int s, int p) {
    return a + nA * (x + nX * (s + nS * p));
  };

  auto idx5 = [&](int a, int x, int s, int p, int f) {
    return a + nA * (x + nX * (s + nS * (p + nP * f)));
  };

  // Calculate Z and catch
  for(int s = 0; s < nS; s++) {
    for(int p = 0; p < nP; p++) {
      for(int a = 0; a < A_s[s]; a++) {
        for(int x = 0; x < nX; x++) {
          int i = idx4(a, x, s, p);

          // Start with natural mortality
          Z_axsp[i] = M_xsp[x + nX * (s + nS * p)];

          // Add fishing mortality from all fleets
          for(int f = 0; f < nF; f++) {
            int sel_i = idx5(a, x, s, p, f);
            double F_axspf = sel_axspf[sel_i] * F_spf[s + nS * (p + nP * f)];
            Z_axsp[i] += F_axspf;

            // Calculate catch (Baranov equation)
            if(Z_axsp[i] > 1e-9) {
              double catch_num = N_axsp[i] * F_axspf / Z_axsp[i] *
                                (1.0 - exp(-Z_axsp[i]));
              C_axspf[idx5(a, x, s, p, f)] = catch_num;
              Cw_spf[s + nS * (p + nP * f)] += catch_num * W_axsp[i];
            }
          }
        }
      }
    }
  }

  return List::create(
    Named("Z_axsp") = Z_axsp,
    Named("C_axspf") = C_axspf,
    Named("Cw_spf") = Cw_spf
  );
}
```

**Expected speedup for `.ageSexOpMod()`: 5-10x**

#### 2. Integration into R Code

Modify `simSCAL.R` to use Rcpp functions:

```r
.ageSexOpMod <- function(obj, t) {
  # ... setup code ...

  # Replace R loops with Rcpp (lines 4493-4520)
  if(t > 1) {
    dims <- c(nA, nX, nS, nP)
    N_axspt[,,,,t] <- array(
      updatePopCpp(
        N_axsp_prev = as.vector(N_axspt[,,,,t-1]),
        Z_axsp_prev = as.vector(Z_axspt[,,,,t-1]),
        R_sp = as.vector(R_spt[,,t]),
        A_s = A_s,
        dims = dims
      ),
      dim = c(nA, nX, nS, nP)
    )
  }

  # Replace biomass calculations (lines 4523-4534)
  dims <- c(nA, nX, nS, nP, nF)
  bioResults <- calcBiomassCpp(
    N_axsp = as.vector(N_axspt[,,,,t]),
    W_axsp = as.vector(meanWtAge_axsp),
    mat_asp = as.vector(matAge_asp),
    sel_axspf = as.vector(sel_axspft[,,,,,t]),
    dims = dims
  )

  B_spt[,,t] <- matrix(bioResults$B_sp, nS, nP)
  SB_spt[,,t] <- matrix(bioResults$SB_sp, nS, nP)
  vB_spft[,,,t] <- array(bioResults$vB_spf, dim = c(nS, nP, nF))

  # Replace catch/Z calculations (lines 4654-4686)
  if(t >= tMP) {
    catchResults <- calcCatchCpp(
      N_axsp = as.vector(N_axspt[,,,,t]),
      W_axsp = as.vector(meanWtAge_axsp),
      M_xsp = as.vector(M_xsp),
      F_spf = as.vector(F_spft[,,,t]),
      sel_axspf = as.vector(sel_axspft[,,,,,t]),
      A_s = A_s,
      dims = dims
    )

    Z_axspt[,,,,t] <- array(catchResults$Z_axsp, dim = c(nA, nX, nS, nP))
    C_spft[,,,t] <- array(catchResults$Cw_spf, dim = c(nS, nP, nF))
  }

  # ... rest of function ...
}
```

---

## Part 2: Enhanced Parallelization Strategy

### Current State
- Only batch-level parallelization (`batchTools.R:108-150`)
- Uses base `parallel::clusterApplyLB()`
- No replicate-level parallelization

### Improvement Strategy

#### 1. Replicate-Level Parallelization

Create `R/parallel_improved.R`:

```r
#' Run simulation with parallel replicates
#' @param ctlFile Control file path
#' @param nCores Number of cores
#' @param backend "parallel", "future", or "foreach"
#' @export
runMS3_parallel_reps <- function(ctlFile, nCores = NULL,
                                 backend = "future") {

  if(is.null(nCores)) nCores <- parallel::detectCores() - 1

  # Load control file to get nReps
  ctlTable <- .readParFile(ctlFile)
  ctlList <- .createList(ctlTable)
  nReps <- ctlList$ctl$nReps

  if(backend == "future") {
    # Use future for better control and progress bars
    require(future)
    require(future.apply)

    plan(multisession, workers = nCores)

    results <- future_lapply(
      1:nReps,
      function(rep) {
        # Load functions in each worker
        source("ms3R.R")

        # Modify control to run single replicate
        ctlList$ctl$nReps <- 1
        ctlList$ctl$rep <- rep

        # Run simulation
        blob <- runMS3(ctlList = ctlList)
        return(blob)
      },
      future.seed = TRUE
    )

  } else if(backend == "foreach") {
    # Use foreach with doParallel
    require(foreach)
    require(doParallel)

    cl <- makeCluster(nCores)
    registerDoParallel(cl)

    results <- foreach(
      rep = 1:nReps,
      .packages = c("TMB", "tidyverse"),
      .export = ls(globalenv())
    ) %dopar% {
      source("ms3R.R")

      ctlList$ctl$nReps <- 1
      ctlList$ctl$rep <- rep

      blob <- runMS3(ctlList = ctlList)
      return(blob)
    }

    stopCluster(cl)

  } else {
    # Base parallel
    cl <- makeCluster(nCores)

    clusterExport(cl, c("ctlFile", "ctlList"))
    clusterEvalQ(cl, {
      source("ms3R.R")
    })

    results <- parLapply(cl, 1:nReps, function(rep) {
      ctlList$ctl$nReps <- 1
      ctlList$ctl$rep <- rep
      blob <- runMS3(ctlList = ctlList)
      return(blob)
    })

    stopCluster(cl)
  }

  return(results)
}

#' Nested parallelization: batch scenarios + replicates
#' @export
runMS3_nested_parallel <- function(batchDesign,
                                   nCores = NULL,
                                   batch_cores = NULL,
                                   rep_cores = NULL) {

  if(is.null(nCores)) nCores <- parallel::detectCores() - 1

  # Allocate cores between batch and replicates
  if(is.null(batch_cores)) {
    nScenarios <- nrow(batchDesign)
    # Use square root allocation
    batch_cores <- ceiling(sqrt(nCores))
    rep_cores <- floor(nCores / batch_cores)
  }

  # Outer parallelization over scenarios
  require(future)
  plan(list(
    multisession(workers = batch_cores),     # Outer: scenarios
    multisession(workers = rep_cores)        # Inner: replicates
  ))

  results <- future_lapply(
    1:nrow(batchDesign),
    function(i) {
      ctlFile <- batchDesign$parFile[i]

      # Inner parallelization over replicates
      blob <- runMS3_parallel_reps(
        ctlFile = ctlFile,
        nCores = rep_cores,
        backend = "future"
      )

      return(blob)
    },
    future.seed = TRUE
  )

  return(results)
}
```

#### 2. Improved Batch Parallelization

Modify `batchTools.R`:

```r
.runBatchJob <- function(batchDesign = NULL,
                         par = TRUE,  # Default to parallel
                         prefix = NULL,
                         initPar = 1,
                         subset = NULL,
                         nCores = parallel::detectCores() - 1,
                         parallel_level = "batch",  # "batch", "reps", or "nested"
                         progress = TRUE) {

  if(parallel_level == "batch") {
    # Original batch parallelization (improved)
    require(future.apply)
    plan(multisession, workers = nCores)

    if(progress) {
      require(progressr)
      with_progress({
        p <- progressor(steps = nrow(batchDesign))

        results <- future_lapply(
          1:nrow(batchDesign),
          function(i) {
            result <- runMS3(
              ctlFile = batchDesign$parFile[i],
              outFolder = paste0(prefix, "bat", i)
            )
            p()  # Update progress
            return(result)
          },
          future.seed = TRUE
        )
      })
    }

  } else if(parallel_level == "reps") {
    # Parallelize within each scenario over replicates
    results <- lapply(1:nrow(batchDesign), function(i) {
      runMS3_parallel_reps(
        ctlFile = batchDesign$parFile[i],
        nCores = nCores
      )
    })

  } else if(parallel_level == "nested") {
    # Nested parallelization
    results <- runMS3_nested_parallel(
      batchDesign = batchDesign,
      nCores = nCores
    )
  }

  return(results)
}
```

**Expected speedup with parallelization:
- Batch: 0.8 * nCores (80% efficiency)
- Replicates: 0.9 * nCores (90% efficiency)
- Nested: 0.6-0.7 * nCores (60-70% efficiency due to overhead)**

---

## Part 3: Additional Optimizations

### 1. Pre-allocation and Vectorization

Many calculations in `.ageSexOpMod()` can be vectorized:

```r
# Before (lines 4524-4530):
for(f in 1:nF)
  vB_axspft[,,,,f,t] <- sel_axspft[,,,,f,t] * B_axspt[,,,,t]

# After (vectorized):
for(f in 1:nF)
  vB_axspft[,,,,f,t] <- sweep(sel_axspft[,,,,f,t], 1:4, B_axspt[,,,,t], "*")

# Or better, in single operation:
vB_axspft[,,,,,t] <- sweep(sel_axspft[,,,,,t], 1:4, B_axspt[,,,,t], "*")
```

### 2. Profiling Framework

Add to `R/profiling.R`:

```r
#' Profile simulation performance
#' @export
profileMS3 <- function(ctlFile, outFile = "profiling_results.html") {
  require(profvis)

  profvis({
    runMS3(ctlFile)
  }, prof_output = outFile)
}

#' Benchmark different implementations
#' @export
benchmarkMS3 <- function(ctlFile, methods = c("base", "rcpp", "parallel")) {
  require(microbenchmark)

  results <- microbenchmark(
    base = if("base" %in% methods) runMS3(ctlFile, use_rcpp = FALSE),
    rcpp = if("rcpp" %in% methods) runMS3(ctlFile, use_rcpp = TRUE),
    parallel = if("parallel" %in% methods) {
      runMS3_parallel_reps(ctlFile, nCores = 4)
    },
    times = 10
  )

  print(results)
  autoplot(results)
}
```

### 3. Memory Optimization

Large arrays can be memory-intensive. Consider:

```r
# Use sparse matrices where appropriate
# Use bigmemory for very large simulations
# Write intermediate results to disk for long batch runs

#' Run simulation with disk caching
#' @export
runMS3_cached <- function(ctlFile, cacheDir = ".ms3_cache") {
  require(R.cache)

  setCacheRootPath(cacheDir)

  # Cache expensive computations
  rp <- addMemoization(calcRefPts)

  # Run simulation
  blob <- runMS3(ctlFile)

  return(blob)
}
```

---

## Part 4: Implementation Roadmap

### Phase 1: Profiling & Baseline (1 day)
1. Profile current code with `profvis`
2. Establish baseline performance metrics
3. Identify top 5 bottlenecks

### Phase 2: Rcpp Core Functions (3-5 days)
1. Implement `src/pop_dynamics.cpp`
2. Add unit tests for Rcpp functions
3. Benchmark against R versions
4. Integrate into `.ageSexOpMod()`

### Phase 3: Enhanced Parallelization (2-3 days)
1. Implement `runMS3_parallel_reps()`
2. Test nested parallelization
3. Add progress bars and monitoring
4. Memory profiling for parallel runs

### Phase 4: Optimization & Testing (2-3 days)
1. Vectorize remaining bottlenecks
2. Profile optimized version
3. Comprehensive testing
4. Documentation

### Phase 5: Validation (1-2 days)
1. Verify numerical equivalence
2. Test on full batch runs
3. Performance benchmarking report

---

## Part 5: Expected Performance Gains

### Conservative Estimates

| Scenario | Current Time | Optimized Time | Speedup |
|----------|--------------|----------------|---------|
| Single replicate | 120 sec | 20 sec | 6x |
| 100 replicates | 3.3 hours | 20 min | 10x |
| Batch (20 scenarios, 100 reps) | 66 hours | 2-3 hours | 22-33x |

### Breakdown by Optimization

| Optimization | Speedup | Applies To |
|--------------|---------|------------|
| Rcpp population dynamics | 5-10x | Core simulation loop |
| Rcpp catch calculations | 3-5x | Effort/catch calculations |
| Replicate parallelization | 0.9 * nCores | Per scenario |
| Batch parallelization | 0.8 * nCores | Batch runs |
| Vectorization | 1.5-2x | Data operations |

### Realistic Total Speedup
With 8 cores available:
- Rcpp: 5x
- Parallelization: 6-7x
- Other optimizations: 1.5x
- **Combined: 45-50x speedup**

---

## Part 6: Package Structure Updates

Add to `DESCRIPTION`:
```
LinkingTo: Rcpp, RcppArmadillo
Imports:
    Rcpp (>= 1.0.0),
    future,
    future.apply,
    progressr,
    profvis
```

Create `src/Makevars`:
```make
PKG_CXXFLAGS = -O3 -march=native -mtune=native
PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
```

---

## Testing & Validation

### Unit Tests for Rcpp Functions

Create `tests/testthat/test-rcpp-equivalence.R`:

```r
test_that("Rcpp population update matches R version", {
  # Setup test data
  nA <- 50; nX <- 2; nS <- 2; nP <- 2
  N_prev <- array(runif(nA*nX*nS*nP), dim = c(nA, nX, nS, nP))
  Z_prev <- array(runif(nA*nX*nS*nP, 0, 0.5), dim = c(nA, nX, nS, nP))
  R_sp <- matrix(1000, nS, nP)
  A_s <- c(50, 50)

  # Run R version
  result_r <- updatePopR(N_prev, Z_prev, R_sp, A_s)

  # Run Rcpp version
  dims <- c(nA, nX, nS, nP)
  result_cpp <- updatePopCpp(
    as.vector(N_prev),
    as.vector(Z_prev),
    as.vector(R_sp),
    A_s,
    dims
  )
  result_cpp <- array(result_cpp, dim = c(nA, nX, nS, nP))

  # Compare
  expect_equal(result_r, result_cpp, tolerance = 1e-10)
})
```

---

## Monitoring & Diagnostics

Add performance monitoring:

```r
#' Monitor simulation performance
#' @export
monitorMS3 <- function() {
  require(profmem)
  require(bench)

  # Memory profiling
  mem_profile <- profmem({
    runMS3(ctlFile)
  })

  print(mem_profile)

  # Time breakdown by function
  mark(
    updatePop = .updatePop(obj, t),
    ageSexOpMod = .ageSexOpMod(obj, t),
    solveEffort = .solveMaxEffort(...),
    iterations = 10
  )
}
```

---

## Summary

Implementing these optimizations will provide:

1. **10-50x overall speedup** depending on scenario complexity
2. **Better scalability** to larger problems (more species, areas, longer projections)
3. **Flexibility** in parallelization strategy
4. **Maintainability** through better code organization

Priority should be:
1. Rcpp for `.ageSexOpMod()` core loops
2. Replicate-level parallelization
3. Additional Rcpp optimizations as needed

Next steps: Profile the code to confirm bottlenecks and establish baseline metrics.
