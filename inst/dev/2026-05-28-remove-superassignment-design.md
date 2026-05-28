# Design: Remove `<<-` Parent-Scope Assignment

**Date:** 2026-05-28  
**Branch strategy:** Two PRs — `refactor/superassignment-pt1` (utilities/boostmtree/generic.predict), then `refactor/superassignment-pt2` (vimp) after pt1 merges  
**Baseline:** main @ `51c94cc`, 144 tests passing  
**Tracking:** vault `Projects/boostmtree.md` item 4

---

## Problem

~30 `<<-` sites across 4 files assign into parent scopes from inside `lapply` closures. This makes data flow non-obvious, complicates reasoning about state, and is flagged in the 2026-03-24 code review. All sites fall into four patterns, each with a clean functional replacement.

---

## Pattern catalogue

| Pattern | Files | Sites | Replacement |
|---------|-------|-------|-------------|
| Counter increment | utilities.R | 1 | `seq_along` / direct index — eliminate counter |
| List-building | boostmtree.R, generic.predict, vimp | ~8 | `x <- lapply(...)` returning the value |
| Running-sum accumulation | generic.predict, vimp | ~10 | `Reduce("+", lapply(...))` |
| Matrix population | vimp | ~9 | `do.call(rbind, lapply(...))` or `sapply` |

No `for` loops. All replacements use functional R idioms.

---

## PR #3 — `refactor/superassignment-pt1`

Files: `utilities.R`, `boostmtree.R`, `generic.predict.boostmtree.R`  
Sites: 11 lines (1 + 5 + 5 active; line 618 of generic.predict is commented-out dead code, excluded)

### utilities.R (1 site, line 625)

**Current:** `count <<- count + 1` inside `lapply(which.null, ...)` that merges failed mclapply results back into `result`. The counter tracks position in `result.lapply`.

**Replacement:** `which(which.null == i)` gives the position directly — no counter needed.

```r
# before
result <- lapply(X, function(i) {
  if (any(i == which.null)) {
    count <<- count + 1
    result.lapply[[count]]
  } else {
    result.mclapply[[i]]
  }
})

# after
result <- lapply(X, function(i) {
  if (any(i == which.null)) {
    result.lapply[[which(which.null == i)]]
  } else {
    result.mclapply[[i]]
  }
})
```

### boostmtree.R (4 sites, lines 1208, 1400, 1402, 1404–1405)

**Site 1 — gamma.i.list (line 1208):** `gamma.i.list[[q]][[m]][[i]] <<- gamma.matx.i` inside an `lapply(1:n, ...)` that currently discards its return value.

Replacement: return `gamma.matx.i` from the closure; assign into `gamma.i.list[[q]][[m]]` outside.

```r
# before
lapply(1:n, function(i) {
  ...
  gamma.i.list[[q]][[m]][[i]] <<- gamma.matx.i
  ...
})

# after
gamma.i.list[[q]][[m]] <- lapply(1:n, function(i) {
  ...
  gamma.matx.i
})
```

**Sites 2–4 — cv.flag block (lines 1400–1405):** `Mopt[q] <<-`, `rmse[q] <<-`, `mu[[q]] <<-` inside `lapply(1:n.Q, function(q) { ...; NULL })`.

Replacement: return a named list from the closure; extract vectors/list outside.

```r
# before
nullObj <- lapply(1:n.Q, function(q) {
  ...
  Mopt[q] <<- min(which(diff.err < tol))   # or M
  rmse[q] <<- err.rate[[q]][Mopt[q], "l2"]
  mu[[q]] <<- lapply(1:n, function(i) mu.cv.list[[q]][[Mopt[q]]][[i]])
  NULL
})

# after
cv.results <- lapply(1:n.Q, function(q) {
  ...
  mopt.q <- if (sum(diff.err < tol) > 0) min(which(diff.err < tol)) else M
  list(
    Mopt = mopt.q,
    rmse = err.rate[[q]][mopt.q, "l2"],
    mu   = lapply(1:n, function(i) mu.cv.list[[q]][[mopt.q]][[i]])
  )
})
Mopt <- sapply(cv.results, `[[`, "Mopt")
rmse <- sapply(cv.results, `[[`, "rmse")
mu   <- lapply(cv.results, `[[`, "mu")
```

### generic.predict.boostmtree.R (6 active sites, lines 526–536, 603)

**Sites 1–2 — beta accumulation (lines 526, 528):** `beta <<- beta.m` / `beta <<- beta + beta.m` inside `lapply(1:Mopt[q], ...)`. The `m==1` case initialises beta to `beta.m` (adding `Ymean` to column 1); subsequent iterations accumulate.

Replacement: compute each `beta.m` in a lapply, then `Reduce("+", beta_m_list)`. The `Ymean` adjustment for `m==1` applies only to column 1 of `beta.m` — add it to the first element of the list before reducing.

```r
# after
beta_m_list <- lapply(1:Mopt[q], function(m) {
  ...
  bm <- t(beta.m.org * nu.vec * Ysd)
  if (m == 1) bm[, 1] <- bm[, 1] + Ymean
  bm
})
beta <- Reduce("+", beta_m_list)
```

**Sites 3–4 — mu.list building (lines 532, 536):** `mu.list[[m]] <<- lapply(1:n, ...)` inside the same outer lapply.

Replacement: return the per-m mu list from the closure; `mu.list` becomes the collected result of `lapply(1:Mopt[q], ...)`.

```r
# after
mu.list <- lapply(1:Mopt[q], function(m) {
  ...
  lapply(1:n, function(i) {
    ...   # mu values for iteration m
  })
})
```

**Site 5 — l_pred_db_hat_Temp cache (line 603):** `l_pred_db_hat_Temp[[q]][[i]][[m]] <<- l_pred_db_Temp` inside a `Reduce("+", lapply(1:Mopt[q], ...))`. The write is a side-effect cache update for the Ordinal monotonicity constraint — the value is also the return for the Reduce.

Replacement: build `l_pred_db_hat_Temp[[q]][[i]]` as the full list from an outer `lapply(1:Mopt[q], ...)`, then pass it into a separate `Reduce` step. This separates the cache-building from the summation.

```r
# after
l_pred_db_hat_Temp[[q]][[i]] <- lapply(1:Mopt[q], function(m) {
  l_pred_db_Temp <- ...
  if (family == "Ordinal" && q > 1) {
    l_pred_db_Temp <- ifelse(
      l_pred_db_Temp < l_pred_db_hat_Temp[[q - 1]][[i]][[m]],
      l_pred_db_hat_Temp[[q - 1]][[i]][[m]],
      l_pred_db_Temp
    )
  }
  l_pred_db_Temp
})
sum_l_pred_db_Temp <- Reduce("+", l_pred_db_hat_Temp[[q]][[i]][1:Mopt[q]])
```

---

## PR #4 — `refactor/superassignment-pt2`

Cut from main after PR #3 merges.  
File: `vimp.boostmtree.R`  
Sites: 17

### oob.list (line 217)

`oob.list[[q]][[m]] <<- oob` inside `lapply(1:Mopt[q], function(m) { ...; oob })`.

Replacement: `oob.list[[q]] <- lapply(1:Mopt[q], function(m) { ...; oob })`

### l_pred_db accumulation (lines 273–276, 422, 572–578, 719 — two code paths)

Running sums `l_pred_db.main.i <<- l_pred_db.main.i + out.main` etc. across `lapply(1:Mopt[q], ...)`.

Replacement: `Reduce` over a list of per-iteration contribution lists:

```r
# after (per subject i, per variable k)
contributions <- lapply(1:Mopt[q], function(m) {
  ...
  list(main = out.main, int = out.int, time = if (k == p) out.time else NULL)
})
l_pred_db.main.i <- Reduce("+", lapply(contributions, `[[`, "main"))
l_pred_db.int.i  <- Reduce("+", lapply(contributions, `[[`, "int"))
if (k == p)
  l_pred_db.time.i <- Reduce("+", Filter(Negate(is.null), lapply(contributions, `[[`, "time")))
```

### vimp matrix population (lines 375, 384, 394, 460, 678, 687, 697, 758 — two code paths)

`vimp.main[k, q] <<- ...`, `vimp.int[k, q] <<- ...`, `vimp.time[q] <<- ...` inside nested `lapply(1:n.Q)(1:p)` returning NULL.

Replacement: return a named list per (k, q) cell; extract into matrices outside:

```r
# after
raw <- lapply(1:n.Q, function(q) {
  lapply(1:p, function(k) {
    ...
    list(main = (err.rate.main - rmse[q]) / rmse[q],
         int  = (err.rate.int  - rmse[q]) / rmse[q],
         time = if (k == p) (err.rate.time - rmse[q]) / rmse[q] else NA_real_)
  })
})
vimp.main <- matrix(sapply(unlist(raw, recursive=FALSE), `[[`, "main"), nrow=p)
vimp.int  <- matrix(sapply(unlist(raw, recursive=FALSE), `[[`, "int"),  nrow=p)
vimp.time <- sapply(raw, function(q_res) q_res[[p]][["time"]])
```

---

## Verification strategy

Each file gets a standalone equivalence script (`/tmp/equiv_<file>.R`) run before editing:

1. Extract the old closure logic into `old_fn()`
2. Express the new functional logic in `new_fn()`
3. Run both on `simLong` data across multiple seeds
4. Assert `identical()` on all outputs

Commit sequence within each PR:
1. Write equivalence script, confirm green
2. Apply edit to source file
3. Run full suite (must be 144+ pass / 0 fail)
4. Commit with equivalence note in message

---

## Out of scope

- Commented-out `<<-` at `generic.predict.boostmtree.R:618` — dead code, leave as-is
- Performance benchmarking — correctness only; `Reduce`/`lapply` are not slower than the current `lapply(...NULL)` pattern
- Any other refactoring in files touched

---

## Success criteria

- Zero `<<-` remaining in R/ (verifiable with `grep -r "<<-" R/`)
- Suite: ≥144 pass / 0 fail after each commit
- Each equivalence script confirms `identical()` before source edit
