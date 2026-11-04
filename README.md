# Induction Motor Design Optimization (IM_Opt)

A MATLAB framework for **maximizing the rated efficiency** of a 5.5 kW three-phase squirrel-cage induction motor. The analytical motor model follows the Boldea design methodology; the problem has **16 design variables** and **10 constraints**. Ten metaheuristic algorithms are compared under eight constraint-handling (penalty) approaches.

This repository contains the **result-producing code only** (model, constraints, algorithms, penalty methods and the driver). Reporting/plotting/tabulation scripts are not included; each run saves raw results as `.mat` files under `sonuclar/`.

## Features
- **10 algorithms:** PSO, TSA, DSA, WOA, GWO, JAYA, FOA (red fox), LFD (Lévy flight distribution), ABC, AAA
- **8 constraint-handling methods:** Deb's feasibility rule + static/dynamic penalty, sum/squared violation, normalized count, Joines–Houck, Liu (linear/squared)
- Single objective: efficiency maximization with feasibility-first selection
- 30 independent runs; per-run best designs and iteration-wise convergence/parameter histories are recorded

## Requirements
MATLAB (R2018b or newer recommended). The code also runs under Octave.

## Quick start
```matlab
% In the repository folder:
run_ALL          % 10 algorithms x 8 methods; saves .mat results to sonuclar/
```
To run a single algorithm/method:
```matlab
% method = 1..8  (1 = Deb's rule)
[eff, x, summary] = motorJaya(30, 40, 16, 0.2, 3, 100000);
```
`motorXXX(maxrun, pop, D, ST, method, max_fes)` is the common signature for all algorithms.

## Repository structure
| File | Purpose |
|------|---------|
| `im_core.m` | Analytical motor model: `[eff, gviol, out] = im_core(x)` (objective + 10 constraints) |
| `im_bounds.m` | Lower/upper bounds of the 16 design variables |
| `penalty.m`, `deb.m` | Eight penalty methods and Deb's selection rule |
| `motor*.m` | Drivers for the 10 metaheuristic algorithms (PSO, TSA/CTSA, DSA, WOA, GWO, JAYA, FOA, LFD, ABC, AAA) |
| `run_ALL.m` | Batch driver: all algorithms × all methods |
| `save_results.m` | Saves per-run results and convergence/parameter histories as `.mat` |
| `im_better.m`, `im_firstbest.m`, `im_bestupd.m`, `im_queue.m`, `im_score.m` | Shared selection / best-update / history helpers |
| `esson_constant.m`, `BmHm.m`, `Kst_Kf.m`, `Kst_Alfa_i.m`, `power_factor.m`, `efficiency_lookup.m`, `conductor_diameter.m` | Table/curve interpolations (range-checked lookups) |

## Design variables (x, 1×16)
Order: `Lamda_coef, KD_coef, Bg_coef, Kst, Jcos_coef, bos[mm], hos[mm], hw[mm], Bts[T], Jb[A/mm^2], Btr[T], hor[mm], bor[mm], rotor_slot[deg], Bcr[T], k`

> Note: internal variable names in the source are in Turkish (e.g., `Lamda_katsayi`, `Jcos_katsayi`, `verim` = efficiency, `ihlal/viol` = constraint violation); the mapping above and the English comments give their meaning.

## Citation
This code accompanies an academic study. If you use it, please cite the related paper. *(Add the full reference / DOI once published.)*

## License
MIT License — see `LICENSE`.
