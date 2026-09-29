# MPM_3D_MATLAB_Hydro-Mechanical_Coupling

> A **3D Material Point Method (MPM)** numerical simulation platform in MATLAB for **long-term
> performance analysis of pavement subgrade structures**, featuring dynamic analysis, contact,
> unsaturated seepage, and **hydro-mechanical coupling** under cyclic traffic loading.

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021b+-orange.svg)](https://www.mathworks.com/)

---

## 📖 Overview

During long-term service, **gradual moisture accumulation** in the subgrade causes stiffness
degradation and sustained accumulation of permanent deformation, eventually leading to
pavement cracking. This problem involves:

- **Climate-driven** unsaturated seepage (rainfall infiltration / evaporation);
- **Cyclic traffic loading** and elastic–plastic strain accumulation;
- **Hydro-mechanical coupling** between the moisture field and the stress field.

Conventional FEM software (e.g. ABAQUS) suffers from mesh distortion and "pseudo-hang"
issues (process alive but no result output) beyond 10k–20k load cycles, making long-term
simulation of subgrade performance infeasible.

This project is a **deeply re-engineered fork** of the open-source program
**[AMPLE-MATLAB](https://github.com/wmcoombs/AMPLE)**.
The original 6-module framework has been extended into a dedicated
**PAVMPM platform** for hydro-mechanically coupled subgrade analysis, with
**3,000+ lines of new code** and **>80% of the codebase modified**.

```
Original platform (AMPLE-MATLAB)          This platform (PAVMPM)
├── Dynamic Solver                        ├── ✅ Fully vectorized framework (efficiency +100%)
├── Contact                               ├── ✅ Hydro-mechanical bounding-surface model (UMAT → MPM)
├── Seepage                               ├── ✅ Explicit–implicit adaptive hybrid solver
├── Coupling                              └── ✅ Long-term cyclic-loading capability
├── Interpolation
└── Axisymmetric
```

---

## 📁 Repository Structure

```
MPM_3D_MATLAB_Hydro-Mechanical_Coupling/
├── setup/                          # Model setup and initialization
├── functions/                      # Core MPM solver (original framework)
├── function_BX/                    # ⭐ Vectorized hydro-mechanical coupling modules (_BX)
│                                     #   (shape functions, constitutive integration, MP update)
├── constitutive/                   # ⭐ Constitutive models
│                                     #   (elastic, elasto-plastic, bounding-surface UMATs)
├── function_axial/                 # Axisymmetric analysis utilities
├── plotting/                       # Post-processing and visualization
├── ample_singleSolve_Flow_3Dslope_BX.m   # ⭐ Main entry: 3D slope, coupled flow-mechanics
└── 三维高程数据.txt                  # 3D elevation data (terrain geometry)
```

> Files/modules suffixed with **`_BX`** are the **vectorized hydro-mechanical coupling
> versions** developed in this work; the other folders largely preserve the original
> AMPLE-MATLAB architecture.

---

## 🔧 Key Improvements

### 1. Full Vectorization — one-pass solution for all material points

**Problem.** The original code loops over material points (`for i = 1:N`) to compute shape
functions, strains, constitutive updates and stiffness assembly one by one — prohibitively
slow at ~100k material points.

**Solution.** All point-wise loops were re-written as **high-dimensional matrix operations**
covering all material points in a single pass:

| Data structure | Original (per-point loop) | Vectorized |
|---|---|---|
| Constitutive matrix [D] | `[6,6]` | `[6,6,N]` |
| Stress | `[6,1]` | `[6,1,N]` |
| Strain | `[6,1]` | `[6,1,N]` |
| Shape-function index | point-by-point | `[x,1,N]` batch fetch |

Core techniques:

- **Elastic/elasto-plastic encoder matrix**: points are flagged with `0/1` so that
  `D[:,:,encoder==1] = D_ep` and `D[:,:,encoder==0] = D_e` — no per-point branching;
- **`pagemtimes` batch matrix multiplication**: stress update
  `STRESS = STRESS + pagemtimes(D, dStrain)` is executed page-wise in parallel, **no `for` loops**.

**Benchmarks** (MATLAB Profiler, 3 steps):

| Stage | Original | Vectorized | Speedup |
|---|---|---|---|
| Internal force & [K] assembly | 19.0 s | 7.0 s | **+173%** |
| Constitutive integration | 12.5 s | 3.3 s | **+280%** |
| Time per step | 23 s | 11 s | **Overall >100%** |

### 2. Hydro-Mechanical Bounding-Surface Model (UMAT → MPM)

**Problem.** The built-in constitutive models cannot describe the behavior of unsaturated
subgrade soils under **varying moisture content, stress state and compaction degree**.

**Solution.** The group's validated unsaturated **bounding-surface plasticity model** was
ported into the MPM framework with full vectorization:

- **Reversible deformation (resilient modulus)** — semi-empirical model explicitly
  accounting for moisture, compaction and stress state:

  E = k₀·pₐ·(p*/pₐ)^k₁·(qᵣ/pₐ + 1)^k²·(q/pₐ + 1)^k³·(K̃^k₄)

  Compared with NCHRP 1-28A (R = 0.45), the new model reaches **R = 0.97**;

- **Irreversible deformation (permanent deformation)** — bounding-surface model with
  loading / memory / bounding surfaces, capturing **OCR effects** and **cyclic loading effects**;
- **Hysteresis** — main wetting/drying lines and scanning curves of the soil-water
  characteristic curve (SWCC).

**Verification.** Permanent-deformation curves from the MPM platform agree closely with
ABAQUS UMAT + FEM results across different moisture contents and compaction degrees
(OMC, PYL130, and 87–96% compaction cases).

### 3. Explicit–Implicit Adaptive Hybrid Solver

**Problem.**
- The **implicit solver** uses large time steps and suits long-term analysis but **fails to converge**;
- The **explicit solver** needs no stiffness matrix and is unconditionally stable, but its small
  time step makes one load cycle ~700× slower than implicit.

**Solution.** An **implicit-first with explicit fallback** adaptive switching strategy:

```
 ┌─────────┐
 │ Initial │
 └────┬────┘
      ▼
 ┌──────────────┐   divergence criterion:
 │  Implicit    │   |f_int + f_dct − f_ext| / f_ext > 1/1000  or no convergence ──► ┌──────────┐
 │  Solver      │                                                                    │ Explicit │
 └────┬────┘                                                                    │  Solver  │
      │                                                                          └────┬───┘
      │ N                                                                              │
      ◄──────────────────── recovery criterion: ◄──────── max(v) < 0.01 m/s for 50 steps ── Y
```

- **Implicit → Explicit**: switch when the residual exceeds 1/1000 or convergence fails;
- **Explicit → Implicit**: switch back once the system has sufficiently calmed
  (max velocity < 0.01 m/s for 50 consecutive steps).

**Results.** Displacement histories from both solvers agree within **< 0.01%** error;
the hybrid solver combines implicit efficiency with explicit robustness, pushing the
analysis through **tens of thousands of load cycles**.

### 4. Climate-Driven Unsaturated Seepage

Complete unsaturated flow capability was added:

- **van Genuchten (VG) model** for the SWCC and **VG–Mualem** permeability function,
  with void-ratio dependence;
- **FAO Penman–Monteith** model for potential evapotranspiration ET₀;
- **Climate boundary**: rainfall infiltration / evaporation flux;
- **Groundwater**: water-head boundary with capillary rise;
- Two-way coupling with the mechanical field through **effective stress**:

  p* = (p_t − u_a) + Sᵣ·s,   s = f(σ, Sᵣ, e)

---

## 🧪 Example Cases

### Coupled response of a pavement structure with/without a sand bed

- Model: 8 m subgrade, half cross-section of a bidirectional four-lane road;
- Loading: 100 kN dual-wheel groups, 1 cycle/min, 1,440 cycles;
- Boundary: daily climate cycle of 6 h rainfall + 18 h evaporation, groundwater at h = 3 m;
- Comparison: **Case A** (no sand bed) vs **Case B** (pervious sand bed, k = 1×10⁻⁵ m/s).

| Output | Key finding |
|---|---|
| Seepage field (suction / head change) | The sand bed markedly accelerates moisture redistribution and drainage |
| Plastic strain (volumetric / shear) | Differences concentrate near the bottom of the pavement structure |
| Displacement history | Permanent deformation accumulates steadily with load cycles; the two cases are clearly distinguishable |

These cases demonstrate that the platform captures **climate-driven seepage evolution**
and **cyclic mechanical response** simultaneously, supporting the selection of design
subgrade moisture states.

---

## 🚀 Getting Started

```matlab
% 1. Add the repository to the MATLAB path
addpath(genpath('MPM_3D_MATLAB_Hydro-Mechanical_Coupling'));

% 2. Configure model parameters in setup/ (soil properties, climate, load history, terrain)

% 3. Run the main script (implicit-first, explicit-fallback switching is automatic)
ample_singleSolve_Flow_3Dslope_BX;
```

**Requirements:** MATLAB R2021b or newer (relies on `pagemtimes` and other tensor intrinsics).

---

## 🗺️ Roadmap

- [ ] **Physics-Informed Neural Networks (PINN)** — GCN-Transformer-ResNet architecture to
      accelerate MPM computations and predict nodal displacement/stress time series;
- [ ] Additional benchmark cases for accuracy verification;
- [ ] Coupling with the national-scale **equilibrium subgrade moisture content** prediction
      (GeoStudio seepage + GEP empirical formula, R² = 0.91) to provide the design worst-case
      moisture state as model input.

---

## 📚 References

1. Dunatunga S., Kamrin K. (2017). Continuum modeling and simulation of granular flows
   through their many phases. *Physical Review Fluids*. — original AMPLE-MATLAB.
2. Fan H., et al. (2024). An extended resilient modulus model considering moisture,
   compaction and stress state for subgrade soils.
3. Zhang et al. (2025). Semi-empirical resilient modulus model of unsaturated subgrade soils.
4. NCHRP Report 1-28A (2004). Laboratory determination of resilient modulus for flexible
   pavement design.
5. van Genuchten M.T. (1980). A closed-form equation for predicting the hydraulic
   conductivity of unsaturated soils.
6. Allen R.G., et al. (1998). FAO Penman-Monteith equation for computing crop evapotranspiration.
---

## ✉️ Contact

- **Author**: FAN Haishan
- **Affiliation**: Department of Civil and Environmental Engineering,
  The Hong Kong Polytechnic University
- **Email**: haishan.fan@polyu.edu.hk
