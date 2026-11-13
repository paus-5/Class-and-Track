# Class-and-Track — Python port

Python port of the MATLAB code in this repository: OTU classification into AOB/NOB
(MILP/MIQP) and optimal tracking of OTU abundances with a nitrification chemostat model.
Everything uses free solvers (HiGHS, SCIP, IPOPT via CasADi).

## Setup

```bash
conda env create -f python/environment.yml      # creates env "class-and-track"
conda activate class-and-track
cd python
```

On Windows, if `python` crashes inside scipy when called directly, run each script via
`conda run --no-capture-output -n class-and-track python <script>.py`.

All scripts read `../Data`, `../Classification/MAT_files` and `../Tracking/MAT_files`, write
results to `results/` (pickles / npz) and figures to `figures/`.

## 1. Check the port against MATLAB (optional)

| Command | What it does |
|---|---|
| `python validate_against_matlab.py` | open-loop simulation and one ASRE step vs the saved MATLAB Reactor A runs |
| `python validate_classification.py` | re-solves the saved classification instances; same AOB/NOB split as MATLAB |

## 2. Reactor A: tracking with its own classification

```bash
python run_direct_reactor_A.py -1 0.5 1e-5        # direct method (collocation + IPOPT), ~10 min
bash run_all.sh                                    # ASRE runs for comparison (slow) + direct run
python plot_reactor_A.py -1 0.5 1e-5               # cost / OTU / control / substrate figures
```

Arguments are `vmin vmax lambda`: control bounds on v = u − 1 and the control weight
(Q = I, R = λI). Expected (days 184–313, v ∈ [−1, 0.5], λ = 1e-5): open-loop J ≈ 0.700,
best ASRE iterate J ≈ 0.095, direct J ≈ 0.047.

## 3. Reactor A tracked with Reactor B's classification

```bash
python classify_reactor.py B l2                    # MIQP (SCIP), ~1 min -> results/classification_B_l2.npz
python compare_classifications.py                  # OTU-by-OTU comparison with A's classification
python run_direct_reactor_A.py -1 0.5 1e-5 --class B   # ~3 min
python plot_reactor_A.py -1 0.5 1e-5 --class B
python compare_substrates.py                       # s1, s2, s3 fit: classification A vs B
```

Only the AOB/NOB labels and yields come from Reactor B; abundances, initial state, inputs
and kinetics are Reactor A's. Expected substrate RMSE (direct method): classification A
0.14 / 0.35 / 0.28, classification B 0.19 / 1.09 / 0.97 (s1 / s2 / s3).
`B.xls` starts at day 5; the Python code uses absolute days for B (MATLAB shifted B by 5 days).

## 4. Lotka–Volterra proof of concept (article, "Proof of concept")

```bash
python plot_lv_dynamics.py                         # synthetic gLV data vs the article's figure
python run_lv_poc.py [noise] [seed]                # simulate -> track -> identify A
```

Case study 1: 1 AOB + 1 NOB, μ = (0.77, 1.07), K = (0.7, 1.3), 1/y = (3.98, 16.13),
A = [[1, 0.1], [−1.5, 1]], D = 0.22, s_in 1.25 → 1.95 at day 150, 300 days, daily samples.
D = 0.22 reproduces the article's images (`tex/Images/proof_of_concept/250309_POC_*`);
the article text says 0.24 and its table says K₂ = 0.3, but neither matches the figures.

`run_lv_poc.py` then
1. tracks the synthetic biomasses with the direct method (R = diag(1e-4, 1e-5) as in the
   article; a few seconds),
2. identifies A by least squares on u(t) − 1 = A x(t) over days 5–280 (the article's
   two-step approach),
3. estimates A directly in the ODE by collocation, with several initial guesses (can take
   well over 10 minutes).

Expected without noise: step 2 gives A ≈ [[0.95, 0.00], [−1.50, 1.00]]. a₁₂ is not
recovered because x₂ ≤ 0.026, so a₁₂x₂ is below the error of the tracking control.

## Code layout

| File | Content |
|---|---|
| `class_and_track/model.py` | nitrification chemostat model, control v = u − 1 on each OTU's growth |
| `class_and_track/data.py` | data loading, classification inputs, tracking-problem setup |
| `class_and_track/classification.py` | AOB/NOB MILP (HiGHS) / MIQP (SCIP) |
| `class_and_track/asre.py` | ASRE tracking (port of the thesis method) and the cost evaluator |
| `class_and_track/direct.py` | direct method: Radau collocation + IPOPT (the problem is stiff; explicit RK4 fails) |
| `class_and_track/lv.py` | gLV model, interaction regression and direct estimation of A |
