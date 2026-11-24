# Two-dimensional Morse clusters — Table 1

Run `main` from `code_review/` for all calculations and automatic figure
export. To run only this task, use the command below.

For particle positions `x_i` in two dimensions, the Morse energy is

```math
E(q) = sum_{i<j} [exp(-2a(r_ij-r0)) - 2exp(-a(r_ij-r0))],
r_ij = ||x_i-x_j||, a=6, r0=1.
```

The experiment measures index-one saddle success rates for N=5–10 using
HiSD with 1–4 directions, SR-HiSD and NPSS at thresholds 0.10 and 0.01.
Each method uses the same 500 initial states per cluster size, for 21,000
searches. Translations and rotations are removed when counting the physical
Hessian index.

## Reproduce Table 1

From `code_review/`, run in MATLAB:

```matlab
run_task('morse2d');  % full calculation and automatic plotting
```

| Subtask | Calculation |
|---|---|
| `original_methods` | Run HiSD with 1–4 directions and SR-HiSD |
| `npss` | Run NPSS at both thresholds |
| `plot_table1` | Combine success rates and export Table 1 |

For separate steps, use the same RunId:

```matlab
options=struct('RunId','morse2d_run');
run_task('morse2d','original_methods',options);
run_task('morse2d','npss',options);
run_task('morse2d','plot_table1',options);
```

To redraw saved numerical results:

```matlab
run_task('morse2d','all',struct('Mode','reference'));
```

Use a new RunId when repeating a calculation. `Quick=true` with compute mode
runs a short execution check using two N=5 initial states and 50 steps.
The scripts can also be run directly from `morse2d/scripts/`.

## Results

Outputs are written to `morse2d/results/<RunId>/`:

| Output | Content |
|---|---|
| `original_methods/{result.mat,trials.csv,rates.csv}` | HiSD and SR-HiSD endpoints, trial diagnostics and rates |
| `npss/{result.mat,trials.csv,rates.csv}` | NPSS endpoints, trial diagnostics and rates |
| `plot_table1/{table1.csv,table1.png,table1.pdf}` | Seven method rows and 42 success percentages corresponding to Table 1 |

The full table is checked against `data/paper/table1.csv`. A successful
endpoint has RMS gradient ≤1e-8, gradient norm ≤1e-5, a connected cluster
and exactly one negative physical Hessian eigenvalue below -1e-6.
Calculation logs are stored in `morse2d/logs/`.
