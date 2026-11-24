# Landau–Brazovskii — Figure 3

The zero-mean field phi is periodic on [0,L)^3, with mean energy

```math
E(phi) = mean[ xi^2/2*((Delta+q0^2)phi)^2
             + tau/2*phi^2 - gamma/6*phi^3 + phi^4/24 ].
```

LBSD uses the original direction space. SR-HiSD projects out translations
and uses shared FFT/force caches with one generalized-Gram Ritz update per
outer iteration. Both methods use semi-implicit position steps and BB2.

The original experiments used MATLAB R2023a on shared, CPU-only nodes with
two 2.90 GHz Intel Xeon Gold 6326 processors (16 physical cores each).
Each Slurm job used four physical cores, 16,000 MB RAM and four MATLAB
computational threads. Reported search times are elapsed wall-clock times
for the stated scope. Different MATLAB versions and hardware may lead to
small differences in iteration counts and search times.

## Reproduce

Run `main` from `code_review/` for all experiments and automatic plotting.
For the LB task alone:

```matlab
run_task('lb');
```

| Case | Method | Calculation |
|---|---|---|
| BCC → HEX | LBSD | Exactly 4 → 3 → 2 → 1; no direct 4 → 1 route |
| BCC → HEX | SR-HiSD | FFT-cached `sr1490067`, exact supplied input |
| BCC → disorder | LBSD | Native four-direction search, with index-2 descent if needed |
| BCC → disorder | SR-HiSD | FFT-cached `sr32410`, exact supplied input |

All cases use N=128, xi=q0=1, RMS tolerance 1e-7 and at most 80,000 iterations
per search attempt. SR uses four threads, seed 20261028, block tolerance 1e-8,
tangent tolerance 1e-6 and a 43,200-second search budget.

| Case | (tau, gamma) | L | Initial step | BB range | Kick |
|---|---|---|---:|---|---:|
| HEX | (-0.3, 0.8) | 24*sqrt(2)*pi | 0.24 | [0.2, 2] | 0.01 |
| Disorder | (0.05, 0.6) | 32*sqrt(2)*pi | 0.1 | [0.02, 0.2] | 0.02 |

HEX LBSD reads the saved `phi_init` and four-column `V_init` from
`data/input/bcc2hex_lbsd.mat` without rebuilding the frame. It follows
4 → 3 → 2 → 1 with a starting downward kick of 0.001. The kick can increase
to 0.002 or 0.004 only when both signs are ineffective or return to the
verified parent. Each accepted stage uses independent warm/cold original
mean-zero Hessian spectra with 12 Ritz pairs and a ±1e-7 zero band.
Fresh SR endpoints require an independently checked raw-Hessian negative
direction before bilateral GF.
This one-direction check does not establish the complete Morse index.
HEX endpoints must also pass sixfold Fourier geometry and independent
same-parameter phase matching.

## Steps and results

| Subtask | Outputs in `results/<RunId>/<subtask>/` |
|---|---|
| `bcc2hex` | Fields, attempts, timing, `transiti1.png/pdf` and available nucleus plots |
| `bcc2disorder` | Fields, attempts, timing, `transiti2.png/pdf` and available nucleus plots |
| `nuclei` | Measurements and individual PNG/PDF images from same-parent GF endpoints |
| `plot_figure3` | `figure3.png` and `figure3.pdf` combining both cases |

Each case runs and plots independently. To run separate steps, keep the same
RunId, for example `run_task('lb','bcc2hex',struct('RunId','lb_run'))`.
Calculation logs are in `logs/`; SR history and checkpoints accompany the
computed fields. Recorded four-case timings and original job logs are in
`logs/imported_20261120/` and `logs/imported_hex1490068_20261120/`.
The selected four-case summary is `logs/timing.csv` and `logs/timing.json`.
Search, kernel, verification and process-wall times have separate columns.
SR bilateral GF timings are in `logs/new_sr_gf_20261120_timing.csv`.

The recorded HEX LBSD index-4 upward search took 26,987 iterations and
56,449.427268 s (15.68 h) of solver time, with 14,418.061355 s of separate
spectral diagnostics. Its sequential descent took 9,134 and 8,856 iterations,
then 80,000 unsuccessful iterations on the positive final branch and 1,729
on the successful negative branch. The complete recorded sequential route
used 126,706 iterations. These are supplied run measurements; 15.68 h is
not the full route or a runtime prediction for another machine.

To draw supplied saved states:

```matlab
run_task('lb','all',struct('Mode','reference'));
```

Saved results include both SR states and their own bilateral GF endpoints in
`data/reference/<case>_srhisd/gf/`. Figure 3 contains all four same-parent
nucleus outlines. Full computation runs new GF from each new endpoint.
A short execution check is `run_task('lb','all',struct('Quick',true))`.

## Plot colors

Images and colorbars use the 256-color linear `lb_phase` map:

| Position | RGB |
|---|---|
| 0.00 | (0.04, 0.10, 0.36) |
| 0.45 | (0.35, 0.78, 1.00) |
| 0.55 | (0.25, 0.68, 0.18) |
| 1.00 | (1.00, 0.94, 0.05) |

HEX uses range [-1,1] and ticks [-1,0,1]; disorder uses [-0.4,0.8] and ticks
[-0.4,0,0.4,0.8]. Values are clipped to the endpoint colors. LBSD arrows are
blue (#0075AB), SR-HiSD red (#BD354A), and GF gray (#71777B). White nucleus
ellipses use the same saddle slice and its matching GF endpoint fields.
