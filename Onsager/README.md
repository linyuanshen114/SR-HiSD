# Onsager — Figure 2

Run `main` from `code_review/` for all calculations and automatic figure
export. To run only this task, use the command below.

Compute the SR-HiSD pathway with physical indices **5 → 4 → 2 → 0**
and plot the four orientational densities.

The probability density rho is defined on the unit sphere with integral
rho = 1. On the 2,562-point grid,

```math
E(rho) = w*sum(rho_i*log(rho_i)) + w^2/2*sum(K_ij*rho_i*rho_j),
K_ij = 80*sqrt(1-(p_i.p_j)^2),  w = 4*pi/2562.
```

A constant-centered kernel sets the energy offset. SR-HiSD removes mass
and rotational orbit directions; both raw and projected residuals must be
at most 1e-7.

## Reproduce

From `code_review/`, run in MATLAB:

```matlab
run_task('Onsager');  % full calculation and automatic Figure 2 export
```

To redraw saved results, run
`run_task('Onsager','all',struct('Mode','reference'))`.
A short execution check uses `struct('Quick',true)` as the third argument.
A complete calculation requires base MATLAB with its JVM and about 4 GB RAM.

| Subtask | Result |
|---|---|
| `search` | SR-HiSD pathway, `result.mat` and `nodes.csv` |
| `plot_figure2` | `figure2.png`, `figure2.pdf` and four single-state images in `states/index{5,4,2,0}.png` |

For separate steps, use the same RunId:

```matlab
options=struct('RunId','onsager_run');
run_task('Onsager','search',options);
run_task('Onsager','plot_figure2',options);
```

Use a new RunId when repeating a calculation. `options.InputDir` can select
a saved run for plotting.

## Results

Outputs are in `Onsager/results/<RunId>/<subtask>/` and calculation logs
are in `Onsager/logs/`.

| Physical index | Energy shown in Figure 2 |
|---|---:|
| 5 | 28.648747 |
| 4 | 25.115151 |
| 2 | 23.384271 |
| 0 | 4.608069 |

The figure shows log(rho) on each sphere. Equivalent density orientations
are related by the model's rotational symmetry.
