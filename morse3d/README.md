# Three-dimensional Morse cluster — Figure 1

Run `main` from `code_review/` for all calculations and automatic figure
export. To run only this task, use the command below.

For eight particle positions `x_i` in three dimensions, the Morse energy is

```math
E(q) = sum_{i<j} [exp(-2a(r_ij-r0)) - 2exp(-a(r_ij-r0))],
r_ij = ||x_i-x_j||, a=6, r0=1.
```

The experiment compares HiSD and SR-HiSD, searches for saddle states, and
connects states by gradient descent. The physical Hessian index excludes
three translations and three rotations. HiSD uses four ascent directions
and SR-HiSD uses one physical direction in the index-one comparison.

## Reproduce Figure 1

From `code_review/`, run in MATLAB:

```matlab
run_task('morse3d');  % full calculation and automatic plotting
```

| Subtask | Calculation |
|---|---|
| `compare_methods` | Run HiSD and SR-HiSD from identical perturbations |
| `landscape` | Search upward and follow both gradient-descent branches |
| `plot_figure1` | Draw the landscape and particle configurations |

For separate steps, use the same RunId:

```matlab
options=struct('RunId','morse3d_run');
run_task('morse3d','compare_methods',options);
run_task('morse3d','landscape',options);
run_task('morse3d','plot_figure1',options);
```

To redraw saved numerical results:

```matlab
run_task('morse3d','all',struct('Mode','reference'));
```

Use a new RunId when repeating a calculation. `Quick=true` with compute mode
runs a short execution check rather than the complete landscape.
The scripts can also be run directly from `morse3d/scripts/`.

## Results

Outputs are written to `morse3d/results/<RunId>/`:

| Output | Content |
|---|---|
| `compare_methods/compare_rows.csv` | HiSD and SR-HiSD comparison records |
| `landscape/result.mat` | State configurations, energies, indices and transitions |
| `landscape/{landscape_nodes.csv,landscape_edges.csv}` | State and transition tables |
| `plot_figure1/{figure1.png,figure1.pdf}` | Landscape corresponding to Figure 1 |
| `plot_figure1/node_images/node_<id>_<label>.png` | Individual configurations of all 17 states, with index and energy labels |

The full landscape contains 17 states and 106 recorded transitions. State
labels, physical indices and six-decimal energies are checked against
`data/paper/figure1_states.csv`. Accepted states have gradient norm <1e-5
and returned residual <5e-8. Calculation logs are stored in `morse3d/logs/`.
