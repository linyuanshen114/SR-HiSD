# SR-HiSD reproduction

Symmetry-reduced high-index saddle dynamics (SR-HiSD) searches for saddle
points of a prescribed Morse index in energy landscapes with continuous
symmetries. It formulates HiSD on the quotient space of symmetry-equivalent
states and lifts the dynamics to the original variables, with ascent
directions orthogonal to the symmetry orbits. This package reproduces the
Morse-cluster, Onsager and Landau–Brazovskii experiments in Table 1 and
Figures 1–3, including saddle searches, solution landscapes and transition
pathways.

## SR-HiSD dynamics

For a symmetry-invariant energy $E(x)$, let $P_x$ be the orthogonal projector
onto the horizontal space $\mathcal H_x=T_x(G\cdot x)^\perp$. With horizontal,
orthonormal directions $v_1,\ldots,v_k$, the index-$k$ SR-HiSD dynamics are

```math
\begin{aligned}
\dot x
&=-\left(I-2\sum_{j=1}^{k}v_jv_j^\top\right)\nabla E(x),\\
\dot v_i
&=-\left(I-v_iv_i^\top-2\sum_{j=1}^{i-1}v_jv_j^\top\right)
P_x\nabla^2E(x)P_xv_i
+DP_x[\dot x]v_i,\qquad i=1,\ldots,k.
\end{aligned}
```

Here $vv^\top u=v\langle v,u\rangle$ uses the model's inner product. The
position evolves uphill along the $k$ ascent directions and downhill in
their orthogonal complement. The reduced Hessian $P_x\nabla^2E(x)P_x$
excludes symmetry directions, while $DP_x[\dot x]v_i$ maintains the
horizontal constraint $P_xv_i=v_i$ as the state changes.

## Numerical method

The horizontal projector is constructed from infinitesimal symmetry
generators using tolerance-based Gram–Schmidt orthonormalization. These
generators represent rigid motions for Morse clusters, rotations for the
Onsager model and periodic translations for the LB model.

The general scheme combines a reflected-gradient Euler position step
with bounded Barzilai–Borwein (BB2) step sizes. The previous ascent frame
is projected onto the horizontal space at the updated state, and projected
LOBPCG/Rayleigh–Ritz updates approximate the $k$ smallest reduced-Hessian
modes while maintaining orthonormality. The LB implementation uses Fourier
pseudospectral discretization and semi-implicit position steps. Upward and
downward searches connect saddles of different indices; bilateral gradient
flows from index-1 saddles identify the connected phases or minima.
Model parameters and implementation details are given in each task's README.

## Recommended lightweight package

For reproduction, download [SR_HiSD_code_review_light.zip](SR_HiSD_code_review_light.zip)
and extract it. This self-contained lightweight package includes the
numerical code, required inputs, saved results, production logs and figures.
Run `main.m` inside the extracted `code_review_light/` folder. The commands
below apply from either package root: `code_review/` or `code_review_light/`.

## Reproduce all experiments

Open `main.m` in MATLAB and press **Run**, or run it from the package root:

```matlab
main
```

This runs all numerical experiments and automatically generates Table 1
and Figures 1–3. Each task saves its data and images in its own `results/`
folder. No separate plotting command is needed.

MATLAB R2025b with the JVM is supported. Onsager requires about 4 GB RAM;
full LB calculations require more memory and can take a long time.

## Experiments and outputs

All paths below are relative to the package root.

| Task | Article result | Main output |
|---|---|---|
| `morse2d` | Table 1: HiSD, SR-HiSD and NPSS | `morse2d/results/<RunId>/plot_table1/table1.csv`, `table1.png`, `table1.pdf` |
| `morse3d` | Figure 1: Morse landscape | `morse3d/results/<RunId>/plot_figure1/figure1.png`, `figure1.pdf` |
| `Onsager` | Figure 2: index 5 → 4 → 2 → 0 pathway | `Onsager/results/<RunId>/plot_figure2/figure2.png`, `figure2.pdf` |
| `lb` | Figure 3: BCC → HEX and BCC → disorder with critical nuclei | `lb/results/<RunId>/plot_figure3/figure3.png`, `figure3.pdf` |

One set of included figures and tables is stored in each task's
`results/paper/` folder. Required numerical inputs and saved states are in
`data/`; original production logs and timings are in `logs/`.

Each run gets a unique `RunId`. Logs are saved in
`<task>/logs/<RunId>_<subtask>.log`. `review_figures` lists the images
exported by `main`; `review_reports` contains each step's result paths.
Morse uses HiSD and SR-HiSD; LB uses LBSD and SR-HiSD.

## Run one task

Each command calculates the selected task and automatically plots it:

```matlab
run_task('morse2d');
run_task('morse3d');
run_task('Onsager');
run_task('lb');
```

## Run separate steps

Use the same options and `RunId` for related steps:

```matlab
options=struct('Mode','compute','RunId','onsager_run');
run_task('Onsager','search',options);
run_task('Onsager','plot_figure2',options);
```

| Task | Subtasks, in execution order |
|---|---|
| `morse2d` | `original_methods`, `npss`, `plot_table1` |
| `morse3d` | `compare_methods`, `landscape`, `plot_figure1` |
| `Onsager` | `search`, `plot_figure2` |
| `lb` | `bcc2hex`, `bcc2disorder`, `nuclei`, `plot_figure3` |

## Redraw saved results or check execution

```matlab
[reports,figures]=run_task('all','all',struct('Mode','reference'));
run_task('all','all',struct('Mode','check'));
run_task('all','all',struct('Mode','compute','Quick',true));
```

Reference mode redraws saved results. Check mode checks inputs. Quick mode
uses reduced budgets and does not produce the complete article results.
LB saved results include matching bilateral GF endpoints and all four
critical-nucleus measurements. Compute mode generates new GF endpoints.
Set `Visible=true` in the options to keep figure windows open.

Each task has separate `functions/`, `scripts/`, `data/`, `logs/` and
`results/` folders. Model settings and detailed outputs:
[Morse2D](morse2d/README.md), [Morse3D](morse3d/README.md),
[Onsager](Onsager/README.md), [LB](lb/README.md).
