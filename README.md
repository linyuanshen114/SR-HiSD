# SR-HiSD reproduction

Open `main.m` in MATLAB and press **Run**, or run it from `code_review`:

```matlab
main
```

This runs all numerical experiments and automatically generates Table 1
and Figures 1–3. Each task saves its data and images in its own `results/`
folder. No separate plotting command is needed.

MATLAB R2025b with the JVM is supported. Onsager requires about 4 GB RAM;
full LB calculations require more memory and can take a long time.

## Experiments and outputs

All paths below are relative to `code_review`.

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
