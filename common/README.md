# Running the experiments

Run `main` from `code_review/` for the full calculations and automatic export
of Table 1 and Figures 1–3.

Use `run_task(task,subtask,options)` to select a task, step or calculation mode:

```matlab
run_task('morse3d');
run_task('Onsager','all',struct('Mode','reference'));
```

See the [reproduction guide](../README.md) for all selections.
The helpers select scripts, set paths, and save logs and run metadata.
Calculations are in each task's `scripts/`; mathematical functions are in
its `functions/`.

Results: `<task>/results/<RunId>/<subtask>/`.
Logs: `<task>/logs/<RunId>_<subtask>.log`.
Each step lists pictures in `figure_files` and `outputs`.
`main` collects their paths in `review_figures`; `run_task` returns the paths
as its second output.
