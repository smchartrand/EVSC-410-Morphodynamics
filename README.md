# EVSC 410 JupyterLite teaching model — v8 (progress + optional output exports)

Browser-oriented teaching prototype; **not validated for restoration design**. Width is spatially uniform in each run. Requires only NumPy and Matplotlib at runtime.

## Deploy

1. Create a GitHub repository and upload this directory (including `.github/workflows/deploy.yml`).
2. In repository **Settings → Pages**, select **GitHub Actions** as the source.
3. Push to `main`, then wait for the workflow to build and deploy JupyterLite.
4. Open the published Pages URL, then `EVSC410_Experiment_1.ipynb` (now includes model equations and the discretized Exner balance).

The JupyterLite/Pyodide kernel runs in the student's browser; no Colab or server-side Python. This build has not been verified in a real browser session; test on a low-powered machine before classroom use. First load may require downloading Python packages and caching them.

## Numerical safeguards and efficiency

- Analytical sediment transport slope sensitivity replaces repeated finite-difference evaluations.
- Explicit conservative finite-volume Exner with adaptive time step and positive-slope rejection.
- Preflight displays hydraulic Courant number (informational) and diffusive morphodynamic stability number (time-step limiting).
- Store only requested snapshots; not every integration step.
- Classroom limits: 500 cells, 60 stored profiles, 150,000 steps. Initial work-size warning above 250 cells or 30,000 estimated steps.
- No spatial width variation, backwater, sorting, or varying discharge.

## Testing

`python tests/test_model.py`

## Caution

Stability is not convergence. The model uses a wide-channel normal-flow closure and cannot simulate zero/adverse slopes. The downstream boundary extrapolates the final interior bed slope. The feed boundary is imposed, and bed evolution is conservative for the stated scheme. Refinement and external benchmark validation remain necessary.

## Teaching background

The notebook now includes normal-depth derivation, Shields stress, the single-grain Meyer-Peter–Müller-type relation implemented in the code, Exner sediment continuity, and the actual finite-volume discretization. The core numerical method is unchanged; an asynchronous wrapper yields to the notebook interface to update progress.

### Citations
The notebook now contains author–year citations and a references section. The numerical scheme is unchanged; a progress display is added to the simulation cell.

## Simulation progress
The experiment cell displays one progress bar across the three sediment-feed scenarios. It uses IPython HTML display (no ipywidgets dependency) and `await simulate_with_progress(...)` to yield periodically so Pyodide can repaint the interface. The original synchronous `simulate(...)` API remains available.


## Optional model outputs (v8)

Set `SAVE_OUTPUT = True` in the notebook setup cell and run all cells. The notebook writes `EVSC410_Output/scenario_0.7`, `scenario_1.0`, and `scenario_1.3`, each with:

- `model_results.csv`: bed elevation, water depth, mean velocity, bedload transport per unit width, Shields stress, and bed slope at all `save_count` output times and cell centres. Face quantities are averaged to cell centres for CSV reporting; the underlying result arrays retain their original face values.
- `bed_elevation.png`: bed profiles at the requested output times.
- `sediment_transport.png`: bedload transport profiles at the requested output times.

`EVSC410_Output.zip` is also created. In JupyterLite, right-click the ZIP in the file browser and choose **Download**. The virtual filesystem is browser-local; saving there does not write to the student's local downloads folder. Re-running replaces the same-named files. No D50 or D84 histories are generated because this model uses one fixed grain size.

Both `content/teaching_model.py` and `content/EVSC410_Experiment_1.ipynb` must be updated together.
