# EVSC 410 JupyterLite teaching model — v3

Browser-oriented teaching prototype; **not validated for restoration design**. Width is spatially uniform in each run. Requires only NumPy and Matplotlib at runtime.

## Deploy

1. Create a GitHub repository and upload this directory (including `.github/workflows/deploy.yml`).
2. In repository **Settings → Pages**, select **GitHub Actions** as the source.
3. Push to `main`, then wait for the workflow to build and deploy JupyterLite.
4. Open the published Pages URL, then `EVSC410_Experiment_1.ipynb`.

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
