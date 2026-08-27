# NBTRP Web

Browser-only implementation of the published NBTRP XGBoost `survival:cox` risk score.

## Run locally

A web server is required because the model and example matrix are loaded with `fetch`:

```powershell
python -m http.server 8080 --directory NBTRP-web
```

Then open <http://localhost:8080>.

## Reproducibility

- Source model: `NBTRP-source/model/xgb.fit_model_new.rds`
- Converted with R `xgboost` 1.7.8.1 to `nbtrp-model.json`
- Browser predictor implements the 83 published regression trees and Cox objective transform.
- `test.js` compares all 223 browser predictions for the authors' E-MTAB-8248 example with R predictions.
- Maximum absolute difference in the current validation run: `< 2.3e-7`.

Run the comparison after regenerating R predictions:

```powershell
node NBTRP-web/test.js
```

## Important input behavior

The original R workflow scales each gene across the uploaded cohort (`scale(t(...))`). This website reproduces that cohort-wise Z-score, so at least two samples are needed. The authors' published example matrix contains 14 genes and omits two model features (`CLEC2D`, `UNG`); the original R indexing supplies those as missing values, and XGBoost uses the learned missing-value branches. The website reproduces that behavior but complete six-gene input is preferred.

The score is a relative Cox risk output, not a disease probability. Research use only.

Original NBTRP code/model copyright belongs to its authors and is distributed under the repository's MIT License; see `LICENSE`.
