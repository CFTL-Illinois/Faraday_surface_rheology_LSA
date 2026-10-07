# LSA_Faraday_waves

MATLAB code for the Floquet linear-stability calculations and figure reproduction.

## Folder structure

- `src/` — solver routines
- `scripts/Final_data.m` — recomputes all data and writes `data/JFM_final_data.mat`
- `scripts/Final_plot.m` — reproduces the manuscript figures from `data/JFM_final_data.mat`
- `data/` — saved numerical data
- `figures/` — generated figures

## Quick start

Open MATLAB in the repository root and run:

```matlab
run('scripts/Final_data.m')
run('scripts/Final_plot.m')
```

The scripts add `src/` to the MATLAB path automatically.

If `data/JFM_final_data.mat` is already present, figures can be regenerated without recomputing:

```matlab
run('scripts/Final_plot.m')
```

Generated figures are written to `figures/`.
