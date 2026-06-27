# GC32 Gaia-HST-RV/PM Dynamical-State Index Reproducibility Package

This archive contains the source inputs, derived tables, figures, code, logs, and metadata for the manuscript:

**A Consensus Gaia-HST-Radial-Velocity Dynamical-State Index for Galactic Globular Clusters**

The package is structured for a GitHub release and a Zenodo DOI deposit.

## Directory structure

- `data/input` - source CSV tables used by the MATLAB pipeline.
- `outputs/tables_csv` - machine-readable result tables used in the manuscript and supplementary material.
- `outputs/supplementary_tables.xlsx` - human-readable workbook version of the supplementary tables.
- `outputs/logs` - execution log and analysis summary from the archived run.
- `outputs/workspace` - archived MATLAB workspace from the run.
- `figures` - manuscript and supplementary figures in PDF and PNG formats.
- `code/matlab` - MATLAB analysis and figure-generation scripts.
- `code/python` - package checksum validator.
- `manuscript` - manuscript PDF and supplementary material.
- `editorial` - journal-upload text files.
- `metadata` - file inventory, checksums, source mapping, and package metadata.
- `environment` - runtime notes.

## Reproduce the numerical analysis

Open MATLAB in the package root and run:

```matlab
addpath(genpath('code/matlab'));
run('code/matlab/run_gc32_dsi_analysis_pipeline.m');
```

The recomputed outputs are written to `outputs/recomputed`.

## Regenerate figures from archived tables

```matlab
addpath(genpath('code/matlab'));
make_gc32_dsi_manuscript_figures('outputs/tables_csv', ...
    'data/input/13_gc32_structural_validation.csv', ...
    'outputs/recomputed/figures', false);
```

## Validate package integrity

```bash
python code/python/validate_package.py
```

The validator checks every file listed in `metadata/sha256_manifest.csv`.

## Main data blocks

The strict sample contains 32 Galactic globular clusters with Gaia EDR3, HST proper-motion, and RV/PM velocity-dispersion data blocks. The archived input layer contains 389 HST proper-motion profile rows, 3232 Gaia rotation/dispersion rows, and 839 RV/PM velocity-profile rows. APOGEE chemistry and structural parameters are included as auxiliary validation layers and are not used to construct the DSI.

## License notes

- Code in `code` is provided under the MIT License.
- Derived tables, figures, and documentation are provided under CC BY 4.0 unless a source-data provider requires additional attribution or reuse terms.
- Public survey and catalogue products retain the citation and reuse requirements of their original providers.

## DOI and GitHub release

Before public release, update these fields after creating the GitHub repository and Zenodo record:

- `CITATION.cff`: `repository-code`, `url`, and DOI fields.
- `.zenodo.json`: DOI/version fields if assigned manually.
- `editorial/data_and_code_availability_statement.txt`: Zenodo DOI and GitHub release URL.
