# Versions

All notable changes to **pd-study-report** are listed here. Version numbers follow [Semantic Versioning](https://semver.org/): `MAJOR.MINOR.PATCH`.

To check which version you have:

```bash
python main.py --version
```

To update, download the newest `main.py` (or run `git pull`) and replace the old file. There is nothing else to install.

---

## 1.0.0 – 2026-10-20

First public release.

### Added

- **HTML report:** one self-contained file that works offline, with light and dark screen themes and a print layout.
- **PDF report** (`--pdf`): printed with a headless Edge, Chrome or Chromium on A4 with page numbers. The appendix is expanded in the PDF.
- **Study design:** samples, quan channels, reporter-ion m/z, sample types, study factors, grouping, ratios and quan method, read from the `.pdResult` and from `.pdStudy` files.
- **TMT channel loading:** summed raw reporter abundance per channel over master proteins that pass the filters.
- **Methods summary (draft):** a methods paragraph generated from the stored workflow parameters.
- **Identification results:** totals, high-confidence counts and counts passing PD filters for PSMs, peptide groups, proteins, master proteins and protein groups, plus master proteins by protein FDR confidence.
- **Raw files:** instrument, acquisition dates and sizes, with per-fraction MS/MS spectra, filtered PSMs and ID rate.
- **Key settings:** database search, PSM validation, peptide/protein level and quantification settings in one table.
- **Workflows:** processing and consensus workflows with template, run time from the log, user, computer and node tree.
- **Processing log:** PD warnings and selected log lines (database versions, Percolator, search time).
- **Software versions:** the PD version that created the file, later format upgrades, node versions and the file schema history.
- **Appendix:** every visible parameter of every node.
- **Command-line options:** `paths`, `-o/--output`, `--pdf`, `--anonymize`, `--fast`, `--title`, `--open`, `--version`.
- **`--anonymize`:** hides user names, computer names and folder paths.

### Changed

- The output is now an HTML/PDF report instead of console text. The console shows only progress and a short summary.
- Results are listed before study files. When both describe the same study, the study design is shown once.

## 0.1.0 – initial script

- Printed the PD software version, `SchemaInfo`, `Workflows` and `WorkflowInputFiles` rows, version attributes found in the stored workflow XML, and unfiltered row counts to the console.

---

## Compatibility

| File type | PD version | Status |
|---|---|---|
| `.pdResult` | 2.2 (also after an upgrade by 3.3) | Tested |
| `.pdResult` / `.msf` | other 2.x and 3.x | Expected to work. Missing tables or columns are skipped. |
| `.pdStudy` | 2.x / 3.x (XML) | Tested with 2.2 |
| `.msf` | 1.x | Version information only |

| Python | Status |
|---|---|
| 3.12 | Tested |
| 3.8 – 3.11 | Expected to work (standard library only) |
