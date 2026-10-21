# pd-study-report

**Turn Proteome Discoverer result files into a clean HTML / PDF report you can send to collaborators.**

`pd-study-report` reads Thermo Proteome Discoverer (PD) `.pdResult`, `.msf` and `.pdStudy` files and writes a single, self-contained report covering:

- the study design
- the raw files
- every search, validation and quantification setting
- identification numbers
- software versions

It also drafts a methods paragraph you can start from.

- **No installation of packages:** pure Python standard library, a single file.
- **Read-only:** your PD files are opened in read-only mode and never changed.
- **Works on big files:** a 32 GB `.pdResult` with 2.4 million PSMs takes about 20 seconds.
- **Shareable:** one `.html` file (styles inlined, works offline) and an optional A4 `.pdf`.

---

## What the report contains

| Section | What you get |
|---|---|
| **Overview** | Key numbers at a glance (raw files, MS/MS spectra, PSMs, peptide groups, protein groups, high-confidence master proteins), the PD version, and who processed the data, when, and on which computer |
| **Study design** | Samples, quan channels and reporter-ion m/z, sample types, study factors, grouping and ratios, and the quantification method (e.g. TMT11) |
| **TMT channel loading** | Summed raw reporter abundance per channel with bars, so under- or over-loaded channels stand out |
| **Methods summary (draft)** | A methods paragraph generated from the stored parameters: search engine, databases, enzyme, tolerances, modifications, FDR thresholds and quantification settings |
| **Identification results** | Totals, high-confidence counts and counts passing the PD filters for PSMs, peptide groups, proteins, master proteins and protein groups |
| **Raw files** | Instrument, acquisition dates, file sizes, and per-fraction MS/MS spectra, filtered PSMs and ID rate with bars |
| **Key settings** | Database search, PSM validation, peptide/protein level and quantification settings in one table |
| **Workflows** | Processing and consensus workflows with template name, run time and node tree |
| **Processing log** | Warnings from PD plus useful log lines (database versions, Percolator results, total search time) |
| **Software versions** | The PD version that created the file, any later format upgrade, every node version, and the schema history |
| **Appendix** | Every visible parameter of every workflow node |

When the report contains a `.pdResult` together with its `.pdStudy`, the study design and settings are shown once, in the result section. If you only have a `.pdStudy`, the report shows the settings of its latest analysis.

---

## Requirements

- **Python 3.8 or newer.** Nothing to `pip install`.
- **For PDF output (optional):** Microsoft Edge, Google Chrome or Chromium. Edge comes with every Windows 10/11 PC. Without one of these, open the HTML report and use *Print → Save as PDF*.
- **Supported files:**
  - `.pdResult` and `.msf` from PD 2.x / 3.x (SQLite)
  - `.pdStudy` (XML)
  - other or older formats get a basic version listing

---

## Installation

```bash
git clone https://github.com/<your-user>/pd-study-report.git
cd pd-study-report
```

Or just download [`main.py`](main.py). It is the whole tool.

Check that Python works:

```bash
python --version          # Windows may also use: py --version
python main.py --version  # pd-study-report 1.0.0
```

---

## Usage

```bash
# 1. All PD files in the folder where main.py lives
python main.py

# 2. All PD files in another folder, with a PDF as well
python main.py "D:\Studies\MyStudy" --pdf

# 3. Specific files only
python main.py "D:\Studies\MyStudy\MyStudy.pdResult" "D:\Studies\MyStudy\MyStudy.pdStudy"

# 4. Hide user names, computer names and folder paths before sharing outside your group
python main.py "D:\Studies\MyStudy" --pdf --anonymize

# 5. Very large files: skip the slow peptide/PSM counts
python main.py "D:\Studies\MyStudy" --fast

# 6. Choose the output file and title, and open the report when done
python main.py "D:\Studies\MyStudy" -o reports/mystudy.html --title "Project X – TMT time course" --open
```

By default the report is saved next to the input files as `PD_report_<study name>.html` (and `.pdf`).

### Options

| Option | Meaning |
|---|---|
| `paths` | Files and/or folders to read (default: the folder of `main.py`) |
| `-o`, `--output` | Path of the HTML report |
| `--pdf` | Also write a PDF (A4, page numbers, appendix expanded) |
| `--anonymize` | Hide user names, computer names and folder paths |
| `--fast` | Skip counting peptide groups and PSMs (the slowest step on huge files) |
| `--title` | Report title (default: the study name) |
| `--open` | Open the report in your browser when it is ready |
| `--version` | Show the version |

You can also change the defaults at the top of `main.py` (the `settings` block), for example `ANONYMIZE = True` or `LIST_ALL_TABLES = True`.

---

## Example output

Console (file and study names shortened):

```text
> python main.py "D:\Studies\MyStudy" --pdf
Reading MyStudy.pdStudy (0.5 MB) ...
Reading MyStudy_TMT11.pdResult (32,133.3 MB) ...
  counting peptide groups ...
  counting PSMs per raw file (can take a minute for large files) ...

======================================================================
SUMMARY
  MyStudy_TMT11.pdResult: 2.2.0.388
  MyStudy.pdStudy: study file

HTML report: D:\Studies\MyStudy\PD_report_MyStudy.html
PDF report:  D:\Studies\MyStudy\PD_report_MyStudy.pdf
```

From the report, the *Identification results* section of a 24-fraction TMT11 experiment:

| Level | Total | High confidence | Passing PD filters |
|---|---:|---:|---:|
| MS/MS spectra in result | 1,842,848 | – | – |
| PSMs | 2,449,011 | 625,297 | 616,579 |
| Peptide groups | 1,342,076 | 255,965 | 255,965 |
| Master proteins | 10,302 | 9,284 | – |
| Protein groups | 10,709 | – | 10,302 |

The auto-generated *Methods summary (draft)* for the same file:

> 24 raw files acquired on Q Exactive HF Orbitrap were processed with Proteome Discoverer 2.2.0.388. Spectra were searched with Sequest HT against custom.fasta; Homo sapiens (SwissProt TaxID=9606_and_subtaxonomies) (v2017-06-07); MQ_contaminants.fasta, using Trypsin (Full) with up to 2 missed cleavages and a peptide length of 6–144 amino acids. Precursor and fragment mass tolerances were 10 ppm and 0.02 Da. Static modifications: TMT6plex / +229.163 Da (Any N-Terminus); Carbamidomethyl / +57.021 Da (C). Dynamic modifications: Oxidation / +15.995 Da (M); TMT6plex / +229.163 Da (K); TMTK8 / +237.177 Da (K); Label:13C(6)15N(4) / +10.008 Da (R); Acetyl / +42.011 Da (protein N-terminus). PSMs were validated with Percolator (strict / relaxed target FDR 0.01 / 0.05, based on q-Value). The strict target FDR was 0.01 at peptide level and 0.01 at protein level. Only peptides with at least high confidence were kept, and proteins needed at least 1 peptide sequence(s). Protein groups were formed with the strict parsimony principle. Reporter ions (TMT11_Beta) were quantified with an integration tolerance of 20 ppm (co-isolation threshold 50%, average reporter S/N threshold 5, peptides used: Unique + Razor). Abundances – normalization: Total Peptide Amount; scaling: On All Average.

---

## Before you share a report

- **Personal data:** without `--anonymize`, the report includes the PD user name, the processing computer's name and full file paths. Workflow, template, study and raw-file names are never changed, so rename them in PD if they contain anything private.
- **Methods draft:** the methods paragraph is generated mechanically. Check it against the settings table before using it in a manuscript.
- **Filtered numbers:** "Passing PD filters" counts reflect the filters that were stored in the result when it was processed.

---

## How it works

`.pdResult` and `.msf` files are SQLite databases. `.pdStudy` files are XML. The script opens them read-only and reads:

| Source in the file | Used for |
|---|---|
| `Workflows.WorkflowXML` | Workflow nodes and all of their parameters |
| `WorkflowInputFiles` | Raw files, instrument, acquisition dates |
| `AnalysisDefinition` and the `.pdStudy` `StudyDefinition` | Samples, channels, factors, quan method, file sets |
| `WorkflowMessages` | Run times, database versions, warnings |
| `SchemaInfo` | The PD version that created the file and later format upgrades |
| `TargetPsms`, `TargetPeptideGroups`, `TargetProteins`, `TargetProteinGroups`, `MSnSpectrumInfo` | Counts, per-file statistics and channel loading |

Opening a `.pdResult` in a newer PD version upgrades its file format. The report shows this under *Software versions* (for example "created with 2.2.0.388, format last upgraded by PD 3.3.1.110"). The processing itself was done with the version that created the file.

### Known limitations

- It has been tested with PD 2.2 results, including one later opened in PD 3.3. Other PD 2.x / 3.x versions should work: missing tables or columns are skipped rather than causing errors. Please open an issue if a section is empty for your files.
- PD 1.x `.msf` files are not parsed in detail; only version information is listed.
- The channel-loading table assumes that abundance arrays follow the sample order of the study. This is the case for PD's standard reporter-ion workflows.
- The PDF is printed by a headless Edge, Chrome or Chromium. If none is installed, print the HTML report from any browser.

---

## Contributing

Issues and pull requests are welcome, especially with examples of PD versions or workflows (label-free, DIA, cross-linking, ...) that are not covered yet. Please never attach real patient data or unpublished results to an issue. A screenshot of the empty section and the PD version are usually enough.

---

## Disclaimer

This is an independent community tool. It is not affiliated with, endorsed by or supported by Thermo Fisher Scientific. Proteome Discoverer is a trademark of Thermo Fisher Scientific.
