# Sugarcane RNA-seq Metadata Scripts and Databases

This repository contains a collection of scripts and Jupyter notebooks designed to process, filter, and analyze metadata from sugarcane RNA-seq experiments. It includes tools for fetching missing DOIs/PMIDs, curating ontology, and generating visualizations.

## Repository Structure

```text
.
├── figures/            # Output visualization files (PNG, SVG)
├── raw/                # CSV datasets and metadata files
├── scripts/            # Python scripts and Jupyter notebooks
├── .gitignore          # Git ignore rules
└── README.md           # This file
```

## Data (`raw/`)

The `raw/` directory contains various stages of the sugarcane RNA-seq metadata:

- `databases_rna_sugarcane.csv`: Sugarcane database of RNA-seq metadata.
- `sugarcane_omics_metadata_cleaned.csv`: Cleaned version of the omics metadata.
- `sugarcane_rna_complete.csv`: Full RNA-seq dataset including initial metadata.
- `sugarcane_rna_filtered.csv`: RNA-seq metadata after preliminary filtering.
- `sugarcane_rna_relevant_filtered.csv`: Final filtered dataset containing only relevant columns.
- `sugarcane_with_pmid.csv`: Dataset updated with PubMed IDs.
- `sugarcane_with_pmid_with_articles_with_dois.csv`: Dataset containing PMIDs, article information, and DOIs.

## Scripts (`scripts/`)

### Python Scripts

These scripts are designed for automated data retrieval and database management:

- `get_all_info.py`: Fetches metadata (BioProject IDs, DOIs, PMIDs) and integrates it into the dataset.
- `create_json_db.py`: Processes CSV files to create or update a JSON-based metadata database.
- `get_doi_from_url.py`: Uses Selenium to scrape and extract DOIs from various web URLs.
- `get_scholar.py`: Scrapes Google Scholar to find missing DOIs and PMIDs for BioProjects.

### Jupyter Notebooks

These notebooks are used for interactive data processing and visualization:

- `create_figures_rna_sugar.ipynb`: Generates figures for publications, including plots for:
  - Country distribution
  - Growth conditions
  - Run counts per BioProject
  - Tissue types
  - Treatment types
- `filter_and_ontology.ipynb`: Handles manual filtering of metadata and curation of ontology terms (Genotype, Tissue, Treatment, Soil, Rainfall, etc.).
- `filter_data.ipynb`: Filters metadata and files based on the distribution of runs per BioProject.
- `filter_sugar.ipynb`: Performs the initial filtering of the sugarcane RNA-seq dataset.
- `load_rna.ipynb`: Loading and initial preparation of RNA-seq data.
- `load_rna_dois.ipynb`: Exploratory analys in the dataset with DOIs.

## Figures (`figures/`)

- `figure_panel.png`: A multi-panel visualization of the sugarcane metadata.
- `map.png`: Geographic distribution of the data sources.
- `overview.svg` / `overview.png`: Summary overview of the dataset characteristics.

## Setup and Requirements

The scripts in this repository are compatible with **Python 3.10+** (with some requiring **3.13+**).

### Using `uv` (Recommended)

This repository is managed with [`uv`](https://github.com/astral-sh/uv). It includes a `pyproject.toml` and a `uv.lock` file for absolute reproducibility of the environment.

To run a script using the project's environment:

```bash
uv run scripts/get_all_info.py [arguments]
```

To run Jupyter notebooks with the project's environment:

```bash
uv run jupyter notebook
```

### Manual Installation

If you prefer using standard `pip`, you will need to install the dependencies:

- `pandas`
- `selenium`
- `webdriver-manager`
- `rich`
- `regex`
- `...`
- `requests`

```bash
pip install pandas selenium webdriver-manager rich regex ... requests
```

> [!NOTE]
> We strong recommend to use `uv` to deal with all dependencies.

> [!NOTE]
> Note: Selenium scripts require a compatible web driver (e.g., ChromeDriver), which is managed automatically by `webdriver-manager` in most cases.

## Citation

If you use this repository or its datasets in your research, please cite it as:

```text
soon
```

Please check the `CITATION.cff` file for machine-readable citation information.

## License

This project is licensed under the MIT License.
