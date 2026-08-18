# Sugarcane RNA-seq Metadata Scripts and Databases

This repository contains a collection of scripts and Jupyter notebooks designed to process, filter, and analyze metadata from sugarcane RNA-seq experiments. It includes tools for fetching missing DOIs/PMIDs, curating ontology, and generating visualizations.

## Repository Structure

```text
.
├── figures/            # Output visualization files (PNG, SVG, Gephi)
├── raw/                # CSV datasets, metadata files, and network graphs
├── scripts/            # Python scripts, Jupyter notebooks, and Bash scripts
├── .gitignore          # Git ignore rules
├── CITATION.cff        # Citation information in CFF format
├── LICENSE             # MIT License file
├── README.md           # This file
├── datapackage.json    # Data package specification
├── pyproject.toml      # Python project configuration
└── uv.lock             # Lockfile for environment dependencies
```

## Data (`raw/`)

The `raw/` directory contains various stages of the sugarcane RNA-seq metadata:

- `databases_rna_sugarcane.csv`: Sugarcane database of RNA-seq metadata.
- `figure2_ontology_terms.graphml`: Network graph of ontology terms in GraphML format.
- `figure2_ontology_terms.pdf` / `figure2_ontology_terms.png`: Static visualizations of the ontology terms network (Figure 2).
- `figure2_supplementary_interactive.html`: Interactive web visualization of the ontology terms network.
- `graph_sugarcane.graphml`: Network graph of sugarcane metadata connections in GraphML format.
- `libtype.csv`: Dataset containing library type and strandedness information from Salmon outputs.
- `ontology_table.csv`: Reference mapping table of curated ontology terms.
- `SraRunTable.csv`: Metadata table downloaded from NCBI SRA containing run attributes.
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
- `create_graph_fig.ipynb`: Generates network graph figures and interactive visualizations for ontology metadata.
- `filter_and_ontology.ipynb`: Handles manual filtering of metadata and curation of ontology terms (Genotype, Tissue, Treatment, Soil, Rainfall, etc.).
- `filter_data.ipynb`: Filters metadata and files based on the distribution of runs per BioProject.
- `filter_sugar.ipynb`: Performs the initial filtering of the sugarcane RNA-seq dataset.
- `load_rna.ipynb`: Loading and initial preparation of RNA-seq data.
- `load_rna_dois.ipynb`: Exploratory analys in the dataset with DOIs.
- `r1.ipynb`: Revision analysis notebook for generating updated plots, library type analysis, and completeness metrics.

### Bash scripts

These scripts are used to run the `nf-core` pipelines and compress the files.

- `pipeline_rnaseq_salmon.sh`: run the `nf-core` `fetchngs` and `rnaseq` using a `.csv` file with SRA IDs.
- `compress_runs.sh`: extract the information from the complete results and save into a `.tar.xz` file.
- `compress_runs_parallel.sh`: Multi-threaded parallel version of `compress_runs.sh` for archiving results.
- `extract_libtype.sh`: get the libType from Salmon output.

### Gephi Project Files

These project files are used for interactive graph visual editing:

- `graph_r1.gephi`: Gephi workspace file containing network graph layout for revision analysis.

## Figures (`figures/`)

- `barplot.png`: Bar plot showing metadata category distributions.
- `completeness_plot.png` / `figure_completeness.png` / `figure1_completeness.png`: Visualizations of dataset completeness across metadata attributes.
- `figure_panel.png`: A multi-panel visualization of the sugarcane metadata.
- `graph.gephi` / `graph.png` / `graph.svg` / `graph2.png` / `graph2.svg`: Network graph representations showing relationships among metadata terms.
- `graph_r1.png` / `graph_r1.svg`: Revision network graph visualizations.
- `growth_conditions.png`: Breakdown plot of plant growth condition categories.
- `interface.png`: Screenshot overview of the interactive visualization interface.
- `map.png`: Geographic distribution of the data sources.
- `overview.svg` / `overview.png`: Summary overview of the dataset characteristics.
- `run_counts.png`: Distribution plot of sequencing run counts per BioProject.
- `tissue_distribution.png`: Breakdown plot of tissue types represented in the dataset.
- `treatment_distribution.png`: Breakdown plot of experimental treatments represented in the dataset.

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
> Selenium scripts require a compatible web driver (e.g., ChromeDriver), which is managed automatically by `webdriver-manager` in most cases.

## Citation

If you use this repository or its datasets in your research, please cite it as:

```text
soon
```

Please check the `CITATION.cff` file for machine-readable citation information.

## License

This project is licensed under the MIT License.
