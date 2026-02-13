
## Docker environment

This project uses Docker to provide a fully reproducible environment for the exam pipeline (scRNA-seq + scATAC-seq, Steps 3–14). The analysis relies on packages with non-trivial system dependencies (e.g., HDF5, Bioconductor genomics stack, Seurat/Signac) and on document generation via R Markdown. Running everything in a container ensures that the same versions of system libraries and R packages are available across machines, avoiding “it works on my computer” issues.

## Data availability

Input datasets are not included in this repository.  
To reproduce the analysis, place the required input files under `/workspace/data`, export the data and run the report.

## Why Docker?

Docker is recommended here because it:

- **Guarantees reproducibility**: the same OS, system libraries, and R/Python packages are used every time.
- **Avoids dependency conflicts**: genomics packages often require specific system libraries (HDF5, SSL, XML, Cairo, compression libraries).
- **Simplifies setup**: no manual installation of dozens of packages on the host.
- **Keeps the host clean**: the container isolates the project runtime from the local R/Python installations.
- **Supports the full workflow**: interactive development in **JupyterLab** and report generation with **R Markdown**.

## What the Dockerfile provides

The Docker image is built from `rocker/r2u:24.04` (Ubuntu 24.04 + R with r2u), and installs:

- System tools and libraries required by single-cell workflows (HDF5, compression libs, samtools/tabix, etc.)
- CRAN packages as Ubuntu binaries via r2u (e.g., **Seurat**, **Signac**, **Matrix**, **rmarkdown**, **devtools**, …)
- Bioconductor binaries (e.g., **DropletUtils**, **SingleCellExperiment**, **GenomicRanges**, …)
- A Python virtual environment with **JupyterLab**
- The **R kernel** for Jupyter (`IRkernel`)

To reduce memory spikes from BLAS/OpenMP parallelism, the image sets:
`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, etc.



## Reproducibility with Docker

To ensure a consistent environment across machines, the full analysis was executed inside a Docker container. The image includes all required system dependencies and R/Python packages, so the workflow can be reproduced without manual setup on the host system.

From the repository root (where the `Dockerfile` is located), build the image:

```bash
docker build -t exam_r:latest .
````

Run the container, exposing JupyterLab on port `8888` and mounting the current directory into `/workspace` so that inputs/outputs persist on the host:

```bash
docker run -it -p 8888:8888 -v "$PWD":/workspace --name exam_cont exam_r:latest
```

After starting the container, open JupyterLab in a browser using the URL shown in the terminal logs (typically `http://localhost:8888` with an access token).

To restart the same container in later sessions:

```bash
docker start exam_cont
```
To restart the container:
```bash
docker stop exam_cont
```


## Reproducibility and project organization

This project was developed to ensure full reproducibility of the analysis. The original exploratory code (scripts/notebooks) was refactored into a structured R package so that each analysis step can be re-run deterministically, with clear inputs/outputs and documented functions.

Concretely, the development workflow was:

- start from generic analysis code (scripts/notebooks),
- refactor the workflow into modular R functions (one functions per pipeline step),
- add function documentation (help pages) and export rules,
- define package metadata (`DESCRIPTION`) and exports/imports (`NAMESPACE`),
- create a vignette (`vignettes/`) describing the pipeline and rationale,
- create a complete R Markdown report (`Report.Rmd`) that runs the full analysis end-to-end and renders an HTML output with results and plots.

## Install and use the R package

All pipeline steps are implemented as an R package (`examPipeline`) with documented functions, a vignette, and a full R Markdown report.

### 1) Install the package from source (inside the container)

The package source is located in:

* `/workspace/examPipeline` (contains `DESCRIPTION`, `NAMESPACE`, `R/`, `man/`, `vignettes/`)

#### Option A — Install functions only (fast install)

```bash
R -q -e "install.packages('/workspace/examPipeline', repos = NULL, type = 'source')"
```

#### Option B — Install including vignette compilation (recommended)

This compiles the `.Rmd` file in `vignettes/` and makes it available via `browseVignettes()`.

```bash
R -q -e "devtools::install('/workspace/examPipeline', upgrade = 'never', build_vignettes = TRUE)"
```

---

### 2) Load the package

```r
library(examPipeline)
```

---

### 3) Access function documentation (help pages)

Open the help page for a specific function:

```r
?step3_select_scrna_h5
```

or:

```r
help("step3_select_scrna_h5")
```

List all exported functions:

```r
ls("package:examPipeline")
```

Open the package documentation index:

```r
help(package = "examPipeline")
```

---

### 4) Access the vignetteb - If accessing the compiled HTML manually (e.g., in Jupyter)

```r
library(examPipeline)

html <- system.file("doc", "pipeline_steps_3_14.html", package = "examPipeline")

IRdisplay::display_html(
  paste(readLines(html, warn = FALSE), collapse = "\n")
)
```

Note: the vignette HTML is available only if the package was installed with `build_vignettes = TRUE` or if it was previously built and included in `inst/doc`.


### 4) Render the full analysis report (R Markdown)

The file `Report.Rmd` is the end-to-end execution document of the pipeline.  
It runs all analysis steps in the correct order and reports the corresponding outputs, including:

- **Inputs and paths** used in the Docker environment (`/workspace/data`, `/workspace/refs`, `/workspace/results`).
- **Step-by-step execution** of the pipeline (Steps 3–14), calling the package functions.
- **Intermediate artifacts** created at each step (e.g., exported MTX bundles, serialized `.rds` objects, reference files).
- **Sanity checks and summaries**, such as:
  - matrix/object dimensions (genes × cells / features × cells),
  - number of retained cells per sample,
  - lists of files produced in each output folder.
- **Plots and figures** generated by the workflow (e.g., UMAP visualizations and other step-specific plots), embedded in the HTML output.
- **Short result comments** after each step, describing what the output represents and how it is used downstream.

The report is rendered to HTML to provide a single, readable deliverable that combines code execution, outputs, and interpretation in one place.

To generate it inside the container:

```bash
R -q -e 'rmarkdown::render("Report.Rmd", output_format="html_document")'

````

This generates an HTML report (by default in the same directory as `Report.Rmd` unless an output directory is specified).


## Online resources

- Source code (GitHub):  
  https://github.com/francescacoriasco/exam

- Docker Hub repository:  
  https://hub.docker.com/repositories/francescacoriasco
