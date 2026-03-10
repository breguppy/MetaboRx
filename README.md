# QCcorrection
Shiny app for vetting data quality and correction signal drift in metabolomics data

### Requirements

 - R >= 4.4.1
 
## Installation

#### Install Via GitHub
To install via GitHub, use the remotes package:
```r
install.packages("remotes")

```
You will also need to install impute via BiocManager
```r
install.packages("BiocManager")
BiocManager::install("impute")
```

To install the QCcorrection package:
```r
# Install QCcorrection (core functionality)
remotes::install_github("breguppy/QCcorrection", dependencies = TRUE)

# make sure required packages are installed
QCcorrection::check_required_dependencies()

# install optional dependencies
QCcorrection::install_optional_dependencies()
```

## To run App
```r
QCcorrection::run_app()
```
## Data Input and Information Requirements

### 1. Acceptable file formats

- `.csv`
- `.xls`
- `.xlsx`

**Note:**  
Raw data must be located on the **first sheet** of `.xls` or `.xlsx` files.

---

### 2. Required data formatting

- **Rows = samples** (can be in any order)
- **Columns = non-metabolite columns and metabolites** (can be in any order)

#### Non-metabolite columns

- **Sample column (required):**  
  Column containing unique sample identifiers.

- **Batch column (optional):**  
  Column containing batch information if samples were acquired in batches.

- **Class column (required):**  
  Column indicating sample type.  
  Must contain QC samples labeled as `NA`, `QC`, `Qc`, or `qc`.  
  If blank samples are present, they must be labeled as `blank`.

- **Injection order column (required):**  
  Column indicating the order in which samples were injected.

- **Additional metadata columns (optional):**  
  Any remaining non-metabolite columns must be explicitly specified.

---

### 3. Injection order requirements

Data (excluding blank samples) must **begin and end with QC samples** when sorted by injection order.

---

### 4. Internal standard metabolites

Internal standard metabolites must have column names beginning with:

- `ISTD`
- `ITSD`

### Example Raw Data Stucture
<img align="center" src="https://github.com/breguppy/QCcorrection/blob/main/www/example_data_structure.png">

# Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/breguppy/QCcorrection/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/breguppy/QCcorrection/issues).

#### [Pull Requests](https://github.com/breguppy/QCcorrection/pulls) are welcome for bug fixes, new features, or enhancements.

## Citation

Coming Soon
