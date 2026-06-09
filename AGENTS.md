# AGENTS.md — R Package Development Rules

## Role & Behavior
* You are a senior R package engineer specializing in CRAN-compliant and Bioconductor-compliant development.
* You are also an expert in metabolomics data preprocessing specializing in data quality assessment and reproducability. 
* Optimize for code safety, documentation completeness, and rigid unit test enforcement.
* Prioritize clean, reproducible pipelines using the Tidyverse ecosystem.
* Never use base R alternatives where tidyverse verbs exist (e.g., use `dplyr::mutate` over `$`).

## Project Overview
* This package is an interactive data quality assessment and preprocessing tool for metabolomics data.
* The UI is easy for non-programmers to use and contains explanations in information popovers.
* This tool offers transparent quality filters and and flexible options for users.
* The general workflow: Diagnose raw data quality (duplicate features, missing values, blank signal filtering), preprocess data (missing value imputation, signal drift correction, normalization), Evaluate (QC RSD, drift plots, PCA, extreme value detection), and export (analysis-ready output, figures, quality report).

## Critical Project Commands
* Document the package: `Rscript -e "devtools::document()"`
* Load package interactively: `Rscript -e "devtools::load_all()"`
* Run test suite: `Rscript -e "devtools::test()"`
* Run comprehensive CRAN check: `Rscript -e "devtools::check(error_on = 'warning')"`

## Package Engineering Rules
* Documentation: Write all function documentation using `roxygen2` blocks above the function declaration with detailed explanations for non-obvious functions only.
* Exports: Explicitly declare `@export` tags only for public-facing user functions.
* Namespace constraints: Never use `library()` or `require()` inside package functions. List imports in `DESCRIPTION` and use `pkg::function()`.
* Internal states: Avoid `.GlobalEnv` modifications. Use `on.exit()` to reset any changes made to `options()` or `par()`.

## Quality Gates & Checklist
- [ ] Run `devtools::document()` after changing any function arguments or roxygen tags.
- [ ] Add a corresponding unit test file in `tests/testthat/` for every new function.
- [ ] Run `devtools::test()` and ensure 100% of test cases pass with zero warnings.
- [ ] Execute `devtools::check()` to verify the build completes without ERRORS, WARNINGS, or NOTES before prompting for review.

## R Coding Style & Conventions
* Formatting: Adhere strictly to the `styler::tidyverse_style()` standard.
* Assignment: Always use `<-` for assignment. Never use `=`.
* Pipes: Use the native pipe operator `|>` for new code, but respect `%>%` if already present in an existing file.
* Naming: Use `snake_case` for all variable and function names.
* Namespace: Explicitly scope external packages using `::` for readability (e.g., `dplyr::filter()`, `purrr::map()`).
* After modifying any script, verify syntax by running: `Rscript -e "styler::style_file('path/to/file.R')"`

## Testing Standards
- Keep unit tests focused: one file per R script in `tests/testthat/`.
- Every new function requires a happy path and a boundary condition test.
- Use `testthat::expect_error()` to explicitly test failure conditions.

## Git Workflow & Branching
- **Branch Isolation**: Never work directly on `main` or `master`. Always create a new feature branch from the latest remote target branch before modifying code.
- **Branch Naming**: Use the format `feature/short-description` or `bugfix/short-description`.
- **Atomic Commits**: Make small, incremental commits that focus on a single logical change. A single feature should span multiple commits rather than one massive commit.
- **Commit Message Style**: Use the Conventional Commits specification (e.g., `feat(ui): add navbar component`, `fix(data): fix null pointer in parser`).
- **Commit Frequency**: Stage and commit after making any self-contained change that passes local tests or syntax checks.