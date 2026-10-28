# GitHub and journal release

## Status of this snapshot

The repository structure, Python scripts, supplied inputs, regenerated figures and numerical reports are prepared for upload. The checked-in workflow is configured for GitHub; no claim is made that a remote workflow has already run. The public repository URL, archival DOI and a publication DOI are intentionally absent until they exist.

Before describing the scientific data package as final, resolve the first Aquistore survey date documented in DATA_PROVENANCE.md, reconcile the precise Sleipner t0=1998.48 with the manuscript table, and confirm the intended reuse licence with the rightsholders. All current numeric inputs are preserved and all choices are recorded.

## Create a repository

Unzip the package, open a terminal in its root, and run:

```bash
git init -b main
git add .
git commit -m "Add reproducible CO2 plume figures and verification"
```

The commit uses your configured Git identity. If GitHub CLI is installed and authenticated, this command creates a private review repository in your account:

```bash
gh repo create co2-plume-reproducibility --private --source . --push
```

Alternatively, create an empty GitHub repository through its website and use the push instructions supplied by GitHub. Choose public visibility when the authors are ready to release the work.

## Reproducibility check

```bash
python -m pip install -r requirements.txt
python scripts/reproduce.py --recompute --validate
```

Inspect the generated plots and the reports in results/. The automated tests assess numerical/software reproducibility; they do not resolve the documented field-source discrepancy.

## Cite and archive

CITATION.cff contains manuscript-author names, software title and version. Add the actual `repository-code` URL after creating the remote. Create a versioned GitHub release and archive that release in a DOI-granting repository when appropriate. Add the actual DOI to CITATION.cff after it is assigned, and cite the archived version in the manuscript's Data Availability statement.

A suitable data-availability description once the repository is accessible is:

> The Python scripts, numerical history, processed seismic-radius table, and verification results required to reproduce Figures 2–4 are available in the archived software release cited here. The field-radius analysis starts from processed image-derived measurements; their source publications and the numerical methods are documented in the repository.

Use the real release citation in the paper. This snapshot contains no fabricated DOI or repository link.
