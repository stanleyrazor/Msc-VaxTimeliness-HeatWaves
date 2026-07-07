# lshtm_msc_project_template

The goal of lshtm_msc_project_template is to provide a template [Quarto](https://quarto.org) book for a MSc project at [LSHTM](https://www.lshtm.ac.uk)

It provides a title page, with all the required components configurable through `_quarto.yml`, and a framework for writing up various sections of a project.

As a demonstration the quarto templates here have been rendered to a pdf file in [_book/[Candidate-Number]_[MSc]_[Year-of-Submission]_Project.pdf](https://github.com/bakenzua/lshtm_msc_project_template/raw/master/_book/%5BCandidate-Number%5D_%5BMSc%5D_%5BYear-of-Submission%5D_Project.pdf)

## Install

Either clone this git repository:

```         
git clone https://github.com/bakenzua/lshtm_msc_project_template.git
```

or download the [zip file](https://github.com/bakenzua/lshtm_msc_project_template/archive/refs/heads/master.zip) and extract somewhere appropriate on your pc.

## Configure

Configuration of the title page can be made by editing the metadata in `_quarto.yml`

Editable options are indicated via comments.

``` yaml
################################################
#   Title page metadata - edit accordingly     #
################################################
candidate-number: "007"
msc-course: "Medical Statistics"
academic-year: "2023-24"
page-count: "0"
project-length: "Standard" # Extended (MSc IID only)
################################################
```

#### Output file

Extra care should be taken with the `output-file` option as this does not compile from the other options:

``` yaml
# edit line below...
  output-file: "[Candidate Number]_[MSc]_[Year of Submission]_Project"
```

Replacing

-   `[Candidate Number]` with the candidate number, e.g. 1234

-   `[MSc]` with MSc title, e.g. MedSats

-   `[Year of Submission]` with year of submission, e.g. 2023

Will produce an output file: `1234_MedStats_2023_Project.pdf` as required

## Use

The project includes the following quarto files that are applicable for different sections of a project.

-   index.qmd

    -   Title page

    -   Table of contents

    -   Abstract

    -   Acknowledgements

-   intro.qmd

    -   Introduction

    -   Aims and Objectives

-   methods.qmd

-   results.qmd

-   discussion.qmd

-   references.qmd

-   appendices.qmd

Please refer to the [Quarto documentation](https://quarto.org/docs/guide/) for guidance on incorporating workings and write up into these files.

Note that each of the files is rendered within it's own R environment, meaning data frames and functions are not passed from e.g. `methods.qmd` to `results.qmd`. If there is a need for common objects between files, some custom caching or common set-up code must be implemented. One method might be to implement a `{targets}` pipeline with this book as a quarto target, which loads parent objects as required. Alternatively concatenate different chapters into one file, i.e. write up methods and results section in the same `qmd` file

The easiest way to compile the template to a pdf, is via the render button in RStudio.

## Citations and References

The `references.bib` file is a bibliography file as [described in Quarto documentation](https://quarto.org/docs/authoring/citations.html#bibliography-files). This can be exported from a reference manager, e.g. [Zotero](https://www.zotero.org/). To cite an article from the bibliography in the text, use \@ followed by the citation identifier, as per the standard pandoc [citation syntax](https://quarto.org/docs/authoring/citations.html#sec-citations). An example is given in `intro.qmd`. Cited articles will then automatically be compiled to the references page.

The `vancouver.csl` file is a [citation style](https://citationstyles.org/) definition which uses the Vancouver referencing and citation system. Alternative citation and reference styles can be downloaded from <https://github.com/citation-style-language/styles>. The `csl: vancouver.csl` field of `_quarto.yml` should be amended accordingly.

## Attribution

-   This template was inspired by the [latex template for LSHTM MSc projects](https://www.overleaf.com/latex/templates/lshtm-msc-project-template/ddvdhnnrcwbz) by Robert Greener.
-   It borrows heavily from [this](https://cameronpatrick.com/post/2023/07/quarto-thesis-formatting) post by Patrick Cameron.
-   The LSHTM Logo (`lshtm_logo.jpg` and `lshtm_logo_small.jpg`) remains property of the London School of Hygiene and Tropical Medicine
    -   Everything else is MIT licensed

## Thanks

I would love to know if you find this project useful, please star the repository or [buy me a coffee](https://buymeacoffee.com/bakenzua).

If you have any issues or find any faults with this repo, please (in order of preference 😃) 

-   Fix it and submit a pull request
-   Raise an issue above
-   Contact me directly
