# Meta-analysis: urbanization and biodiversity

This repository contains the R scripts use in the following study:

Endika Blanco-Urdillo, Alfredo Sánchez-Tójar, Juan D. Ibáñez-Álamo. Under review. *Methodological approaches to study the effect of urbanization on biodiversity: a review and meta-analysis*.

For any further information, please contact: [Alfredo Sánchez-Tójar](https://scholar.google.co.uk/citations?hl=en&user=Sh-Rjq8AAAAJ&view_op=list_works&sortby=pubdate), email: alfredo.tojar@gmail.com

## Scripts and data:

`data_preparation_urbanization.R`: script to prepare data for multilevel meta-analyses and meta-regressions on the relationship between urbanization and biodiversity. The code is used to exclude studies not providing enought information to be included in the analyses, estimate effect sizes, code variables, etc. `data_preparation_urbanization.R` needs the dataset `meta_AST_v2.csv` as input, and it creates the `meta_final_reduced.csv`as output, which will be then used by the following script.

`meta-analysis_urbanization.R`: script to run the actual multilevel meta-analyses and meta-regressions testing the relationship between urbanization and biodiversity, and exploring the importance of methodology in the study of such relationship and temporal, geographical and other patterns. The analyses are run using the R package '[metafor](http://www.metafor-project.org/doku.php/metafor)'. In addition, the script contains code for testing publication and time-lag bias, and for generating publication-ready plots for the multilevel meta-regressions run. `meta-analysis_urbanization.R` needs the dataset `meta_final_reduced.csv` as input. 

### Notes:

29th March 2019: The script available now will likely be split in two: processing data script and analytical script.

30th April 2019: The script was split in two, `data_preparation_urbanization.R` and `meta-analysis_urbanization.R`, to increase future re-usability.
