# DCoNA: tool for fast Differential Correlation Network Analysis
DCoNA is a statistical tool that allows one to identify pair interactions, which correlation significantly changes between two conditions. DCoNA was designed to test the hypothesis for a predefined list of source and target pairs ("Network" regime). However, DCoNA can also be used in the complete- network regime when the list is not given ("Exhaustive" regime). In this regime, DCoNA tests the hypothesis for all possible pairs of molecules from expression data.
Aside from the hypothesis testing, DCoNA can be used to test that significantly altered correlations of a particular source molecule are overrepresented among all significantly changed correlations. Also, DCoNA can compute mean, median, and other quantiles of z-statistics associated with a particular molecule and its targets to determine a trend in correlation changes.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents

- [Installation](#installation)
  - [Installation using pip](#installation-using-pip)
  - [Downloading example dataset](#downloading-example-dataset)
- [Usage](#usage)
  - [Example](#example)
  - [Available functions](#available-functions)
    - [`dcona.ztest`](#dconaztest)
    - [`dcona.zscore`](#dconazscore)
    - [`dcona.hypergeom`](#dconahypergeom)
  - [Data structure for CLI launch](#data-structure-for-cli-launch)
  - [Network and exhaustive regimes](#network-and-exhaustive-regimes)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## Installation

### Installation using pip
```
pip install dcona
```
### Downloading example dataset
You can try DCoNA on [TCGA-PRAD test dataset](https://eduhseru-my.sharepoint.com/:f:/g/personal/azhiyanov_hse_ru/Eo6INCepVSBDogyS5E9q-PkBDz_n_QDKUBf9kDcqHllpBw?e=kJdTdQ)



## Usage
You can use DCoNA either as Python-module or as a command-line tool.

### [Example](example/example.ipynb)
Detailed description of functions with data example and test launch.

### Available functions
#### `dcona.ztest`
**It tests the hypothesis on correlation equiavalence between pairs of genes**
``` python
dcona.ztest(data_df, description_df, reference_group, experimental_group, correlation='spearman', alternative='two-sided', interaction=None, repeats_number=None, output_dir=None, process_number=None)
```
* Command-line usage:
  ``` bash
  dcona ztest config.json
  ```

#### `dcona.zscore`
**It aggregates correlation changes of source molecule with all its targets.**  
``` python
dcona.zscore(data_df, description_df, reference_group, experimental_group, correlation='spearman', score='mean', alternative='two-sided', interaction=None, repeats_number=None, output_dir=None, process_number=None)
```
* Command-line usage:
  ``` bash
  dcona zscore config.json
  ```

#### `dcona.hypergeom`
**It groups pairs with changed correlations by the source molecules and finds overrepresented groups using the hypergeometric test.**  
``` python
dcona.hypergeom(ztest_df, alternative='two-sided', oriented=True, output_dir=None)
```
* Command-line usage:  
  You should launch `ztest` and then `hypergeom` with the same config file.
  ``` bash
  dcona hypergeom config.json
  ```

### Data structure for CLI launch
To run the tool in command line you need the following data:

* [`config.json`](example/configs/config.json) containing data filenames and tool usage parameters
```json
{
	"data_path": "./example/data/data.csv",
	"description_path": "./example/data/description.csv",
	"interaction_path": "./example/data/interactions.csv",
	"output_dir_path": "./../output/",
	
	"reference_group": "Normal",
	"experimental_group": "Tumor",

	"correlation": "spearman",
	"alternative": "two-sided",
	"score": "mean",
	"repeats_number": 500,
	"process_number": 2
}
```
Both relative and absolute file paths can be used.

Data description:

* `data_path` : `data.csv` contains an expression table. Rows of the table should be grouped by genes, miRNAs, isomiRNAs and other items. Columns of the table are grouped by patients taken from two different groups.

  Structure of `data.csv` :

  |        | sample_1 | ...  | sample_n |
  | :----: | :------: | :--: | :------: |
  | gene_1 |  1.2345  | ...  |  1.2345  |
  |  ...   |   ...    | ...  |   ...    |
  | gene_n |  1.2345  | ...  |  1.2345  |

  

* `description_path` : `description.csv` divide patients into **two** non-intersecting groups (e.g. `Normal` and `Tumor` patients). It is assumed that a patient does not belong to the both groups simultaneously.

  Structure of `description.csv`:

  |  Sample  |    Group    |
  | :------: | :---------: |
  | sample_1 | condition_1 |
  |   ...    |     ...     |
  | sample_n | condition_2 |

  Column names have to be exactly `Sample` and `Group`.

* `interaction_path` (*optional*): `interaction.csv` contains source/target pairs - correlations will be computed among this pairs (in `network` mode). You should delete this line from the config file if you want to launch an `exhaustive` mode.

  Structure of `interaction.csv`:

  |    Source     |    Target     |
  | :-----------: | :-----------: |
  | source_gene_1 | target_gene_2 |
  |      ...      |      ...      |
  | source_gene_n | target_gene_n |

  Column names have to be exactly `Source` and `Target`.

* `output_dir_path` is a path to an output directory.

Usage parameters:

* `reference_group`, `experimental_group` are names of the patient groups.

* `correlation` : `spearman` or `pearson`, defines the type of correlation that will be used in the tool.

* `alternative` : `two-sided`, `less` or `greater`. 

  TODO: describe the parameter meaning in `ztest` and `zscore` regimes.



### Network and exhaustive regimes

DCoNA has two working regimes:

* Network (interactions) regime - performs calculations only on given gene pairs. Requires an `interaction.csv` file.
* Exhaustive (all vs all) regime - generates all possible gene pairs from genes listed in `data.csv` and performs calculations. An `interaction.csv` file is not needed.
