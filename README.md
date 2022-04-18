# DCoNA: tool for fast Differential Correlation Network Analysis
TODO: what is DCoNA and why should you use it?
TODO: put a short feature description.



<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li><a href="#installation">Installation</a></li>
      <ul>
          <li><a href="#installation-as-a-command-line-tool">Install1</a></li>
          <li><a href="#installation-from-source-code">Install2</a></li>
      </ul>
    <li><a href="#usage">Usage</a></li>
      <ul>
          <li><a href="#data-structure">Data structure</a></li>
          <li><a href="#Working-modes">Working modes</a></li>
            <ul>
              <li><a href="#ztest">Ztest</a></li>
              <li><a href="#hypergeom">Hypergeom</a></li>
              <li><a href="#zscore">Zscore</a></li>
            </ul>
          <li><a href="#Network-and-exhaustive-regimes">Network and exhaustive regimes</a></li>
      </ul>
  </ol>
</details>








## Installation

### Installation as a command line tool
TODO: installation from pip3

### Installation from source code
Start with cloning the repo. It is really crucial to clone it with all the submodules as follows:
```
git clone --recurse-submodules git@github.com:zhiyanov/DCoNA.git
```
#### Linux (Ubuntu)
TODO: change the example below
First, install the necessary dependencies:
```
sudo apt-get install -y g++ make cmake python3-dev python3-pip python3-numpy
sudo pip3 install cython pot
```
#### OS X
All instructions are the same, except that you need to install dependencies differently, through a combination of `pip3` and `brew`.

### Downloading TCGA-PRAD dataset
TODO: put the link here.

## Usage

### Data structure
To run the tool you need the following data
* `config.json` data filenames and tool usage parameters:
```json
{
	"data_path": "~/input_directory/data.csv",
	"description_path": "~/input_directory/description.csv",
	"interaction_path": "~/input_directory/interaction.csv",
	"output_dir_path": "~/output_directory",
	
	"reference_group": "Normal",
	"experimental_group": "Tumor",

	"correlation": "spearman",
	"alternative": "two-sided",
	"score": "mean",
	"repeats_number": 800,
	"process_number": 64,

	"fdr_threshold": 0.05
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

* `interaction_path` : `interaction.csv` (optionally) contains source/target pairs: correlations will be computed among this pairs (in `network` mode). You should delete this line from the config file in `exhaustive` mode.

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

### Working modes

#### Ztest

TODO: description

Usage:

```bash
# For network regime
python3 ~/network/ztest.py config.json
# For exhaustive regime
python3 ~/exhaustive/ztest.py config.json
```

#### Hypergeom

TODO: description

Usage:

```bash
# For network regime
python3 ~/network/hypergeom.py config.json
# For exhaustive regime
python3 ~/exhaustive/hypergeom.py config.json
```

#### Zscore

TODO: description

Usage:

```bash
# For network regime
python3 ~/network/zscore.py config.json
# For exhaustive regime
python3 ~/exhaustive/zscore.py config.json
```

### Network and exhaustive regimes

DCoNA has two working regimes:

* Network (interactions) regime - performs calculations only on given gene pairs. Requires an `interaction.csv` file.
* Exhaustive (all vs all) regime - generates all possible gene pairs from genes listed in`data.csv` and performs calculations. An `interaction.csv` file is not needed.
