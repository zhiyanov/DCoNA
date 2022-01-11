# DCoNA: tool for fast Differential Correlation Network Analysis
TODO: what is DCoNA and why shold you use it?
TODO: put a short feature description.

## Installation instructions

### Installation as a comand line tool
TODO: installation from pip3

### Installation from source code
Start with cloning the repo. It is really crucial to clone it with all the submodules as follows:
```
git clone --recurse-submodules git@github.com:zhiyanov/DCoNA.git
```
#### Linux (Ubunt)
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

## Usage examples

### Data structure
To run the tool you need the following data
* `config.json` data filenames and tool usage parameters:
```
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

	"fdr_treshold": 0.05
}
```
Data description:
* `data_path` : `data.csv` contains an expression table. Rows of the table should be grouped by genes, miRNAs, isomiRNAs and other items. Columns of the table are grouped by patients taken from two different groups.
* `description_path` : `description.csv` devide the patients on two non-intersecting groups (e.g. `Normal` and `Tumor` (`Tumor`) patients). It is assumed that a patient does not belong to the both groups simultaneously.
* `interaction_path` : `interaction.csv` (optionally) contains source/target pairs: correlations will be computed among this pairs (in `network` mode). You should to delete this line from the config file in `exhaustive` mode.
* `output_dir_path` is a path to an output directory.
Usage parameters:
* `reference_group`, `experimental_group` are names of the patien groups.
* `correlation` : `spearman` or `pearson`, defines the type of correlation that will be used in the tool.
* `alternative` : `two-sidede`, `less` or `greater`. TODO: describe the parameter meaning in `ztest` and `zscore` regimes.

### Regime of given interactions
As a comand line tool you can run th
```

```

### All vs all regime

