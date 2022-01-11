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
* `data.csv` contains an expression table. Rows of the table should be grouped by genes, miRNAs, isomiRNAs and other items. Columns of the table are grouped by patients taken from two different groups.
* `description.csv` devide the patients on two non-intersecting groups (e.g. `Normal` and `Tumored` patients). It is assumed that a patient does not belong to the both groups simultaneously.
* (optionally) `interactions.csv` contains source/target pairs: in the network mode correlations will be computed only

### Regime of given interactions
As a comand line tool you can run th
```

```

### All vs all regime

