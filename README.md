SBROME
======

Synthetic Biology Reusable Optimization Methodology

# PAMDB
A Parts Database with Consensus Parameter Estimation for Synthetic Circuit Design

### What is PAMDB?
PAMDB aims to bridge the gap between published circuits and lack of universal quantitive parameters for their parts. PAMDB includes parameter values inferred from hundreds of papers under a single unified model. PAMDB can be used to rank parts based on their activity, acquire parameter values and quantify their uncertainty, simulate circuits and visualize the part universe for the published synthetic biology work so far

### Dependencies
* Gcc
* Matlab


### Installation
```
git clone https://github.com/ameenetemady/MyCommon.git
git clone https://github.com/DeepPep/DeepPep.git
```

### Running
* Step1: prepare a directory containing your input files (with exact names):

  * ```identification.tsv```: tab-delimeted file:  **column1**: peptide, **column2**: protein name, **column3**: identification probability
  * ```db.fasta```: reference protein database in fasta format.

* Step2: ```python run.py directoryName```

Upon completion, ```pred.csv``` will contain the predicted protein identification probabilities.

### Benchmark Datasets
There are [7 example datasets](https://github.com/DeepPep/public/tree/master/data) (used for benchmarking in the paper). Each dataset is generated from MS/MS raw files using TPP pipeline. For example, to run the [18Mix benchmark dataset](https://github.com/DeepPep/public/tree/master/data/18mix), simply run the following:

```
python run.py data/18Mix
```
### Support

If you have any questions about PAMDB, please contact Linh Huynh (huynh@ucdavis.edu) or Ilias Tagkopoulos (itagkopoulos@ucdavis.edu).

### Citation
 L. Huynh, and I. Tagkopoulos, “A Parts Database with Consensus Parameter Estimation for Synthetic Circuit Design”, ACS Synthetic Biology (2016) [\[link\]](https://pubs.acs.org/doi/abs/10.1021/acssynbio.5b00205)

### Licence
See the [LICENSE](./LICENSE) file for license rights and limitations (Apache2.0).

### Acknowledgement
This work was supported by the NSF CAREER grant #1254205 to IT.

