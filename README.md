# nonlinear-xmcd
Analysis code and data for 
[D. J. Higley *et al.*, "Femtosecond X-ray induced changes of the electronic and magnetic response of solids from electron redistribution,"  *Nature Communications* **10**, 5289 (2019)](https://www.nature.com/articles/s41467-019-13272-5).

This analysis was tested with python 3.6.5 and Jupyter notebook version 4.4.0.

The data analysis performed by the notebooks in this repository use preprocessed versions of the data recorded at LCLS. These preprocessed data files have only the important quantities for each X-ray shot, as extracted from the much bigger raw data. The preprocessed data files are in the /data/preprocessed directory.

## Installation

The analysis uses python 3 and Jupyter notebooks, and can be run with the Anaconda Distribution.

To install, download [Anaconda](https://www.anaconda.com/distribution/#download-section), then follow the instructions of the download page. The installation should take about 15 minutes.

## How to run analysis

To run the analysis, first download this repository.

The analysis is run with the notebooks in the /code/notebooks directory. To perform the analysis and save the results as a source data file, run the demag.ipynb notebook followed by the description_of_analysis.ipynb notebook. The resultant source data file is saved as source_data.xlsx in the /data directory. The notebooks can be run by entering jupyter notebook at a terminal, then navigating to the notebook, then clicking Kernel -> Restart & Run All from the menu at the top of the notebook viewer.

The analysis runs in about 5 minutes.

The expected source data file and notebook outputs are saved in the /expected-output directory.
