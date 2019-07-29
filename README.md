# nonlinear-xmcd
Data analysis code and data for manuscript titled "Femtosecond X-ray induced changes of the electronic and magnetic response of solids from electron redistribution" by D. J. Higley and others. This analysis was tested with python 3.6.5.

The data analysis performed by the notebooks in this repository uses 'preprocessed' versions of the data recorded at LCLS. These 'preprocessed' data files save only the important quantities for each X-ray shot. Thus, the 'preprocessed' data files are much smaller than the raw data files. The 'preprocessed' data files are in the /data/preprocessed directory.

## Installation

The analysis uses python 3 and Jupyter notebooks, and can be run with the Anaconda Distribution.

To install, download [Anaconda](https://www.anaconda.com/distribution/#download-section), then follow the instructions of the download page.

## How to run analysis

The analysis is run with the notebooks in the /code/notebooks directory. To perform the analysis and save the results as a source data file, run the demag.ipynb notebook followed by the description_of_analysis.ipynb notebook. The source data file is saved as source_data.xlsx in the /data directory.
