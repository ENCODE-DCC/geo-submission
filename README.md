# geo-submission

## Installation
1. Clone this repository and `cd` into it.
2. Make sure you have `conda` installed.
3. Run the following command: `conda env create -f environment.yml` . This will create a new `conda` environment called `geo-submission` in the location `/anaconda3/envs/geo-submission`

## Dependencies
* Python >= 3.4
  * Using previous Python versions may add trailing whitespace to the json outputs, resulting in spurious diffs
* pandas 0.24.2 (`conda install pandas==0.24.2`)
  * Newer versions have removed pandas.compat.StringIO
* requests (`conda install requests`)
