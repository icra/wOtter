# wOtterEU

General description, ...

## Setup Python environment
To create the environment with all packages needed, simply run
```
conda env create -f environment.yml
conda activate otter
pip install -e .
```
Then, download the compressed data [here](https://mega.nz/file/QTV11QgL#6sxcA-bBNfcroTg2EThsSVnhGBskYgNjX5bpahGG_3E), and unzip the file at the main directory.
Finally, create a folder at the main directory called ```results```.

## Run locally

Before starting to do anything, make sure you have activated the conda environment. Please, keep in mind that a computer with minimum 16GB of RAM is needed.

```
conda activate otter
```
### Add your own raster files and generate graph model from scratch (OPTIONAL)
Work in progress ...

### Generate pollutant map
 In the `notebooks/run_graph.ipynb`  notebook, you will find how to visualize the results of the simulation on a dynamic map. You can open the file by simply running
```
jupyter notebook
```
and selecting the file.
Run all the cells, and a  ```html ``` file will be created at the ```results/``` folder.  This file contains a map created with the [kepler.gl](https://kepler.gl/) library. This map is not finalized, meaning that you have to open the file with any browser and edit it to show all the information as you want. 
