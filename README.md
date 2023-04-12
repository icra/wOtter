# wOtterEU
<img width="889" alt="Captura" src="https://user-images.githubusercontent.com/49274979/211218395-eb230286-abc4-40c5-acee-9fafcd025839.PNG">

General description, ...


## Open example
This should probably be your first step. Go to ```examples/``` and open the .html file to visualize the map created with the current source code configuration. This map contains all the European WWTPs and the pollution level of the rivers.

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
### Use your own raster files and generate graph model from scratch (OPTIONAL)
If you want to perform all the steps of the simulation, modify the source code or change the raster files, you will have to download the data [here](https://mega.nz/file/JLdjTS4Z#2kgA3S60tn-bpMvbrulXUgN4Bae5RKR4KCS70Me4hH4) and unzip it inside your previously downloaded and unzipped ```data/``` folder (Warning: the unzipped file weighs 40GB, and the simulation takes about 3 hours). Then, from your root folder, run
```
python .\runall.py
```
Once the simulation is over, generate the pollutant map.

### Generate pollutant map
 In the `notebooks/run_graph.ipynb`  notebook, you will find how to visualize the results of the simulation on a dynamic map. You can open the file by simply running
```
jupyter notebook
```
and selecting the file.
Run all the cells, and some files will be created at the ```results/``` folder. These files have to be uploaded to the [kepler.gl](https://kepler.gl/) web page, which is the tool used to create this map. Then, customize your map following the ```instructions/``` folder.
