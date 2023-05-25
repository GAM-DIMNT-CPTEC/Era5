# Era5 

Conjunto de scripts para baixar e tratar a reanálise do Era5.

* `get_era5_single_level_parameters.py`: script para baixar as variáveis de superfície (2d), escreve as variáveis e período em um único arquivo grib;
* `get_era5_pressure_levels_parameters.py`: script para baixar as variáveis atmosféricas (3d), escreve as variáveis, níveis e período em um único arquivo grib;
* `merge_atm_sfc_Era5.ipynb`: notebook que mostra como separar os arquivos grib baixados em 1 arquivo grib por tempo (utiliza algumas bibliotecas do Python - Xarray, Zarr e cfgrib para tratar os arquivos em datasets).

**Nota:** para a utilização dos scripts `get_era5_single_level_parameters.py` e `get_era5_pressure_levels_parameters.py`, é necessária a configuração da [CDS API](https://cds.climate.copernicus.eu/api-how-to).
