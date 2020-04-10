---
layout: page
title: "Welllogspy"
permalink: /welllogspy/
---

# Welllogspy - Module

Welllogspy is a module of Reservoirpy package to visualize Oil And Gas Well logs in order to make interpretetion and petrophysics analysis. 

It supports the reading and exporting `.las` files throughout the [Lasio package](https://github.com/kinverarity1/lasio) as well as the the visualization feature with the [Matplotlib package](https://matplotlib.org/). It is included basic equations for petrophysical calculation such as VShale, posority, Water saturation, permeability, flow capacity etc. 

## Getting Started

### Define Well attributes
The recomended workflow starts with defining the main attributes of a well throught the module [wellpy](https://scuervo91.github.io/reservoirpy/wellpy). This module contains a `class` named `reservoirpy.wellpy.path.wells` in which can be attached the basic attributes such as name, coordinates, elevation, formation tops and perforations. 

Let's define the first attributes of a well
```python
from shapely.geometry import Point

#Create the well object
name  = 'well-011'
rte   = 458.2 # Rotary table Elevation in ft
surf_coord = Point(872565.89,1055067.55) # Surface Coordinates defined as Point  
crs = 'EPSG:3117' # Coordinates reference System of the surface Coordinates
deviation = pd.read_csv('deviation.csv', header=[0])   #
print(deviation.head())

      md   inc    azi
0    0.0  0.00    0.0
1  193.0  0.06    0.0
2  375.0  0.12    0.0
3  559.0  0.33  347.6
4  651.0  0.36  340.5
```
Let's define some formations tops with the `class` `reservoirpy.wellpy.path.tops` which is a subclass of a `GeoDataFrame` to support spatial and geographic information. Hence it can be defined with the same sintaxis of `pandas.DataFrame` 

```python
from reservoirpy.wellpy import path as ph

#Tentative Formations Tops
tops = ph.tops([
            ['fm1',11289,11317,20,170,1.5],
            ['fm2',11318,11965,10,170,2],
            ['fm3',11966,12080,50,190,0.12],
            ['fm4',12081,12222,50,190,0.12],
            ['fm5',12223,12233,50,190,0.1],
            ['fm6',12234,12520,15,170,0.1]], 
            columns=['formation','md_top','md_bottom','gr_sand','gr_shale','rw'])
print(tops)

  formation  md_top  md_bottom  gr_sand  gr_shale    rw
0       fm1   11289      11317       20       170  1.50
1       fm2   11318      11965       10       170  2.00
2       fm3   11966      12080       50       190  0.12
3       fm4   12081      12222       50       190  0.12
4       fm5   12223      12233       50       190  0.10
5       fm6   12234      12520       15       170  0.10
```

Let's import some logs files in `.las` format with the class `reservoirpy.welllogspy.log.log` which is a subclass of `lasio.LASFile`. This package allows you to easily manipulate `.las` files. The full documentation of `lasio` can be checked [here](https://lasio.readthedocs.io/en/latest/). By adding the keyword arg `find_mnemonics=True` it search the mnemonics in the `.las` file in a mnemonics list with the purpose of classifing and modify its name by adding a prefix. 

```python
from reservoirpy.welllogspy.log import log 

#Upload logs
masterlog = log('master.las')
openlog = log('FE_BORAL-2_RES_GR_NUKES 27MAR2020_RC_MD.las', find_mnemonics=True)

Mnemonic:  ROPA  =>  ROP_ROPA
Mnemonic:  AFRMC  =>  RM_AFRMC
Mnemonic:  AFRSC  =>  RS_AFRSC
Mnemonic:  AFRDC  =>  RD_AFRDC
Mnemonic:  ALCDLC  =>  RHOB_ALCDLC
Mnemonic:  ALDCLC  =>  DCOR_ALDCLC
Mnemonic:  ALPELC  =>  PEF_ALPELC
Mnemonic:  TNPS  =>  NPHI_TNPS
```
Here `log` found some aliases and added a prefix for easily indentify the curves. For example it changed the name of `ALCDLC` to `RHOB_ALCDLC` 

If you have a list of mnemonics you could provide one.

```python
# list of mnemonics 
mn = pd.read_csv('mnemonics.csv')
openlog = log('FE_WELL-1_RES_GR_NUKES_RC_MD.las', find_mnemonics=True, mnemonics=mn)
```

### Create the well object

Now you can create the well object ```reservoirpy.wellpy.path.well```

```python
#Create the well
w1 = ph.well(name='well-1', 
             rte=rte, 
             surf_coord=surf_coord, 
             crs=crs, 
             deviation=deviation, 
             tops=tops, 
             masterlog=masterlog, 
             openlog=openlog)
```

Just by adding a deviation in a `pandas.DataFrame` with the columns `md` `inc` `azi` it automatically estimates the whole trajectory and saves it as a `geopandas.GeoDataFrame`. 

```python
print(w1.survey.head())

        inc    azi         tvd       tvdss  northing_off  easting_off  \
md                                                                      
0.0    0.00    0.0    0.000000  458.200000      0.000000     0.000000   
193.0  0.06    0.0  192.999965  265.200035      0.101055     0.000000   
375.0  0.12    0.0  374.999732   83.200268      0.386939     0.000000   
559.0  0.33  347.6  558.998219 -100.798219      1.097142    -0.113784   
651.0  0.36  340.5  650.996554 -192.796554      1.628348    -0.267154   

           northing        easting      dleg  \
md                                             
0.0    1.055068e+06  872565.890000  0.000000   
193.0  1.055068e+06  872565.890000  0.031088   
375.0  1.055068e+06  872565.890000  0.032967   
559.0  1.055068e+06  872565.855319  0.116497   
651.0  1.055068e+06  872565.808571  0.056709   

                                          geometry  
md                                                  
0.0                   POINT (872565.89 1055067.55)  
193.0          POINT (872565.89 1055067.580801428)  
375.0          POINT (872565.89 1055067.667939113)  
559.0  POINT (872565.8553186734 1055067.884408979)  
651.0  POINT (872565.8085714115 1055068.046320477)  
```


