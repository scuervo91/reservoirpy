# Reservoirpy - Python

Package to perform Oil and Gas reservoir Analysis in Python. 

## Introduction.
This package is being designed to provide Petroleum Engineering tools in Python. It is part of a project which proposes to make basic but powerful engineering software packages that cover the main topics of the Oil and Gas development phases which could be applied to any case study by suitable engineers.

Currently 'reservoirpy' contains (so far) the next modules:

* `wellpy` : Module that contains the Well python object definition. A well object have attributes like coordinates, elevations, perforations, formation tops and survey. It can calculate the the full survey (`tvd, tvdss, offset, coordinates, shape points`) by givin the `md,inc,azi`. It is done through the [Wellpathpy package](https://github.com/Zabamund/wellpathpy)
(Under development)

* `welllogspy` : Module that allows to make Well Logs visualization, interpretation petrophysical estimation. It allows you to visualize you `.las` files in a fully customizable way, perform petrophysical calculations, analysis and export the results. It is built on top of [Matplotlib package](https://matplotlib.org/) to visualize and [Lasio package](https://github.com/kinverarity1/lasio) to handle the .las files. (Under development)

* `wellproductivitypy`: Module that allows to perform well productivity, decline curve analysis and forecast. A declination object can be defined from a real production data and enables to make easily forecast ultil a defined time or economic limit. Make oil inflow object to plot Oil productivity index based on reservoir properties and get values of flow, bottom hole pressure and drawdown on different conditions. Artifitial lift systems performance under development (ESP, Jet Pump). (Under development)

* `pvtpy`: Module that allows you to build black oil pvt tables based on commonly used correlation. Initially, the correlation used are those published on [Correlaciones Numéricas PVT, Carlos Banzer] (Under development)

The package currently is in permanent development and open to suggestions to enhance the program.  As the code has been written so far by a code enthusiastic Petroleum Engineer I hope to learn as much as possible to get better and usefull programs.

## Instalation

So far there is no a final version to download through [PyPi](https://pypi.org/) however it can be downloaded and used in two ways:

1. Clone the repository
* Clone the repository wherever you want in your pc.
```
git clone https://github.com/scuervo91/reservoirpy.git
```

* When using the python script, you must add the repository to `sys.path` to be able to import it
```python
import sys
sys.path.append('path_to_repo/reservoirpy')
```

2. Install on your environment 
You can install on you environment either `conda` or `virtualenv`. By this way you can import the repo directly
```
pip install git+https://github.com/scuervo91/reservoirpy.git
```

## Documentation

* [Reservoirpy](https://scuervo91.github.io/reservoirpy/)
    *[Welllogspy](https://scuervo91.github.io/reservoirpy/welllogspy)
    *[Wellpy](https://scuervo91.github.io/reservoirpy/wellpy)
    *[Wellproductivitypy](https://scuervo91.github.io/reservoirpy/wellproductivitypy)
    *[Pvtpy](https://scuervo91.github.io/reservoirpy/pvtpy)
