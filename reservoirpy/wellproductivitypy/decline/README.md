# WellLogs - Python

## Introduction.
This package is being designed among others to provide Petroleum Engineering tools in a modern programming language. This package is part of the project 7G which  proposes to make basic but powerful engineering software packages that cover the main topics of the Oil and Gas development phases which could be applied to any case study by suitable engineers.

There are five topics in which the project is going to be focused on:

<br>-Geoscience* (Current Package)
<br>-Reservoir
<br>-Production
<br>-Economics
<br>-Integration

<br> The package will always be in permanent development and open to suggestions to enhance the program. As the code has been written so far by a code enthusiastic Petroleum Engineer I hope to learn as much as possible to get better and usefull programs.

## WellLogspy Description
WellLogspy is a package to visualize Oil And Gas Well logs in order to make interpretetion and petrophysics analysis via Matplotlib package.

<br> The visualization of logs is based on '''tracks''' subpackage. Each log track is implemented through a single function that plot the desired data. For example, the Gr Track is composed generally by GammaRay and Sp logs, so grtrack function can plot the GammaRay and SP logs with optional features like add Well tops (Formations or Units), Gr clean and Gr shale to Vshale estimations, etc...

<br> The petrphysics analysis is contained in the '''petrophysics''' subpackage which allows you to estimate the basic petrophysics properties such like Shale volume, effective porosity, water saturation, permeability and flow capacity. A Montecarlo simulation can be implemented for uncertainty analysis on the original resources in place. 