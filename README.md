# Nuclear Data Analysis Starter Pack

This repository is dedicated to the development of a simple Python application (at least at the moment) with useful tools and a GUI (Graphical User Interface) for scientific data analysis (mainly, but not exclusively, in the area of Nuclear Physics, due to the author's needs).

For now, only a few simple features are available (again, due to the author's needs), such as:
- Open a spectrum. This spectrum will be plotted in a MatPlotLib figure, along with the respective background and the spectrum without the background.
- Estimate the background of the spectrum with the SNIP (Statistics-sensitive Non-linear Iterative Peak-clipping) algorithm. Since this algorithm depends on some parameters, it is possible for the user to change said parameters to obtain the desired background estimation.
- Select a region of interest (ROI) of the spectrum and to calculate the number of counts of that ROI. It is possible to select whether this total number of counts concerns the raw spectrum or the spectrum without the background.

If you see potential in this project (even if at the moment, it doesn't support your needs in scientific data analysis), I encorage you to contribute with suggestions and discussions (and/or pull-requests if you have the time and patience for it).

## Plans for the future

As it was said before, the project is in an early-development phase and the idea is to extend the range of application of this application to other types of analysis and scientific areas. For this reason, there is already some ideas of what I would like to add to the project (technical and interface-wise).

### Future technical developments
- Write a user-guide document, where the implementation of the technical details will be explained, in addition to the "how to use" section
- Add data fitting methods
- Add other baseline removal algorithms
- Add a Nuclear Physics menu where it will be possible to calculate things such as Rutherford cross-sections and the kinematics of nuclear reactions (in order to know where to find some peaks in the spectrum), from the configuration of the experimental setup (angle of detection, calibration curve, etc...). This point can be extended to other areas as well

### Future interface developments
- Add support to open multiple spectra at once and an easy way to select the spectra showed in the figure
- Add drag and drop support
- Add an option to save a desired spectrum (the button already exists, but it doesn't work)
- Add support for different file formats (if you try to use the application, it is possible that you will need to rewrite the raw data file)
- Redo the interface organization when more tools are added (since not everyone will need to use SNIP or to simply calculate counts of a region). To facilitate the access to the tools needed, these will probably be grouped in different "modes" that may be easily selected at startup.
- Make the interface prettier :)

## Python Libraries Used
- Numpy
- Matplotlib
- Pybaselines
- PyQt5 (I know there is already a PyQt6, but I found more documentation of PyQt5 so I sticked to it. If a good reason to migrate appears, I will do so. If you want migrate yourself, I invite you to do a pull-request afterwards)

## References
[1] Charles R. Harris, K. Jarrod Millman, Stéfan J. van der Walt et al. “Array programming with NumPy”. In: Nature 585.7825 (set. de 2020), pp. 357–362. doi: 10.1038/s41586- 020- 2649- 2. URL: https://doi.org/10.1038/s41586-020-2649-2.

[2] J. D. Hunter. “Matplotlib: A 2D graphics environment”. In: Computing in Science & Engineering 9.3 (2007), pp. 90–95. doi: 10.1109/MCSE.2007.55

[3] Donald Erb. pybaselines: A Python library of algorithms for the baseline correction of experimental data.
Version 1.0.0. doi: 10.5281/zenodo.5608581. URL: https://github.com/derb12/pybaselines.