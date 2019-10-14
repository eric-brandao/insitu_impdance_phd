# insitu_impdance_phd

This is a very simple repository to get you started with in situ impedance measurement.
There are two main files.

In insitu_simulation_example.m you can explore the simulation of an in situ impedance measurement. You will define a sound source, receiver, prescribe a surface impedance, calculate the acoustic field in a PU sensor at a given position and retrieve the surface impedance of the sample.

In insitu_measurement_example you can explore an actual measurement. There are two structs you need to import (Calibration-22-Sep-2008 and Measurement-22-Sep-2008-flamex). Thean, Correction factor, transfer functions, and the retrieved impedances are calculated. 
