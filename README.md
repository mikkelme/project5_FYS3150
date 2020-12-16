# FYS 3150 - Project 5
This repository contains the work of Fredrik Hoftun and Mikkel Metzsch Jensen for project 5 in Computational Physics FYS3150.
As our fifth and final project we chose the project exploring the heat equation.

The repository is divided into two folders: `article` and `code`. 

* In the article folder we have the article itself and all the figures included in it. 
* In the code folder we have the main code main.cpp with the additional solver.cpp and solver.h used by the main file. The version of the main.exe outputs a .txt file named after the number of spatial dimensions, the numerical method used, the time-step dt used. For example `1DCN0.1.txt` has 1 spatial dimension, uses the Crank-Nicolson scheme with a time-step `dt = 0.1`. The .txt file includes `N`, number of steps, `Time`, the total time, `dx`, the spatial step-size and `dt`, the time-step, and then it writes the solution for all spatial dimensions for all time-steps. There are several Python-scripts for plotting the computed values. There are two folders `error_plot_data0.01, error_plot_data0.1` that contains an automated Python-script for error plotting. The `Geo` folder contains everything about the second part of the project, the geoscience part. Here there are similar files to the main code directory. In addition Gsolver.cpp contains a method that only writes the last state. The Gmain.exe outputs a .txt file similarly to the one in the main code directory.

To use the compiler just type `make` in the commandline in the `/code` or `/code/Geo` directory.
The compiler will generate a file `main.exe` which takes the following arguments in the command line: time-step size, total time. The `Gmain.exe` takes the arguments: step-size (in km), total time (in billions of years, Gy)
