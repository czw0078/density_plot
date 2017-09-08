# DensityPlot
An electron density difference calculation from orbital functions and visualization by Python Matplot for computational chemistry

## Gallery:
See results: http://www.auburn.edu/~czw0078/intr.php 

## Quick start:
1. CD to the Example1
2. type code
```bash
./DiffDen.py job1_polarized.molden job2.molden -j 0.10 -x "-4.0 4.0 0.03" -y "-4.0 4.0 0.03" --p1 "0.0 0.0 0.0" --p2 "1.0 0.0 0.0" --p3 "0.0 1.0 0.0" -d True
```

3. A PNG will gererated and show in pop-up window

## Tutorials:
Read the slides named Tutorial.pptx and it explains the functionality and theory of the whole program.

## Development:
The main entry point of the package is the DiffDen.py, all the functions needed by the main entry have been put into function_read_files.py and function_calculate_density.py
