#!/bin/bash
g++ zad1.cpp -lgsl -lgslcblas -o zad1
./zad1
# python3 zad1.py
gnuplot plottemp.gp
gnuplot plotlap.gp
convert -delay 100 -loop 0 res/laplacian_*.png res/laplacian.gif
convert -delay 100 -loop 0 res/temperature_*.png res/temperature.gif