# SeismicSaliency_GP2018

## Paper

This directory contains all source codes and data used to produce the results published in the following papers:

M. Shafiq, T. Alshawi, Z. Long, and G. AlRegib, “The role of visual saliency in the automation of seismic interpretation,” Geophysical Prospecting, vol. 66, issue S1, pp. 132-143, Mar. 2018.

M. Shafiq, T. Alshawi, Z. Long, and G. AlRegib, "SalSi: A New Seismic Attribute For Salt Dome Detection," IEEE Intl. Conf. on Acoustics, Speech and Signal Processing (ICASSP), Shanghai, China, Mar. 20-25, 2016.

## How to Run

1. In Matlab, make this directory as the working directory.

2. Make sure you add following directory in your Matlab search path by 
   using "Set Path" and "Add with Subfolders" in the Matlab Home-Environment 
   tab: "...\Seismic_Data\"

3. Run the following Program from "SaliencyForSeismic" folder
   
   >> MainFile_determinSaliencyForSaltDomes.m
   
   to generate the Mat file "salMapAllSlices.mat" required to simulation, 
   already present in Mat Files folder.

4. Run the program:

   >> Main.m

   to reproduce all Figures and Tablesof the paper

## Contact Info

The code was written by 
```
Muhammad Amir Shafiq
Fulbright PhD Scholar
Graduate Research Assistant
School of Electrical and Computer Engineering
Georgia Institute of Technology, Atlanta, GA, USA.
```

If you have found any bugs or have any questions/suggestions, 
please contact

     amirshafiq@gatech.edu 
     amirshafiq@gmail.com

Last Modified: 9-Nov-2016

Copyright 2015
