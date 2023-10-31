# Purpose:
**This code is shared to help understand the methods in the research paper below or potentially apply the methods to explore the mechanism of interstitial sluggish diffusion and chemically-bias diffusion in other CSAs.**
***
## Research Title: Mechanism of  Sluggish and Chemically-Biased Interstitial Diffusion in Concentrated Solid Solution Alloys: “Barrier Lock” and “Component Dominance” 
***
### Research Group: J.J. Kai's Group in MNE, CityU
### Contact email: biaoxu4@cityu.edu.hk 
***
## Abstract:
![image](https://github.com/Jeremy1189/interstitial-diffusion/assets/85468234/a8b7784d-3422-402e-aa04-31aaf2a8a066)

Interstitial diffusion is a pivotal process that governs the phase stability and irradiation response of materials in non-equilibrium conditions. In this work, we study sluggish and chemically-biased interstitial diffusion in Fe-Ni concentrated solid solution alloys (CSAs) by combining machine learning (ML) and kinetic Monte Carlo (kMC), where ML is used to accurately and efficiently predict the migration energy barriers of interstitials on-the-fly. By comprehensively analyzing the diverse migration patterns of interstitial dumbbells, we find that the observed sluggish diffusion and the "Ni-Ni-Ni" biased diffusion in Fe-Ni alloys are ascribed to a unique "Barrier Lock" mechanism, whereas the "Fe-Fe-Fe" biased diffusion is influenced by a "Component Dominance" mechanism. Additionally, a practical AvgS-kMC method, that relies only on the mean of migration patterns’ energy barriers, and an inverse design approach, that optimizes sluggish diffusion properties, are proposed for efficiently and conveniently exploring the interstitial-mediated diffusion behaviors and properties in CSAs.

# Guide for using Code
***
**Description:**
***
**Running Platform**
All the code is developed based on the MATLAB. So if you want to use the code for your own purpose, please make sure that you installed MATLAB already. (The platform of [Online Matlab](https://www.mathworks.com/products/matlab-online.html) or [Octave](https://octave.org/) may also be used to run this code well (not tested yet) ). In addition, the parallel skills are used in this code. We strongly recommend running the code on a machine that has multi-CPU cores or high-performance GPU modules, better to run them on the computer clusters (The reference shell scripts are also shared).  
***



