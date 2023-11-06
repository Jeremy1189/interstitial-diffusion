# Purpose:
**This code is shared to help understand the methods in the research paper below or potentially apply the methods to explore the mechanism of interstitial sluggish diffusion and chemically-bias diffusion in other CSAs.**
***
## Research Title: Sluggish and Chemically-Biased Interstitial Diffusion in Concentrated Solid Solution Alloys: Mechanisms and Methods
***
### Research Group: J.J. Kai's Group in MNE, CityU
### Contact email: biaoxu4@cityu.edu.hk 
***
## Abstract:
![image](https://github.com/Jeremy1189/interstitial-diffusion/assets/85468234/7f15c182-cffb-4c76-9d76-bcabb4fcb671)

Interstitial diffusion is a pivotal process that governs the phase stability and irradiation response of materials in non-equilibrium conditions. In this work, we study sluggish and chemically-biased interstitial diffusion in Fe-Ni concentrated solid solution alloys (CSAs) by combining machine learning (ML) and kinetic Monte Carlo (kMC), where ML is used to accurately and efficiently predict the migration energy barriers of interstitials on-the-fly. By comprehensively analyzing the diverse migration patterns of interstitial dumbbells, we find that the observed sluggish diffusion and the "Ni-Ni-Ni" biased diffusion in Fe-Ni alloys are ascribed to a unique "Barrier Lock" mechanism, whereas the "Fe-Fe-Fe" biased diffusion is influenced by a "Component Dominance" mechanism. Additionally, a practical AvgS-kMC method, that relies only on the mean of migration patterns’ energy barriers, and an inverse design approach, that optimizes sluggish diffusion properties, are proposed for efficiently and conveniently exploring the interstitial-mediated diffusion behaviors and properties in CSAs.

# Guide for using Code
***
**Description:**
***
**Running Platform:** All the code is developed based on MATLAB. So if you want to use the code for your own purpose, please make sure that you installed MATLAB already. (The platform of [Online Matlab](https://www.mathworks.com/products/matlab-online.html) or [Octave](https://octave.org/) may also be used to run this code well (not tested yet) ). In addition, the parallel skills are used in this code. We strongly recommend running the code on a machine that has multi-CPU cores or high-performance GPU modules, better to run them on the computer clusters (The reference shell scripts are also shared).  
***
**Code structure：** There are three folders: **train_ML_model, AvgS_ML_kMC_constant_att_fre, AvgS_ML_Inverse_design**. They are used for different purposes. 
***
**"train_ML_model:"** This folder contains the code for training the ML model. The database for ML training was obtained through LAMMPS based on the interatomic potential of the embedded-atom method (EAM) type developed by Bonny et al. This potential has been widely applied to simulate the defect properties in Ni-containing CSAs, and it can reproduce reasonable and consistent results with density functional theory (DFT) calculations 19,24,37. The initial interstitial dumbbell configuration was created by inserting an atom to form a [100] dumbbell with a lattice atom located in the center of a 10×10×10 face-centered cubic (FCC) lattice. Different compositions of FexNi1-x were considered. Note that [100] dumbbells are the most stable interstitial form in FCC Ni. The migration of the [100] dumbbell to a [001] dumbbell through rotation in the {100}  plane was then simulated through the NEB method.   

The input of the ML model was the local atomic environment(LAE), described by the elemental types in the NN shells around the dumbbell interstitial, as the “INPUT” part shown in the Abstract figure. In this work, the LAE contains atoms within the 1st NN to the 10th NN based on our test. This LAE range is significantly larger than in vacancy cases, consistent with the greater spatial influence of interstitial dumbbells. The output labels were E_b, as shown by the “OUTPUT” part in the Abstract figure. 

![image](https://github.com/Jeremy1189/interstitial-diffusion/assets/85468234/4ab93fc3-c23e-48d9-a3a3-b19d1e44cd37)

**"ML_kMC_constant_att_fre:"** This folder is for showing how the ML-KMC method is applied with constant frequency. Since the first step is to determine the attempt frequency under different components, we need first run the ML-kMC under with constant attemp frequency,and then we use the Huang, W. and X.-M. Bai (2022) method to obtain the attempt frequency under different components. With the obtained attempt frequency, we just need to modify the following code lines:
```matlab
D0_set = [your setting attempt frequency];
D0 = D0_set(num_ratio);
```
![image](https://github.com/Jeremy1189/interstitial-diffusion/assets/85468234/e6736047-dd0b-4cf7-9e85-3deaa09febc0)

**AvgS_ML_Inverse_design:"** This folder is used to show how to use the AvgS_ML method to do an inverse design with minimum or maximum sluggish diffusion properties.
![image](https://github.com/Jeremy1189/interstitial-diffusion/assets/85468234/12f2cb07-83e9-4d6b-8b01-742fba865668)


