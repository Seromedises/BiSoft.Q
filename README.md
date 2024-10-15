# BiSoft.Q
This repository contains matlab code for simplified modelling of the BiSoft.Q bi-directional Pneumatic Artificial Muscle.

For further information, please refer to the following scientific paper:

- Colucci G. et al, A Bi-Directional Linear Artificial Muscle for Assistive Robotic Systems, [...].


https://github.com/user-attachments/assets/2b0d14ee-fb9c-4b50-8b4b-db9917180dd8

# How to use it

The code is developed in Matlab 2023b. Installation of the Parallel Computing Toolbox is required. The modelling library is contained within the BiSoftQ_Class folder (both types A and B). The "DesignScript.m" file explains how to model a BiSoft.Q actuator by giving as imput the following parameters:
- Geometrical
  - $L$, Longitudinal (nominal) length of the actuator when printed ($x=$0)
  - $l_{v,p}$, $l_{c,p}$, the valley and crest longitudinal arc lengths
  - $R_{ie,p}$, the radius of the arc defined in the median cross-sectional plane of the BiSoft.Q
  - $R_{ve,p}$ and $R_{ce,p}$, the radii of fillet arcs in the terminal cross-sectional plane
  - $\displaystyle \frac{R_{im,p}}{R_{om,p}}$, the ratio of the inner valley and outer crest arcs defined in the median cross-sectional plane
  - $N_{p}$, the number of pleats
  - $R_{i,b}$, the bellows radius defined on the median cross-secitonal plane 
	- $\alpha_{b}$, the angle between two subsequent bellows folds
	- $R_{v,b}$ and $R_{c,b}$, the valley and crest fillet radii of the bellows \\ 
	- $N_b$, the number of bellows folds
  - $t$, thickness of the membranes
- Physical
  - $p_{s}$, supply gauge pressure. The same value is assumed for both chambers
  - actuator.Bellows.E, Young's modulus of the bellows material for its stiffness estimation
  - actuator.Bellows.nu, Poisson's ration of the bellows material for its stiffness estimation

https://github.com/user-attachments/assets/a846d4d5-e54a-4632-a96b-7a1f3173f108

