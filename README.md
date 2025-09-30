# Learning Dedalus

[Dedalus](https://dedalus-project.org) is an open-source PDE solver based on modern spectral methods.
It provides a high-level Python interface for creating and solving global spectral discretizations of many PDEs.
It is primarily used for computational fluid dynamics (CFD), but can solve many other types of equations as well.
This repository gathers various tutorials and resources for learning how to use Dedalus for wide varieties of applications.

## Tutorial notebooks

The recommended starting point for learning Dedalus is to go through the [official tutorial notebooks](https://dedalus-project.readthedocs.io/en/latest/pages/tutorials.html#tutorial-notebooks) in the main code documentation.
These document many of the features of the code and explain the class structure and possible workflows using Dedalus objects.
The examples are designed to run nearly instantly on a single core / laptop.

## Example scripts

The main code documentation also contains a [range of example scripts](https://dedalus-project.readthedocs.io/en/latest/pages/tutorials.html#example-scripts) showing the typical usage of Dedalus in an HPC environment (e.g. executing model scripts rather than notebooks).
These examples are designed to run within several minutes on a single core / laptop, but can be easily scaled up to production sized simulations on a cluster.
Many users begin implementing their own models based on one of these scripts.

## Practical sessions for GAFD

The following material was originally developed for the practical numerical sessions for the 2023 CISM school on "Fluid Mechanics of Planets and Stars".
These lectures explore how to get started applying Dedalus to problems related to geophysical and astrophysical fluid dynamics.
You can view the notebooks online at the provided HTML links, open and execute them via Google Colab, or download them to run locally after [installing Dedalus on your computer](https://dedalus-project.readthedocs.io/en/latest/pages/installation.html).

* Lecture 1: Introduction to Spectral Methods & Dedalus
  [[slides]](https://raw.githubusercontent.com/kburns/cism_dedalus_2023/main/lecture_1_compressed.pdf)
* Lecture 2: Basic Dedalus API -- Burgers & KdV Equations
  [[nbviewer]](https://nbviewer.org/github/kburns/cism_dedalus_2023/blob/main/lecture_2_intro_to_dedalus.ipynb)
  [[colab]](https://colab.research.google.com/github/kburns/cism_dedalus_2023/blob/main/lecture_2_intro_to_dedalus.ipynb)
* Lecture 3: Forcing & Analysis -- 2D Turbulence
  [[nbviewer]](https://nbviewer.org/github/kburns/cism_dedalus_2023/blob/main/lecture_3_2d_turbulence.ipynb)
  [[colab]](https://colab.research.google.com/github/kburns/cism_dedalus_2023/blob/main/lecture_3_2d_turbulence.ipynb)
* Lecture 4: Eigenvalue Problems -- Shallow Water & Shear Instability
  [[nbviewer]](https://nbviewer.org/github/kburns/cism_dedalus_2023/blob/main/lecture_4_shallow_water_evp.ipynb)
  [[colab]](https://colab.research.google.com/github/kburns/cism_dedalus_2023/blob/main/lecture_4_shallow_water_evp.ipynb)
* Lecture 5: Boundaries & Tau Terms -- Rayleigh-Benard & Spherical Convection
  [[nbviewer]](https://nbviewer.org/github/kburns/cism_dedalus_2023/blob/main/lecture_5_convection.ipynb)
  [[colab]](https://colab.research.google.com/github/kburns/cism_dedalus_2023/blob/main/lecture_5_convection.ipynb)

[[view repository]](https://github.com/kburns/cism_dedalus_2023)
<br>
[[view webpage]](https://kburns.github.io/cism_dedalus_2023)
