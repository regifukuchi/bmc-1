# Introduction to modeling and simulation of human movement

### Instructor  
- [Marcos Duarte](http://demotu.org/people/marcos-duarte/)
- [Renato Watanabe](http://demotu.org/pessoal/renato/)

### Time and place  
- Tuesdays, 17-19h; Thursdays, 17-19h, lab A2 L003, campus São Bernardo do Campo

### Grades

* [Grades](https://docs.google.com/spreadsheets/d/e/2PACX-1vS6yAX5ZHkzhnij4b_lklDwWtV-0KrqTlEQv2W_X2b-w1woXNmBR6p0Mq-IcV51gw7y0EeaBnIC5Xf0/pubhtml)

## Lecture Schedule

### Lecture 1

[Introduction](ModSim2018_0.pdf)

 * Control of human movement, modeling and simulation  
 * Course Information   
 * [OpenSim software](https://simtk.org/projects/opensim), [Mendeley](https://www.mendeley.com), [Github](https://www.github.com), [Python](https://www.python.org/), [Anaconda](https://www.anaconda.com/) 

Tasks (for Lecture 3):

*  Install OpenSim and do the first three tutorials (Help menu).

Readings (for Lecture 5):

* TSIANOS, G. A .; LOEB, G. E. Muscle and limb mechanics. Comprehensive Physiology, v. 7, n. 2, p. 429-462, 2017.
* ZAJAC, F. E. Muscle coordination of movement: a perspective. J Biomech, v. 26 Suppl 1, n. SUPPL. 1, p. 109–124, 1993. 

### Lecture 2

 * [Python for scientific computing](http://nbviewer.jupyter.org/github/BMClab/bmc/blob/master/notebooks/PythonForScientificComputing.ipynb)
 * [Python Tutorial](https://nbviewer.jupyter.org/github/BMClab/bmc/blob/master/notebooks/PythonTutorial.ipynb)
 
Tasks (for Lecture 4):

- Write a Jupyter notebook to compute and plot the acceleration of the particle with displacement data contained in the file [pezzack.txt](http://isbweb.org/data/pezzack/index.html). Additionally, compare the computed acceleration with the true acceleration (contained in the same data).
 
Readings:

 * Shen, H. (2014). Interactive notebooks: Sharing the code. Nature, 5–6. Retrieved from http://europepmc.org/abstract/med/25373681  
 * Perkel, J. M. (2015). Programming: Pick up Python. Nature, 518(7537), 125–126. https://doi.org/10.1038/518125a
 
### Lecture 3

 * OpenSim tutorials
 
### Lecture 4

 * Python programming  
   + [Python Tutorial](https://nbviewer.jupyter.org/github/BMClab/bmc/blob/master/notebooks/PythonTutorial.ipynb)
 
 Tasks (for Lecture 6):
 
  * Plot the trajectory of a particle during 10 seconds with the gravity force acting on it (gravity acceleration g = -9.81 m/s^2). Consider that the initial velocity of the particle is v0 = 30 m/s and the initial angle with the horizontal plane is 30 degrees. 
 
 
### Lecture 5

 * Discussion about TSIANOS and LOEB (2017) and Zajac (1993)  
   + Animations about muscle contraction: [1](https://youtu.be/GrHsiHazpsw), [2](https://youtu.be/jqy0i1KXUO4), [3](https://youtu.be/s_TRsf6tJsc)  
  

### Lecture 6
 
 * Numerical methods to solve ordinary differential equations. [https://github.com/BMClab/bmc/blob/master/notebooks/newtonLawForParticles.ipynb](https://github.com/BMClab/bmc/blob/master/notebooks/newtonLawForParticles.ipynb),[https://github.com/BMClab/bmc/blob/master/notebooks/OrdinaryDifferentialEquation.ipynb](https://github.com/BMClab/bmc/blob/master/notebooks/OrdinaryDifferentialEquation.ipynb)
 
 Task for this class:
 
 * Write a Jupyter notebook to find the trajectory of a ball considering the air drag proportional to the square root of the ball velocity. Consider that the initial velocity of the particle is v0 = 30 m/s and the initial angle with the horizontal plane is 30 degrees, the gravity acceleration is g = 9.81 m/s^2, the mass of the ball is m = 0.43 kg and the air drag coefficient is b= 0.006 kg m^{1/2}/s{3/2}.
 
 Readings (for Lecture 7):
 
 * NIGG, B. M.; HERZOG, W. Modelling. In: Biomechanics of the Musculo-skeletal System.  p. 501 - 534.  
 
### Lecture 7

  * Discussion about NIGG, B. M.; HERZOG, W. Modelling. In: Biomechanics of the Musculo-skeletal System.  p. 501 - 534.
  
### Lecture 8

  * Models of Kelvin, Voigt and Maxwell.   
  
  Additional reading:
  
  * YAMAGUCHI, Y. T; Dynamic modeling of musculoskeletal motion: A Vectorized Approach for Biomechanical Analysis in Three Dimensions (2001), Sections 2.2.1 and 2.2.2.
  
  Readings (for Lecture 9):
  
  * NIGG, B. M.; HERZOG, W. Modelling. In: Biomechanics of the Musculo-skeletal System.  p. 631 - 634.  
  
 ### Lecture 9
 
  * Implement muscle model of Nigg and Herzog.
 
