# Brain Ultrassound Simulation
Simulation of Ultrassound signals on a 3D Brain model using [K-Wave MATLAB toolbox](http://www.k-wave.org/) .<br>

### Brain Model
The brain model was downloaded from [scalablebrainatlas](https://scalablebrainatlas.incf.org/human/NMM1103).<br>
A rough model with the skull can be created using the CreateSkull.m script, which uses morfological dilation operations.<br>
The simulation medium proprieties can be set using masked operations using the values from the brain model.
![BrainModel](https://i.imgur.com/IaexuvH.png)<br><br>


### Focusing Algorithm
A simple focusing algorithm was created in order to assert wether or not ultrassound signals could be focused, from an array of transducers, on a brain region, through the skull. <br>

The algorithm works by sending an ultrassound pulse from the target, which reaches the transducers at diferent points in time. The diferent travel times are then used to calculate the delays between the ultrassound transducers. In practice, due to interface reflections, the same signal can reach a transducer multiple times, so all but the maximum value are ignored. The algorithm assumes the same travel path in both directions, which is not true in most cases, despite this, it achieves satisfatory results.
![FocusingAlgorithm](https://i.imgur.com/cCXyAqC.png)<br>

The transducers can be "placed" on the top of the skull/brain, by only to defining the number and spacing between array elements.<br>
![TransducersFocus](https://i.imgur.com/P9sezeY.png)<br><br>


### Results
The maximum signal pressure was recorded at every point in the simulation array in order to determine the highest intensity points. VolumeViewer was then used to better visualize the results.
![Results](https://i.imgur.com/ULS00cR.png)<br><br>




### K-Wave Toolbox
Written by Bradley Treeby, Ben Cox, and Jiri Jaros

k-Wave is an open source MATLAB toolbox designed for the time-domain 
simulation of propagating acoustic waves in 1D, 2D, or 3D [1]. The toolbox
has a wide range of functionality, but at its heart is an advanced numerical
model that can account for both linear and nonlinear wave propagation, an 
arbitrary distribution of heterogeneous material parameters, and power law 
acoustic absorption.

The numerical model is based on the solution of three coupled first-order 
partial differential equations which are equivalent to a generalised form 
of the Westervelt equation [2]. The equations are solved using a k-space 
pseudospectral method, where spatial gradients are calculated using a 
Fourier collocation scheme, and temporal gradients are calculated using a
k-space corrected finite-difference scheme. The temporal scheme is exact in
the limit of linear wave propagation in a homogeneous and lossless medium, 
and significantly reduces numerical dispersion in the more general case.

Power law acoustic absorption is accounted for using a linear integro-
differential operator based on the fractional Laplacian [3]. A split-field 
perfectly matched layer (PML) is used to absorb the waves at the edges of 
the computational domain. The main advantage of the numerical model used in 
k-Wave compared to models based on finite-difference time domain (FDTD) 
schemes is that fewer spatial and temporal grid points are needed for 
accurate simulations. This means the models run faster and use less memory. 
A detailed description of the model is given in the k-Wave User Manual and 
the references below.

[1] B. E. Treeby and B. T. Cox, "k-Wave: MATLAB toolbox for the simulation 
and reconstruction of photoacoustic wave-fields," J. Biomed. Opt., vol. 15,
no. 2, p. 021314, 2010.<br>
[2] B. E. Treeby, J. Jaros, A. P. Rendell, and B. T. Cox, "Modeling 
nonlinear ultrasound propagation in heterogeneous media with power law 
absorption using a k-space pseudospectral method," J. Acoust. Soc. Am., 
vol. 131, no. 6, pp. 4324-4336, 2012.<br>
[3] B. E. Treeby and B. T. Cox, "Modeling power law absorption and 
dispersion for acoustic propagation using the fractional Laplacian," J. 
Acoust. Soc. Am., vol. 127, no. 5, pp. 2741-2748, 2010.
