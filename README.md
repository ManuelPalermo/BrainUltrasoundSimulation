# BrainUltrassoundSimulation
Simulation of Ultrassound signals on a Brain3D model using K-Wave


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
no. 2, p. 021314, 2010.<br><\br>
[2] B. E. Treeby, J. Jaros, A. P. Rendell, and B. T. Cox, "Modeling 
nonlinear ultrasound propagation in heterogeneous media with power law 
absorption using a k-space pseudospectral method," J. Acoust. Soc. Am., 
vol. 131, no. 6, pp. 4324-4336, 2012.<br><\br>
[3] B. E. Treeby and B. T. Cox, "Modeling power law absorption and 
dispersion for acoustic propagation using the fractional Laplacian," J. 
Acoust. Soc. Am., vol. 127, no. 5, pp. 2741-2748, 2010.
