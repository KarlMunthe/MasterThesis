# MasterThesis
This repo contains the code used to obtain the numerical results for the one dimensional [Navier-Stokes](https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations)
and [Svärds modified Navier-Stokes equations](https://www.researchgate.net/publication/322328860_A_new_Eulerian_model_for_viscous_and_heat_conducting_compressible_flow)
which will is referred to as the NSS equations.

The folder named DifferenceMatrices contains functions that will return spectral difference matrices approximating the first and second derivative for grids of 
an arbitrary amount of points.

The folder named AcousticCode contains the code needed to compare the Navier-Stokes-Svärd equations with Navier-Stokes equations. SpectralNS and SpectralNSS solve
the acoustic attenuation problem of the NS and NSS equations respectively.

The folder named BoundaryLayerCode coontains the code to solve the [Blasius ODE](https://en.wikipedia.org/wiki/Blasius_boundary_layer#Blasius_equation_-_first-order_boundary_layer) 
and the boundary layer for the two dimensional NSS equations.
