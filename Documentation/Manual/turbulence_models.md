## Following turbulence models are currently implemented:

* Linear eddy-viscosity k-ε models
  * Standard high-Re
  * Low-Re version (Abe, Kondoh and Nagano) with compound wall treatment
* Elliptic-relaxation eddy-viscosity model
  * Linear ζ − f with compound wall treatment
* LES
  * Smagorinsky SGS model
  * Dynamic Smagorinsky SGS model
  * WALE SGS model
* Hybrid LES/RANS
  * Detached eddy simulation (Spalart-Allmaras)
  * Hybrid LES/RANS ζ − f model
* Conventional and advanced treatment of wall boundary conditions (wall integration, wall functions with roughness, compound wall treatment, blending wall function and wall integration approach)

All RANS models can be ran in unsteady mode, thus effectively becoming Unsteady RANS (URANS) model.
