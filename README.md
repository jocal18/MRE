# MRE
Repo for numerical Implementation of Magneto Rheological Elastomer simulations
Tools : 
It is assumed that the user knows Python programming language and Anaconda Virtual Environement
VS Code 
While it is not absolutely required, it is a good IDE and it makes reading/writtng/running the code much easier. But you can still do everything in the shell/terminal. 
Useful ressources: 
Gmsh : Open source 3D finite element mesh generator with a built-in CAD engine and post-processor
      Note for gmsh : Pip install gmsh
      *The gmsh package in conda-forge don't provide the Python API. Please install gmsh using pip install gmsh. If you still prefer to use conda, you can download the python-gmsh package from the conda-forge channel using conda install -c conda-forge gmsh python-gmsh. Here you will find the Python API you're looking for.  https://stackoverflow.com/questions/70947216/python-cant-use-gmsh-after-installing-via-pip-and-conda
Gmshmodel : Gmshmodel provides a Python-API with which all the capabilites of Gmsh can be used within Python
FEniCS: Open-source computing platform for solving partial differential equations (PDEs) with the finite element method (FEM)
•How to install FEniCSx/DOLPHINx
•https://github.com/FEniCS 
Paraview.  Opensource post-processing visualization engine (.xdmf)

Research paper # 1 : A computational framework for magnetically hard and soft viscoelastic magnetorheological elastomers
Paper : https://doi.org/10.1016/j.ijsolstr.2013.08.024
Code : https://zenodo.org/records/5543516
Related files - The hyperlink leads to the original file (reference). The actual file is contained in the folder provided: 
Main simulation code : / MagnetoVisco/
fenics_magneto_visco_final.py: 3D full-field magneto-mechanical  homogenization that accounts for: finite strains, viscoelasticity, incompressibility and mixed stress/magnetic flux control.
aux_homog_nl.py : Process data. Save data as .xdmf for data visualization with paraview
Geometry generation code : /GeometryGeneration/
simpleCubicCell3DSphere.py: Generates a RVE with periodic spheres inclusion. Sphere radius determined based on desired volume fraction. If "periodicityFlags": [0, 0, 0] à No periodic spheres.
randomInclusion3DSphere.py: Generates a RVEwith a randomly positioned N spheres of radius R. volume fraction and R will determine the number of spheres N.
bodyCenteredInclusions.py
Provided geometry to use with : /model/
s30.xmdf : original geometry from research paper. RVE with centered sphere. Volumic Fraction : 0.3
model-RI-VF-07-N_2-R_20.xdmf : Random Inclusion (RI). Volumic Fraction (VF) 0.07. Number of spheres (N) : 2. Sphere radius (R) 0.2.
model-RI-VF-017-N_5-R_20.xdmf
model-RI-VF-027-N_8-R_20.xdmf

Research paper #2 : Influence of magnetic boundary conditions on the quantitative modelling of magnetorheological elastomers
Paper : https://doi.org/10.1016/j.mechmat.2023.104742 
Code : https://zenodo.org/records/8129310  
Related files: 
Main simulation code : fenics_magneto_hMRE_AirDomain.py (main script),
Geometry generation code :
GenerateCylinderWithAIr.py
GenerateXDMFH5.py : will generate .xdmf and .h5 file used in the main code
Provided geometry to use with :
/model/cylinder_with_air_domain.xdmf 
/model/cylinder_with_air_domain_bulky.xdmf
Results from this simulation: 
/results/AirDomain/cylinder_with_air_domain8.xdmf
/results/AirDomain/cylinder_with_air_domain9.xmdf
