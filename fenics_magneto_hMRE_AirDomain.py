from __future__ import print_function
'''
 _____ _   _ _     _                                    
|  ___| | | | |   | |                                   
| |_  | | | | |   | |                                   
|  _| | |_| | |___| |___                                
|_|__ _\___/|_____|_____|       _  _____ ___ ___  _   _ 
/ ___|_ _|  \/  | | | | |      / \|_   _|_ _/ _ \| \ | |
\___ \| || |\/| | | | | |     / _ \ | |  | | | | |  \| |
 ___) | || |  | | |_| | |___ / ___ \| |  | | |_| | |\  |
|____/___|_|  |_|\___/|_____/_/   \_\_| |___\___/|_| \_|
'''
"""
    FEniCS code for the permanent magnet based on a vector-potential formulation
    of the magnetic field and accounting for finite strains and remanent magnetisation of MRE samples.

    Copyright (C) 2023: Miguel Angel Moreno-Mateos, Kostas Danas, Daniel Garcia-Gonzalez

    If using this code for research or industrial purposes, please cite:
    M.A. Moreno-Mateos, K. Danas, D. Garcia-Gonzalez
    Influence of magnetic boundary conditions on the quantitative modelling of magnetorheological elastomers.
    Mechanics of Materials, 2023.
    
    Reference pape : https://doi.org/10.1016/j.mechmat.2023.104742
    Original code : https://doi.org/10.5281/zenodo.8129310

    Code modified by Jolan Cadieux-Langevin for a 2 phases simulations.
    First phase is the slender/cylinder, Second phase is the air domain.
    Commented out parts include a third phase (permanent magnet)

    i. Inputs/outputs parameters
    ii. CHECK BOUNDARY CONDITIONS! 

    ### Introduction
    This code simulates the behavior of a magnetorheological elastomer (MRE) sample under a magnetic field using the FEniCS finite element library. It employs a vector-potential formulation to model the magnetic field and accounts for finite strains and remanent magnetization. The simulation includes two phases: the MRE sample and an air domain (free space). A third phase for a permanent magnet is available but commented out.

    The simulation proceeds in two stages:
    1. Pre-magnetization ramp: applies a remanent magnetic field to the MRE.
    2. Magnetic actuation: subjects the MRE to an external magnetic field.

    The code solves a coupled problem for the displacement field (mechanical) and the magnetic vector potential using a mixed finite element approach. The mechanical behavior is modeled with a Gent hyperelastic model plus volumetric penalization, while the magnetic behavior incorporates material permeability and remanent magnetization.

    Boundary conditions enforce a homogeneous magnetic field via the vector potential and mechanically fix the air domain. Results are saved to an XDMF file for visualization in Paraview.
"""


import fenics as fn
import numpy as np
from numpy import array
#from dolfin import *
import matplotlib.pyplot as plt
from ufl import indices

# ----------------------------------------------------------------------------------------------
#                   TIME SETTINGS 
# ----------------------------------------------------------------------------------------------
# - `T`: Total simulation time (dimensionless).
# - `T_ramp`: Time when the pre-magnetization ramp ends and actuation begins.
# - `steps`: Tracks the number of completed time steps.
# - `tot_steps`: Total number of time steps for the simulation.
# - `t`: Current simulation time.
# - `tiempo`: FEniCS Expression for time, used in time-dependent parameters.
# - `dt`: Size of each time step.
T = 1
T_ramp=0.5
steps = 0
tot_steps = 20
t = 0.0
tiempo = fn.Expression("t", t=t, degree=0)
dt = T / tot_steps # time step size

# ----------------------------------------------------------------------------------------------
#                    MESH READING  
#  Loads a mesh from an XDMF file containing the MRE sample and air domain.
# ----------------------------------------------------------------------------------------------
# - `foldername`: Directory for output files.
# - `name`: Base name for output XDMF files.
# - `mesh_file`: Input mesh file name.
# - `model`: Path to the mesh file.
# - `importMesh`: Function to read the mesh and subdomain markers.
# - `domains`: MeshFunction assigning subdomain indices (0: MRE, 1: air).
# - `dx`: Integration measure over subdomains with quadrature degree 2.
foldername = 'results'                  
name = "cylinder_with_air_domain_bulky.xdmf"       
mesh_file = 'cylinder_with_air_domain_bulky.xdmf'   
model = f'./model/{mesh_file}'      
print('Processing mesh_file :', mesh_file)

# Read the .xdmf mesh file
def importMesh(file):
    ## Import mesh
    mesh = fn.Mesh()
    with fn.XDMFFile(file) as infile: 
        infile.read(mesh)
    mvc = fn.MeshValueCollection("size_t", mesh, 3)
    with fn.XDMFFile(file) as infile: 
        infile.read(mvc, "phase")
    materials = fn.MeshFunction("size_t",mesh, mvc)
    return mesh,materials,mvc

mesh,materials,mvc = importMesh(model)
# Read the domains of the mesh
domains = fn.MeshFunction("size_t", mesh, dim=3)
domains.set_all(0)
domains.array()[:]=materials.array()[:]-1
dx = fn.Measure('dx')(domain=mesh,subdomain_data=domains)
dx = dx(metadata={'quadrature_degree': 2})
print("Unique domain tags:", set(domains.array()))

# ----------------------------------------------------------------------------------------------
#                   DEFINE FUNCTION SPACES 
#       Defines trial, test, and solution functions for the variational problem.
# ----------------------------------------------------------------------------------------------
# **Function spaces:**
# - `P1`: Vector space for displacements (CG, degree 2).
# - `P2`: Vector space for magnetic vector potential (CG, degree 2).
# - `V`: Mixed function space combining displacement and vector potential.
# - `TT`: Tensor space for stresses (DG, degree 0).
# - `VV`: Vector space for fields like magnetic induction (DG, degree 0).
# - `P0`: Scalar space for subdomain visualization (DG, degree 0).

P1 = fn.VectorFunctionSpace(mesh, 'CG', 2) #Displacements
P2 = fn.VectorFunctionSpace(mesh, 'CG', 2) #Magnetic vector-potential
## Use ufl_element() to return the UFL element of the function spaces
P1elem = P1.ufl_element() # Displacements
P2elem = P2.ufl_element() # Magnetic potential
## Define mixed function space specifying underlying finite element
Vele = fn.MixedElement([P1elem,P2elem])
V = fn.FunctionSpace(mesh, Vele)

#Function spaces to interpolate the results and plot in Paraview.
TT = fn.TensorFunctionSpace(mesh,'DG',0) # Displacement Tensor
VV = fn.VectorFunctionSpace(mesh,'DG',0) # Magnetic vector potential
P0 = fn.FunctionSpace(mesh, "DG", 0) # Required for phases visualisation in Paraview
domains_func = fn.Function(P0)
domains_func.vector()[:] = domains.array()

# Functions in the function space V
du = fn.TrialFunction(V)
w = fn.Function(V)
d = w.geometric_dimension()
v = fn.TestFunction(V)

(u, A) = fn.split(w)
(v_u, v_A) = fn.split(v)
(d_u, d_A) = fn.split(du)

# ----------------------------------------------------------------------------------------------
#                   CONSTITUTIVE PARAMETERS 
# Phases: MRE sample - Free space - Permanent magnet
# If permanent magnet included : 
# mu0Hr_max : if (-1) --> mu0Hr_max is negative
#             if (1) --> mu0Hr_max is positive
# ----------------------------------------------------------------------------------------------
# - Lists are ordered: [MRE sample, air domain] (permanent magnet commented out).
# - `G`: Shear modulus [GPa].
# - `jm`: Gent stretchability parameter.
# - `poisson`: Poisson's ratio.
# - `material_magnet`: Relative magnetic permeability.
# - `mu_0`: Vacuum permeability.
# - `mu0Hr_max`: Maximum remanent magnetic field strength (direction via sign).
# - `Hr_mag`: Time-dependent remanent field expression for the MRE.
# - `v_mag`: Direction of remanent field per phase.

# Physical parameters 
G = [2.9e-3, 80e-6]#, 81e3] # G modulus [GPa]
jm=[1.7, 2000]#,2000] # Gent strechability 
poisson=[0.47,0.2]#,0.49]
material_magnet = [1,1]#,1] # Relative magnetic permitivity

# hMRE Sample
mu_0 = 4e-1*np.pi # = 1.257
mu0Hr_max = -30e-03 * (1) # Remanent Magnetic field hMRE. (1)/(-1) is the direction
Hr_mag = fn.Expression('tiempo <= T_ramp ? mu0Hr_max/mu_0*tiempo/T_ramp : mu0Hr_max/(mu_0)',\
   mu0Hr_max=mu0Hr_max, tiempo=tiempo, T_ramp=T_ramp, mu_0=mu_0, degree=0)
v_mag = [(0,0,1),(0,0,0)]# ,(0,0,0)] # v_mag to apply on each phase

# ----------------------------------------------------------------------------------------------
#                   BOUNDARY CONDITIONS 
# BCs has to match the dimension of the air domain 
# fn.Expression has to be used to access x[0] and x[1]
# For an Homogeneous magnetic field, the Curl of the bc_ref1 (magnetic induction) 
# has to be constant. 
#   -e_1, e_2 --> B in (+) z direction
#   e_1, -e_2 --> B in (-) direction
# ----------------------------------------------------------------------------------------------
# **Boundary conditions:**
# - `domain_air`: Identifies air domain boundaries at x, y, z = ±15.
# - `domain_permanentmag`: Defines permanent magnet boundary (unused).
# - `bc_ref1`: Sets vector potential on air boundary for a uniform magnetic field in z-direction.
# - `bc_air`: Fixes displacement on air boundary.
# - `bc_magnets`: Fixes displacement on permanent magnet (unused).

def domain_air(x,on_boundary): # has to match the air domain dimension
    return fn.near(x[0],-15) or fn.near(x[0],15) or fn.near(x[1],-15) or fn.near(x[1],15)\
        or fn.near(x[2],-15) or fn.near(x[2],15)
def domain_permanentmag(x,on_boundary):
    return fn.between(x[0]**2+x[1]**2, (0,10**2)) and fn.between(x[1],(-15,-5))

# Homogenoues magnetic field
bc_ref1 = fn.DirichletBC(V.sub(1),
                         fn.Expression(("0.5*0.3*x[1]",
                                        "-0.5*0.3*x[0]",
                                        "0 + eps"),
                                       degree=1,
                                       mu0Hr_max=mu0Hr_max,
                                       mu_0=mu_0,
                                       eps=fn.DOLFIN_EPS),
                         domain_air, method="pointwise")
#"0.5 * mu0Hr_max / mu_0 * x[1]",
                                        #"-0.5 * mu0Hr_max / mu_0 * x[0]",

# Boundary air mechanically fixed (if air domain included)
bc_air=fn.DirichletBC(V.sub(0),
                      fn.Constant((0.0+fn.DOLFIN_EPS,0.0+fn.DOLFIN_EPS,0.0+fn.DOLFIN_EPS)),
                      domain_air,method="pointwise")
# Permanent magnet fixed (unused)
bc_magnets=fn.DirichletBC(V.sub(0),fn.Constant((0.0+fn.DOLFIN_EPS,0.0+fn.DOLFIN_EPS,0.0+fn.DOLFIN_EPS)),
                          domain_permanentmag,method="pointwise")

bc1 = [bc_air, bc_ref1] 
for bc in bc1:
    print(f"BC DOFs: {bc.get_boundary_values().keys()}")
                                
# ----------------------------------------------------------------------------------------------
#                   WEAK/STRONG FORMULATION AND FUNCTION 
# VECTOR AND TENSOR FIELDS DERIVED FROM VECTORIAL POTENTIAL FIELDS A AND u (magnetic potential and displacement fields).
# ----------------------------------------------------------------------------------------------
# **Weak formulation:**
# - Defines key tensors and fields:
#   - `F`: Deformation gradient.
#   - `CG`, `BG`: Cauchy-Green tensors.
#   - `B`: Magnetic induction (curl of A).
# - `H`: Magnetic field function per phase, includes Maxwell, susceptibility, and remanence terms.
# - `Pmagr`: Magnetic stress contribution from remanent field.
# - `P`: Total first Piola-Kirchhoff stress (mechanical + magnetic).
# - `Res`: Weak form residual combining mechanical, magnetic, and gauge terms.
## 
d = len(u)                      # Spatial dimension
I = fn.Identity(d)                 # Identity tensor
F = I + fn.grad(u)                 # Deformation gradient from current time step
CG = fn.dot(F.T,F)                 # Right Cauchy-Green tensor
BG = fn.dot(F,F.T)                 # Left Cauchy-Green tensor
B = fn.curl(A)                     # Magnetic induction from the vector-potential field

def H(F,B,i):
    mu = material_magnet[i]
    chi = mu - 1
    I5 = fn.dot(fn.dot(F, B), fn.dot(F, B))
    I5root=fn.sqrt(I5)
    dI5dB = 2 * fn.dot(CG, B)
    Hmaxw= 1 / (2 * mu_0 * fn.det(F)) * dI5dB
    x = fn.SpatialCoordinate(mesh)
    Hr=Hr_mag*fn.as_vector(v_mag[i])
    #Hr_1=Hr_permag*fn.as_vector(v_permag[i])
    return -chi / (2 * mu_0 * (1 + chi) * fn.det(F) ** 2) *dI5dB \
           + Hmaxw + fn.det(F)*fn.dot(F,Hr) # + fn.det(F)*fn.dot(F,Hr_1)

def Pmagr(F,i,B):
    h,j,k,l=indices(4);
    mu = material_magnet[i]
    x = fn.SpatialCoordinate(mesh)
    Hr=Hr_mag*fn.as_vector(v_mag[i])
    H_=H(F,B,i)
    HdiadHr=fn.as_tensor(H_[h]*Hr[j],(h,j))
    I1U=fn.tr(F)
    I2U=0.5*(fn.tr(F)**2+fn.tr(fn.dot(F,F)))
    I3U=fn.det(F)
    delta=8*(I1U*I2U-I3U)*I3U
    I=fn.Identity(3)
    dI5erdC=-4/delta * (I1U*fn.dot(fn.dot(F,fn.sym(HdiadHr)),F) - I1U**2 * (fn.dot(F,fn.sym(HdiadHr)) + fn.dot(fn.sym(HdiadHr),F)) + \
            (I1U*I2U-I3U) * (fn.dot(fn.dot(F,fn.sym(HdiadHr)),fn.inv(F)) + fn.dot(fn.dot(fn.inv(F),fn.sym(HdiadHr)),F)) + (I1U**3 + I3U) * fn.sym(HdiadHr) - \
            I1U**2*I2U* (fn.dot(fn.sym(HdiadHr),fn.inv(F)) + fn.dot(fn.inv(F),fn.sym(HdiadHr))) + (I1U**2*I3U + (I1U*I2U-I3U)*I2U)*fn.dot(fn.dot(fn.inv(F),fn.sym(HdiadHr)),fn.inv(F)))
    dCdF=fn.as_tensor(F[k,j]*I[h,l]+F[k,h]*I[j,l],(h,j,k,l))
    dI5erdCdotdCdF=fn.as_tensor(dI5erdC[k,l]*dCdF[k,l,h,j],(h,j))
    Pmagr_= dI5erdCdotdCdF
    return fn.as_tensor(B[j]*Hr[k],(j,k))
    
def P(F, i, B):
    shear = G[i];
    mu = material_magnet[i]
    Finv=fn.inv(F);FinvT=Finv.T
    chi = mu - 1
    I5 = fn.dot(fn.dot(F, B), fn.dot(F, B))
    h, j, k, l = indices(4)
    dCdF = fn.as_tensor(F[k, j] * I[h, l] + F[k, h] * I[j, l], (h, j, k, l))
    BdCdF = fn.as_tensor(B[h] * dCdF[h, j, k, l], (j, k, l))
    dI5dF = fn.as_tensor(BdCdF[j, k, l] * B[l], (j, k))
    beta=2*poisson[i]/(1-2*poisson[i])
    CauchyMech=shear*jm[i]/(jm[i]-(CG[0,0]+CG[1,1]+CG[2,2])+3)*BG
    PmechGent=fn.det(F)*CauchyMech*fn.inv(F.T)
    Pvol= shear*(-fn.det(F)**(-beta)*fn.inv(F).T)
    return PmechGent+Pvol\
                  -chi/(2*mu_0*(1+chi)*fn.det(F)**2)*dI5dF\
                  +Pmagr(F,i,B)

# VARIATIONAL PROBLEM - WEAK FORM
Res_u =\
  fn.inner(P(F,0,B),fn.grad(v_u))*dx(0)\
  + fn.inner(P(F,1,B),fn.grad(v_u))*dx(1) 
  #+ fn.inner(P(F,2,B),fn.grad(v_u))*dx(2) # perma. magnet
Res_A =\
  fn.dot(H(F,B,0),fn.curl(v_A))*dx(0)\
  + fn.dot(H(F,B,1),fn.curl(v_A))*dx(1)
  #+ fn.dot(H(F,B,2),fn.curl(v_A))*dx(2) # perma. magnet

Gauge_constraint = (fn.dot(fn.div(A),fn.div(v_A))/mu_0/1e-7/1e0) * dx

Res = Res_A + Res_u + Gauge_constraint

# Compute directional derivative about w in the direction of du (Jacobian)
Jacobian = fn.derivative(Res, w, du)

# ----------------------------------------------------------------------------------------------
#                   Lagrangian/Eulerian fields (project them accounting for the subdomains 
# ----------------------------------------------------------------------------------------------
# - Projects fields onto appropriate spaces for visualization in Paraview.
# - Handles each subdomain separately.
def split_projectP(P,F,B,VSpace):
    u = fn.TrialFunction(VSpace)
    v = fn.TestFunction(VSpace)
    f_tot = fn.Function(VSpace)
    fn.solve(fn.inner(u,v)*dx == fn.inner(P(F,0,B),v)*dx(0)\
        + fn.inner(P(F,1,B),v)*dx(1), f_tot)
        #+ fn.inner(P(F,2,B),v)*dx(2),f_tot) # perma. magnet
    return f_tot
def split_projectCauchy(P,F,B,VSpace):
    u = fn.TrialFunction(VSpace)
    v = fn.TestFunction(VSpace)
    f_tot = fn.Function(VSpace)
    fn.solve(fn.inner(u,v)*dx ==\
              (1/fn.det(F))*fn.inner(fn.dot(P(F,0,B),F.T),v)*dx(0)\
              + (1/fn.det(F))*fn.inner(fn.dot(P(F,1,B),F.T),v)*dx(1), f_tot)
              #+ (1/fn.det(F))*fn.inner(fn.dot(P(F,2,B),F.T),v)*dx(2), f_tot) # perma. magnet
    return f_tot
def split_projectH(H,F,B,VSpace):
    u = fn.TrialFunction(VSpace)
    v = fn.TestFunction(VSpace)
    b_tot = fn.Function(VSpace)
    fn.solve(fn.inner(u,v)*dx == fn.inner(H(F,B,0),v)*dx(0)\
        + fn.inner(H(F,B,1),v)*dx(1), b_tot)
        #+ fn.inner(H(F,B,2),v)*dx(2),b_tot) # perma. magnet
    return b_tot
def split_projecth(H,F,B,VSpace):
    u = fn.TrialFunction(VSpace)
    v = fn.TestFunction(VSpace)
    b_tot = fn.Function(VSpace)
    fn.solve(fn.inner(u,v)*dx == fn.inner(fn.dot(fn.inv(F.T),H(F,B,0)),v)*dx(0)\
        + fn.inner(fn.dot(fn.inv(F.T),H(F,B,1)),v)*dx(1), b_tot)
        #+ fn.inner(fn.dot(fn.inv(F.T),H(F,B,2)),v)*dx(2),b_tot) # perma. magnet
    return b_tot
def split_projectM(B,FF_sol,H,V3vm):
        u = fn.TrialFunction(V3vm)
        v = fn.TestFunction(V3vm)
        m_tot = fn.Function(V3vm)
        fn.solve(fn.inner(u,v)*dx == fn.inner(fn.dot(fn.dot(FF_sol.T,FF_sol),B)/(mu_0)-H(F,B,0),v)*dx(0)\
          + fn.inner(fn.dot(fn.dot(FF_sol.T,FF_sol),B)/(mu_0)-H(F,B,1),v)*dx(1), m_tot)
          #+ fn.inner(fn.dot(fn.dot(FF_sol.T,FF_sol),B)/(mu_0)-H(F,B,2),v)*dx(2),m_tot) # perma. magnet
        return m_tot
def split_projectm(B,FF_sol,H,V3vm):
        u = fn.TrialFunction(V3vm)
        v = fn.TestFunction(V3vm)
        m_tot = fn.Function(V3vm)
        fn.solve(fn.inner(u,v)*dx == fn.inner(fn.dot(FF_sol,B)/(mu_0)-fn.dot(fn.inv(FF_sol.T),H(F,B,0)),v)*dx(0)\
          + fn.inner(fn.dot(FF_sol,B)/(mu_0)-fn.dot(fn.inv(FF_sol.T),H(F,B,1)),v)*dx(1), m_tot)
          #+ fn.inner(fn.dot(FF_sol,B)/(mu_0)-fn.dot(fn.inv(FF_sol.T),H(F,B,2)),v)*dx(2),m_tot) # perma. magnet
        return m_tot
def B_eul(B,F):
    return 1/fn.det(F)*fn.dot(F,B)


#To integrate fields on the upper surface of the MRE sample (reaction force on plate rheometer).
#boundary_markers = fn.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
#boundary_markers.set_all(0)

# ----------------------------------------------------------------------------------------------
#                           SOLVING PARAMETERS 
# PETSc SNES solver parameters --> Setup Non-linear variational problem
# ----------------------------------------------------------------------------------------------
# **Solver setup:**
# - Configures PETSc SNES solver for the nonlinear problem.
# - Optimizes form compilation for performance.

snes_solver_parameters = {"nonlinear_solver": "snes",
                          "symmetric": True,
                          "snes_solver": {"maximum_iterations": 5,
                                          "report": True,
                                          "line_search": "bt",
                                          "linear_solver": "mumps",
                                          "method": "newtonls",
                                          "absolute_tolerance": 1e-4,
                                          "relative_tolerance": 1e-4,
                                          "error_on_nonconvergence": False}}

# Form compiler options
fn.parameters["form_compiler"]["optimize"] = True
fn.parameters["form_compiler"]["cpp_optimize"] = True

problem1 = fn.NonlinearVariationalProblem(Res, w, bc1, J=Jacobian)
solver_problem1 = fn.NonlinearVariationalSolver(problem1)
solver_problem1.parameters.update(snes_solver_parameters)

# Save results into an .xdmf
file_results = fn.XDMFFile(name)
file_results.parameters["flush_output"] = True
file_results.parameters["functions_share_mesh"] = True
file_results.parameters["rewrite_function_mesh"] = False

# ----------------------------------------------------------------------------------------------
#                                       SOLVING LOOP 
# ----------------------------------------------------------------------------------------------
# - Iterates over time steps, solving the problem and projecting results.
# - Saves fields to XDMF for visualization.

while (t < T):
    # Define boundary conditions
    steps += 1
    # Update time parameters
    tiempo.t = t
    t += dt
    # Solve weak form
    solver_problem1.solve()
    (u, A) = w.split()

    PTensor = split_projectP(P,F,B,TT)
    CauchyTensor = split_projectCauchy(P,F,B,TT)
    HVector=split_projectH(H,F,B,VV)
    hVector=split_projecth(H,F,B,VV)
    BVector=fn.project(B,VV)
    bVector=fn.project(B_eul(B,F),VV)
    MVector=split_projectM(B,F,H,VV)
    mVector=split_projectm(B,F,H,VV)

    #Rename fields
    u.rename("Displacement", "u")
    A.rename("Field", "A") # magnetic vector potential H. Applying Curl will yield induction B
    PTensor.rename("Nominal (P) Stress", "PTensor")
    CauchyTensor.rename("Cauchy (Sigma) Stress", "CauchyTensor")
    BVector.rename("Magnetic B-field", "BVector")
    HVector.rename("Magnetic H-field", "HVector")
    MVector.rename("Magnetic M-field", "MVector")
    bVector.rename("Magnetic b-field", "bVector")
    hVector.rename("Magnetic h-field", "hVector")
    mVector.rename("Magnetic m-field", "mVector")

    print("Domain values before writing:", np.unique(domains.array()))
    domains_func.rename("Phases", "domains")
    
    # Write to .xdmf results file
    file_results.write(domains_func,t)
    file_results.write(u,t)
    file_results.write(A,t)
    file_results.write(PTensor,t)
    file_results.write(CauchyTensor,t)
    file_results.write(BVector,t)
    file_results.write(bVector,t)
    file_results.write(MVector,t)
    file_results.write(mVector,t)
    file_results.write(HVector,t)
    file_results.write(hVector,t)

    #Printouts
    print('Step: ',steps,'/',tot_steps)
    print('dt: ',dt)
    print('Time: ',t,'/',T)

file_results.close()
print("Finished")
