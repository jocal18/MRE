'''
 __    __    ______   .___  ___.   ______     _______  _______ .__   __.  __ 
|  |  |  |  /  __  \  |   \/   |  /  __  \   /  _____||   ____||  \ |  | |  |
|  |__|  | |  |  |  | |  \  /  | |  |  |  | |  |  __  |  |__   |   \|  | |  |  ___________
|   __   | |  |  |  | |  |\/|  | |  |  |  | |  | |_ | |   __|  |  . `  | |  | |___________|
|  |  |  | |  `--'  | |  |  |  | |  `--'  | |  |__| | |  |____ |  |\   | |  | 
|__|  |__|  \______/  |__|  |__|  \______/   \______| |_______||__| \__| |__|
                                                                             
 ________      ___   .___________. __    ______   .__   __.                  
|       /     /   \  |           ||  |  /  __  \  |  \ |  |                  
`---/  /     /  ^  \ `---|  |----`|  | |  |  |  | |   \|  |                  
   /  /     /  /_\  \    |  |     |  | |  |  |  | |  . `  |                  
  /  /----./  _____  \   |  |     |  | |  `--'  | |  |\   |                  
 /________/__/     \__\  |__|     |__|  \______/  |__| \__|                  
                                                                             
  ______   ______    _______   _______                                       
 /      | /  __  \  |       \ |   ____|                                      
|  ,----'|  |  |  | |  .--.  ||  |__                                         
|  |     |  |  |  | |  |  |  ||   __|                                        
|  `----.|  `--'  | |  '--'  ||  |____                                       
 \______| \______/  |_______/ |_______|                                                                                           
 
 '''


"""
FEniCS complete code for magneto-visco-mechanical homogenization: 3D full-field magneto-mechanical 
homogenization that accounts for: finite strains, viscoelasticity, incompressibility and mixed 
stress/magnetic flux control.

Copyright (C) 2021: Sergio Lucarini, Miguel Angel Moreno-Mateos, Kostas Danas, Daniel Garcia-Gonzalez

If using this code for research or industrial purposes, please cite:
S. Lucarini, M.A. Moreno-Mateos, K. Danas, D. Garcia-Gonzalez
Insights into the viscohyperelastic response of soft magnetorheological
elastomers: competition of macrostructural versus microstructural players
International Journal of Solids and Structures, 2022.
Supplementary note - Jolan Cadieux-Langevin 

Ref paper : https://www.sciencedirect.com/science/article/pii/S0020768322004346 
DOI : https://doi.org/10.1016/j.ijsolstr.2022.111981 

Coode modified by Jolan Cadieuxc-Langevin.  
The full original code can be found at : https://zenodo.org/records/7112767 

Overview:
This code implements a 2-phase homogenization framework using FEniCS to simulate the mechanical and magnetic behavior 
of a composite material, such as a magnetorheological elastomer with inclusions. It supports periodic boundary conditions 
and can handle strain-controlled or stress-controlled mechanical loading, with optional magnetoactive, incompressible, 
and viscous material behaviors.

The simulation advances through time steps defined in a timer array, applying mechanical and/or magnetic loading. 
It solves a coupled problem involving displacement, magnetic potential (if magnetoactive), pressure (if incompressible), 
and internal variables (if viscous) using a Newton-Raphson solver with adaptive time stepping. Results are saved to an 
XDMF file for visualization and post-processing.
"""

import numpy as np
import fenics as fn
from fenics import XDMFFile
from aux_homog_nl import importMesh, PeriodicBoundary, post, adjust_xdmf
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------------
#                   OUTPUT AND MESH SETTINGS
# ----------------------------------------------------------------------------------------------
# Define output directory and file names
foldername = 'results'                   # Directory where output files will be saved
output_filename = 'simple_compression_BCC_v4'     # Base name for output XDMF and text files
mesh_file = 'BCC_VF03_v4.xdmf'    # Name of the input mesh file in XDMF format
model = f'./model/homo/{mesh_file}'            # Full path to the mesh file
print('Processing mesh_file :', mesh_file)  # Inform user of the mesh being processed

# ----------------------------------------------------------------------------------------------
#                   TIME SETTINGS
# - `timer`: A list where each row defines a simulation step with [initial time increment, final time, min time increment, max time increment].
# - This allows for multiple loading steps, similar to Abaqus, with adaptive time stepping for convergence.
# - Here, we define a single step for mechanical compression.
# ----------------------------------------------------------------------------------------------

#b_end = 15e-03;  b_rate= 1e-03 # Magnetic rate, Final magnetic field
f_end = 0.1; strain_rate = 0.1 # Final mechanical ramp, Strain rate
timer = []
#t_ramp = b_end/b_rate  #0.3 # also t_end -> t_end=timer[1]
t_ramp = f_end / strain_rate # For compression example
# set steps 
timer.append([t_ramp / 100, t_ramp, t_ramp/10000, 2*t_ramp/100]) # [initial dtime, total time end, dtime min, dtime max]
#timer.append([t_ramp / 100, 2*t_ramp + 100 * t_ramp, t_ramp / 10000, 2 * t_ramp/100]) # Hold timer to extend past the fist timer
timer = np.array(timer)

# ------------------------------------------------------------------------
#              Material Properties
# [(domain), (inclusions)]
# elastic :(shear modulus, bulk modulus = 0, shear modulus of the viscous part, viscosity parameter tau)
# magnetic: (relative magnetic permeability, magnetic saturation)
# ------------------------------------------------------------------------
material_elas = [(1.034e-3, 0, 1.515e-3, 0.217e-3), (81e03, 0, 0, 1)]
material_magnet = [(1, 1), (30, 2.5)]

# ------------------------------------------------------------------------
#              Magnetic Field, Incompressibility, Viscosity
# 
# if magnetoactive -> Vv=fn.FiniteElement('CG', mesh.ufl_cell(), 2)
# if incompressibility -> Vpr=fn.FiniteElement('DG', mesh.ufl_cell(), 0)
# if viscosity -> VFV=fn.VectorElement('DG', mesh.ufl_cell(), 1, dim=6)
# if not viscosity -> FV_sol_t=fn.Constant((0,0,0,0,0,0)) (state variable)
# ------------------------------------------------------------------------

magnetoactive = False
incompressibility = True 
viscosity = False   

# ------------------------------------------------------------------------
#                           SOLVER PARAMETERS
# ------------------------------------------------------------------------
# Tolerances and optimization
maxiter_nw = 10 # Max iterations before reducing time step
tol_nw = 5e-5    #  tolerance for Newton convergence
tol_nw_forces = 1e-3  # Tolerance for force residuals
solver = fn.LUSolver('mumps') 
fn.parameters["form_compiler"]["optimize"] = False # Disable optimization for stability
fn.parameters["form_compiler"]["cpp_optimize"] = False # Disable C++ compilation

# ------------------------------------------------------------------------
#                           MACRO LOADING  
#       0 -> Strain-Controlled / 1 -> Stress-Controlled
#       control[i,j] : (i)-direction along the (j)-normal in the reference configuration.
#       colum notation: [11, 12, 13], [21, 22, 23], [31, 32, 33] (direction, normal)
# - `control`: 3x3 matrix where 0 = strain-controlled, 1 = stress-controlled for each deformation gradient component.
#  - Rows: direction (1, 2, 3); Columns: normal (1, 2, 3) in reference configuration.
# - `control_mag_macro`: Vector where 0 = H-controlled, 1 = B-controlled for magnetic field components.
# - `F_macro`, `P_macro`, `H_macro`, `B_macro`: Expressions defining time-dependent macroscopic loading.
# ------------------------------------------------------------------------

#control = [[1, 0, 0], [1, 1, 0], [1, 1, 0]] # Compression control example
#control = [[1, 0, 0,], [1, 1, 0], [1, 1, 1]] # Magnetorestriction example
control = [[1, 1, 0], [1, 1, 0], [1, 1, 0]] 

# control B or H -> 0-H / 1-B
control_mag_macro=np.array([0,0,0])

F_macro = fn.Expression((("0","0","0"),
                        ("0","0","0"),
                        #("0","0","0")), # Magnetostriction example
                        ("0","0","t<=t_ramp ? "+str(-f_end/t_ramp)+"*t : "+str(-f_end))), # Compression example
                        degree=0, t=0, t_ramp=t_ramp)

P_macro=fn.Expression((("0","0","0"),
                        ("0","0","0"),
                        ("0","0","0")), 
                        degree=0, t=0)

H_macro=fn.Expression((("0.0001","0.0001","0.0001")), degree=0, t = 0)  # keep small value to avoid singularity           
#B_macro=fn.Expression((("0.00","0.00","t<=t_ramp ? "+str(b_end/t_ramp) +  
#             "*t : "+str(b_end))), degree=0, t = 0, t_ramp=t_ramp) # B in z direction
#B_macro=fn.Expression((("0.00","t<=t_ramp ? "+str(b_end/t_ramp)+ 
#             "*t : "+str(b_end), "0.00")), degree=0, t = 0, t_ramp=t_ramp) # B in y direction
B_macro=fn.Expression(("0.00","0.00", "0.00"), degree=0, t = 0, t_ramp=t_ramp) # Compression example. B_macro null 

if not magnetoactive: 
    control_mag_macro=np.array([0, 0, 0])

control_mech_macro=np.array(control)
control_mech_index=np.zeros([3,3])
control_mech_index[np.where(control_mech_macro == 1)[0], 
                   np.where(control_mech_macro == 1)[1]] = np.arange(np.sum(control_mech_macro))

control_mag_index=np.zeros([3]);
control_mag_index[np.where(control_mag_macro==1)] = np.arange(np.sum(control_mag_macro))

# ------------------------------------------------------------------------
#                   Mesh reading and preparation
# - Load the mesh and subdomain markers from the XDMF file using a custom function.
# - `domains`: Assigns subdomain indices (0 = matrix, 1 = inclusions).
# - `vertices`: Defines the unit cell vertices for periodic boundary conditions.
# - `dx`: Integration measure over subdomains with a quadrature degree of 2 for accuracy.
# ------------------------------------------------------------------------
mesh, materials, mvc = importMesh(model)
domains = fn.MeshFunction("size_t", mesh, dim=3)
domains.set_all(0)
domains.array()[:]=materials.array()[:]-1
vertices = np.array([[0, 0.,0],[1, 0.,0],[1, 1.,0],[0, 1.,0],\
                     [0, 0.,1],[1, 0.,1],[1, 1.,1],[0, 1.,1]])
dx = fn.Measure('dx')(domain=mesh, subdomain_data=domains)
dx = dx(metadata={'quadrature_degree': 2}) 
# initialize time increment
DT_macro = fn.Expression("dt", dt=0, degree=0)

# ----------------------------------------------------------------------------------------------
#                       CONSTITUTIVE FUNCTIONS
# derivatives  of the energy functional Psi, constitutive relations
# dPsi/dF - Equation 15. Returns Total first Piola-Kirchhoff stress. 
# Sum of Mechanical, Maxwell ang Magnetism energy. 
# - `P`: Total first Piola-Kirchhoff stress (mechanical + magnetic contributions).
# - `P0`: Mechanical part of the stress (excludes magnetic effects).
# - `B`: Magnetic induction derived from the magnetic field H.
# - `compare`: Viscous flow rule for internal variables (if viscosity is enabled).
# ----------------------------------------------------------------------------------------------


def P(F,p, i, H, FV):
    """Compute total first Piola-Kirchhoff stress for subdomain i."""
    shear,bulk,shearv,visco = material_elas[i]
    mu,ms = material_magnet[i]
    chi=mu-1
    Finv=fn.inv(F);FinvT=Finv.T
    J=fn.det(F)
    FE = fn.dot(F,fn.inv(FV))
    I5=fn.dot(fn.dot(H,Finv),fn.dot(H,Finv))
    I5r=fn.sqrt(I5)
    dI5dF=-2*fn.outer(fn.dot(FinvT,H),fn.dot(fn.dot(H,Finv),FinvT))
    return shear*F-p*J*FinvT+\
      fn.dot(shearv*FE,fn.inv(FV))-\
      J*dI5dF*4e-1*np.pi/2-\
      J*I5*4e-1*np.pi/2*FinvT-\
      J*4e-1*np.pi*ms/2/(I5r+fn.DOLFIN_EPS)*(-1+2/(fn.exp(-2*chi/(ms+fn.DOLFIN_EPS)*I5r)+1))*dI5dF-\
      J*4e-1*np.pi*FinvT*(ms**2/(chi+fn.DOLFIN_EPS)*fn.ln(0.5*(fn.exp(chi/(ms+fn.DOLFIN_EPS)*I5r)+fn.exp(-chi/(ms+fn.DOLFIN_EPS)*I5r))))

# dPsi/dF_mech - Equation 12
def P0(F, p, i, FV):
    """Compute mechanical first Piola-Kirchhoff stress for subdomain i."""
    shear, bulk, shearv, visco = material_elas[i]
    Finv=fn.inv(F)
    FinvT=Finv.T
    J=fn.det(F)
    FE = fn.dot(F,fn.inv(FV))
    return shear*F-p*J*FinvT+fn.dot(shearv*FE,fn.inv(FV))

# Magnetic constitutive equation
# dPsi/dH  - Equation 16
def B(H, i, F):
    """Compute magnetic induction B for subdomain i."""
    mu,ms = material_magnet[i]
    chi=mu-1
    Finv=fn.inv(F);FinvT=Finv.T
    J=fn.det(F)
    I5r=fn.sqrt(fn.dot(fn.dot(H,Finv),fn.dot(H,Finv)))
    dI5dH=2*fn.dot(fn.dot(FinvT,H),FinvT)
    return J*4e-1*np.pi*(1/2+ms/2/(I5r+fn.DOLFIN_EPS)*(-1+2/(fn.exp(-2*chi/(ms+fn.DOLFIN_EPS)*I5r)+1)))*dI5dH

# Viscous flow rule - Section 3.3.3
# Equation 22
# describe the evolution of the internal variable F_viscous
def compare(F, i, FV,DT):
    """Compute viscous flow for internal variables in subdomain i."""
    shear, bulk, shearv,visco = material_elas[i]
    FE = fn.dot(F,fn.inv(FV))
    J=fn.det(F)
    # Viscous flow
    # Pnm1 = fn.dot(shearv*FE-shearv*fn.inv(FE).T,fn.inv(FV))
    Pnm1 = fn.dot(shearv*FE,fn.inv(FV))
    sigman = fn.dot(Pnm1,F.T)/J
    sigmaDev = sigman - 1./3.*fn.tr(sigman)*fn.Identity(3)
    DV = 1/2**0.5/visco*sigmaDev
    CompareFV = DT*fn.dot(fn.dot(fn.inv(FE),DV),F) 
    return fn.as_vector((CompareFV[0,0], CompareFV[1,1], CompareFV[2,2], CompareFV[0,1]+\
     CompareFV[1,0],CompareFV[1,2]+CompareFV[2,1], CompareFV[0,2]+CompareFV[2,0])) 

# ------------------------------------------------------------------------
#                  Defining Functional Space
# fn.VectorElementdefines a vector finite element within the Unified Form Language (UFL) framework. 
# This is used to specify the type of finite element space for a vector-valued field, 
# such as displacement in a mechanics problem. 
# FiniteElement creates a scalar-valued field (i.e., a single-component field). 
# This is often used for quantities like pressure, temperature, or a scalar potential (e.g., magnetic potential).
#
# - Define a mixed function space `V` combining all fields based on simulation flags.
# - `Vu`: Displacement (vector, CG degree 2).
# - `Vv`: Magnetic potential (scalar, CG degree 2, if magnetoactive).
# - `Vpr`: Pressure (scalar, DG degree 0, if incompressible).
# - `VFV`: Viscous variables (vector, DG degree 1, dim=6, if viscous).
# - `Vlambda`, `Vlambdav`: Lagrange multipliers for mechanical and magnetic control.
# ------------------------------------------------------------------------
all_elements=[]

# displacement finite element
Vu = fn.VectorElement('CG', mesh.ufl_cell(), 2, dim=3)
all_elements.append(Vu)

# magnetic potential finite element
if magnetoactive: 
    Vv=fn.FiniteElement('CG', mesh.ufl_cell(), 2)
    all_elements.append(Vv)

# presure potential finite element
if incompressibility: 
    Vpr=fn.FiniteElement('DG', mesh.ufl_cell(), 0)
    all_elements.append(Vpr)

# Presure potential finite element
if viscosity: 
    VFV=fn.VectorElement('DG', mesh.ufl_cell(), 1, dim=6)
    all_elements.append(VFV)

# lagrange multiplyer constant element
if np.sum(control_mech_macro)>0: 
    Vlambda = fn.VectorElement('R', mesh.ufl_cell(), 0, int(np.sum(control_mech_macro)))
    all_elements.append(Vlambda) 

# lagrange multiplier viscous constant element
if magnetoactive and np.sum(control_mag_macro)>0: 
    Vlambdav=fn.VectorElement('R', mesh.ufl_cell(), 0, int(np.sum(control_mag_macro)))
    all_elements.append(Vlambdav)
 
# mixed formulation with periodic pbcs
V = fn.FunctionSpace(mesh, fn.MixedElement(all_elements), constrained_domain=PeriodicBoundary(vertices))

# ----------------------------------------------------------------------------------------------
#                   TEST, TRIAL, AND SOLUTION FUNCTIONS
# ----------------------------------------------------------------------------------------------
"""
Test, Trial, and Solution Functions:
- Define test functions (`W_`), trial functions (`dW`), and solution functions (`W_sol`) in the mixed space.
- Split the solution into individual fields (displacement, pressure, etc.) based on simulation flags.
"""
W_ = fn.TestFunction(V)  # Test functions for the weak form
all_test = fn.split(W_)
dW = fn.TrialFunction(V)  # Trial functions for the Jacobian
W_sol = fn.Function(V)  # Current solution
all_func = fn.split(W_sol)  # Split solution into components
W_sol_t = fn.Function(V)  # Solution at previous time step
deltaW_sol = fn.Function(V)  # Increment in solution
diffW = fn.Function(V)  # Difference between current and previous solution
diffW2 = fn.Function(V)  # Additional buffer for solution difference

iiii = 0  # Index to track subspaces
u_ = all_test[iiii]  # Test function for displacement
u_sol = all_func[iiii]  # Solution for displacement

if magnetoactive:
    iiii += 1
    v_ = all_test[iiii]  # Test function for magnetic potential
    v_sol = all_func[iiii]  # Solution for magnetic potential
else:
    v_sol = fn.Constant(0)  # Zero if not magnetoactive

if incompressibility:
    iiii += 1
    pr_ = all_test[iiii]  # Test function for pressure
    pr_sol = all_func[iiii]  # Solution for pressure
else:
    pr_sol = fn.Constant(0)  # Zero if not incompressible

if viscosity:
    iiii += 1
    FV_ = all_test[iiii]  # Test function for viscous variables
    FV_sol = all_func[iiii]  # Solution for viscous variables
    ifv = 1 * iiii  # Index for viscous variables
else:
    FV_sol = fn.Constant((0, 0, 0, 0, 0, 0))  # Zero vector if not viscous

if np.sum(control_mech_macro) > 0:
    iiii += 1
    lambda_ = all_test[iiii]  # Test function for mechanical Lagrange multipliers
    lambdas = all_func[iiii]  # Solution for mechanical Lagrange multipliers
if magnetoactive and np.sum(control_mag_macro) > 0:
    iiii += 1
    lambdav_ = all_test[iiii]  # Test function for magnetic Lagrange multipliers
    lambdavs = all_func[iiii]  # Solution for magnetic Lagrange multipliers

# Store DOF indices for each subspace
dofindex = [V.sub(ai).dofmap().dofs() for ai in range(iiii + 1)]
# Previous time step solution
all_func_t = fn.split(W_sol_t)
if viscosity:
    FV_sol_t = all_func_t[ifv]  # Previous viscous variables
else:
    FV_sol_t = fn.Constant((0, 0, 0, 0, 0, 0))

# Construct Lagrange multiplier tensors and vectors
lambda_tensor = fn.as_tensor([[lambdas[int(control_mech_index[i, j])] if control_mech_macro[i, j] == 1 else 0 for i in range(3)] for j in range(3)])
lambda__tensor = fn.as_tensor([[lambda_[int(control_mech_index[i, j])] if control_mech_macro[i, j] == 1 else 0 for i in range(3)] for j in range(3)])
lambdav_vector = fn.as_vector([lambdavs[int(control_mag_index[i])] if control_mag_macro[i] == 1 else 0 for i in range(3)])
lambdav__vector = fn.as_vector([lambdav_[int(control_mag_index[i])] if control_mag_macro[i] == 1 else 0 for i in range(3)])

# Construct solution fields
FVV_sol = fn.as_tensor(((FV_sol[0], FV_sol[3] / 2, FV_sol[5] / 2),
                        (FV_sol[3] / 2, FV_sol[1], FV_sol[4] / 2),
                        (FV_sol[5] / 2, FV_sol[4] / 2, FV_sol[2]))) + fn.Identity(3)  # Viscous deformation gradient
FF_sol = F_macro + lambda_tensor + fn.Identity(3) + fn.grad(u_sol)  # Total deformation gradient
H_macro2 = H_macro  # Copy of macroscopic H
B_macro2 = B_macro  # Copy of macroscopic B

if magnetoactive:
    FF_macro = F_macro + lambda_tensor + fn.Identity(3)  # Macroscopic deformation gradient
    HH_sol = H_macro2 + lambdav_vector - fn.grad(v_sol)  # Total magnetic field
    B_macro2 = fn.det(FF_macro) * fn.dot(fn.inv(FF_macro), B_macro)  # Adjust B_macro for consistency
else:
    HH_sol = fn.Constant((0, 0, 0))  # Zero magnetic field

# ----------------------------------------------------------------------------------------------
#                                   WEAK FORM
# - Assemble the residual `form` by summing contributions from:
#  - Mechanical energy (stress equilibrium).
#  - Magnetic energy (if magnetoactive).
#  - Incompressibility constraint (if enabled).
#  - Viscous evolution (if enabled).
#  - Lagrange multipliers for stress/magnetic control.
# - `fint`: List of individual residual components for convergence checking.
# ----------------------------------------------------------------------------------------------
form =  0    
fint=[]

#mechanical energy including F macro and F due to P macro
if magnetoactive: 
    Ptrial_0 = P(FF_sol, pr_sol, 0, HH_sol, FVV_sol)
    Ptrial_1 = P(FF_sol, pr_sol, 1, HH_sol, FVV_sol)
else:
    Ptrial_0 = P0(FF_sol, pr_sol, 0, FVV_sol)
    Ptrial_1 = P0(FF_sol, pr_sol, 1, FVV_sol) 
Ptrial0_0 = P0(FF_sol, pr_sol, 0, FVV_sol)
Ptrial0_1 = P0(FF_sol, pr_sol, 1, FVV_sol)
mech_form0 = fn.inner(Ptrial_0, fn.grad(u_))*dx(0)  
mech_form1 = fn.inner(Ptrial_1, fn.grad(u_))*dx(1) 
form += mech_form0 + mech_form1  # Sum residual
fint.append(mech_form0 + mech_form1)

if magnetoactive:       # magnetic energy including H macro and H due to B macro
    Btrial_0=B(HH_sol,0,FF_sol)
    Btrial_1=B(HH_sol,1,FF_sol)
    mag_form0 = -fn.dot(Btrial_0,fn.grad(v_)) * dx(0)
    mag_form1 = -fn.dot(Btrial_1,fn.grad(v_)) * dx(1)
    form += mag_form0 + mag_form1
    fint.append(mag_form0 + mag_form1)

if incompressibility:   # incompressibility condition
    incompr_form0 = - pr_ * (fn.det(FF_sol) - 1) * dx(0)
    incompr_form1 = - pr_ * (fn.det(FF_sol) - 1) * dx(1)
    form += incompr_form0 + incompr_form1
    fint.append(incompr_form0 + incompr_form1)

if viscosity:       # viscous residual
    Rvisco0 = - 1e-3 * fn.dot(FV_sol - FV_sol_t - compare(FF_sol,0, FVV_sol, DT_macro), FV_) * dx(0)
    Rvisco1 = - 1e-3 * fn.dot(FV_sol - FV_sol_t - compare(FF_sol, 1, FVV_sol, DT_macro), FV_) * dx(1)
    form += Rvisco0 + Rvisco1
    fint.append(Rvisco0 + Rvisco1)
    
# imposition of macroscopic P for 2 subdomains
if np.sum(control_mech_macro)>0:
    therm_form0 = - fn.inner(lambda__tensor, (Ptrial0_0 - P_macro)) * dx(0)
    therm_form1 = - fn.inner(lambda__tensor, (Ptrial0_1 - P_macro)) * dx(1)
    form += therm_form0 + therm_form1
    fint.append(therm_form0 + therm_form1)
    
# imposition of macroscopic B for 2 subdomains
if magnetoactive and np.sum(control_mag_macro)>0:
    therm_formv0 = - fn.inner(lambdav__vector, (Btrial_0 - B_macro2)) * dx(0)
    therm_formv1 = - fn.inner(lambdav__vector, (Btrial_1 - B_macro2)) * dx(1)
    form += therm_formv0 + therm_formv1
    fint.append(therm_formv0 + therm_formv1)

# definiton of the global Jacobian
Jac = fn.derivative(form, W_sol, dW)

# ------------------------------------------------------------------------
#               DIRICHLET BOUNDARY CONDITIONS
# Define a function to identify the origin point (0, 0, 0) on the boundary
# Dirichlet Boundary Conditions:
# - Fix displacement (and magnetic potential, if applicable) at the origin (0, 0, 0) to prevent rigid body motion.
# - Applied pointwise to ensure stability in periodic domains.
# ------------------------------------------------------------------------


def bnd_func(x, on_boundary):
    """Check if a point is near the origin (0, 0, 0)."""
    return fn.near(x[0], 0) and fn.near(x[1], 0) and fn.near(x[2], 0)
bcss = [] 
# Dirichlet BC for displacement: Fix u = (0, 0, 0) at the origin
displacement_value = fn.Constant((0.0 + fn.DOLFIN_EPS, 0.0 + fn.DOLFIN_EPS, 0.0 + fn.DOLFIN_EPS))
bcsu = fn.DirichletBC(
    V.sub(0),            # Subspace for displacement (vector field)
    displacement_value,  # Value to enforce (near-zero vector)
    bnd_func,            # Function identifying the origin
    method="pointwise"   # Apply at specific points
)
bcss.append(bcsu)        # Add to boundary condition list

# Dirichlet BC for magnetic potential (if magnetoactive is True)
if magnetoactive:
    potential_value = fn.Constant(0.0 + fn.DOLFIN_EPS)
    bcsv = fn.DirichletBC(
        V.sub(1),        # Subspace for magnetic potential (scalar field)
        potential_value, # Value to enforce (near-zero scalar)
        bnd_func,        # Function identifying the origin
        method="pointwise"  # Apply at specific points
    )
    bcss.append(bcsv)    # Add to boundary condition list

for bc in bcss:
    print(f"BC DOFs: {bc.get_boundary_values().keys()}")

# ------------------------------------------------------------------------
#                   INITIALIZE VARIABLES AND POSTPROCESSING
# Initialize Variables and Postprocessing:
# - Arrays to store time and macroscopic quantities (deformation gradient F, stress P).
# - Convergence tracking variables to monitor residuals and adjust time steps.
# ------------------------------------------------------------------------
ttt = np.array([0])
H33 = np.array([0]) # Magnetic field Z direction is imposed therefore wants to save it
B33 = np.array([0]) # Magnetic induction Z direction is imposed therefore wants to save it
F11 = np.array([0]); F12=np.array([0]); F13=np.array([0]) 
F21 = np.array([0]); F22=np.array([0]); F23=np.array([0]) 
F31 = np.array([0]); F32=np.array([0]); F33=np.array([0]) 
P11 = np.array([0]); P12=np.array([0]); P13=np.array([0]) 
P21 = np.array([0]); P22=np.array([0]); P23=np.array([0]) 
P31 = np.array([0]); P32=np.array([0]); P33=np.array([0]) 

t=0
dtflag = False
diffW.vector()[:] = 0

qq_t = np.zeros(iiii+1); qq_t[:] = 1e-2
qqave = np.zeros(iiii+1); qqave[:] = 1e-2
qm_t = np.zeros(iiii+1); qm_t[:] = 1e-2
qmave = np.zeros(iiii+1); qmave[:] = 1e-2
allinactive = np.zeros(iiii+1, dtype  ='bool')
lowflux = np.zeros(iiii+1, dtype='bool')

# ------------------------------------------------------------------------
#                       Solution Loop
# - Iterate over time steps defined in `timer`.
# - Solve the nonlinear problem using Newton-Raphson with adaptive time stepping.
# - Save results to XDMF and text files, and plot stress vs. strain.
# ------------------------------------------------------------------------
file_results = XDMFFile(foldername+'/'+output_filename+'.xdmf') ##### ADDED JCL
for kstep in range(len(timer)):
    #initialize residuals
    qmcounter=np.zeros(iiii+1)
    qqcounter=np.zeros(iiii+1)
    qq_t[np.logical_or(qq_t < 1e-2, qm_t < 1e-2)] = 1e-2
    qm_t[np.logical_or(qq_t < 1e-2, qm_t < 1e-2)] = 1e-2

    # set timer
    dtime = timer[kstep, 0]
    t_end = timer[kstep, 1]
    dtmin = timer[kstep, 2]
    dtmax = timer[kstep, 3]
    dtold = 1 * dtime
    diffW2.vector()[:] = 0
  
    #time loop
    inc = 0
    while t < t_end-1e-8:
        inc += 1
        #update time and bcs
        if t + dtime > t_end : 
            dtime = t_end + fn.DOLFIN_EPS-t
        
        t += dtime
        F_macro.t = t
        P_macro.t = t
        B_macro.t = t
        H_macro.t = t
        DT_macro.dt = dtime

        #coupled solver
        iter_nw = 0
        deltaW_sol.vector()[:] = 1
        diffW.vector()[:] = 1
        error_nwall1 = np.ones(iiii+1)
        error_nwall2 = np.ones(iiii+1) #init errors

        while iter_nw < maxiter_nw:
            iter_nw += 1 
            print('Line 458. iter_nw :', iter_nw)       
            A, b = fn.assemble_system(Jac, -form, bcss) # assemble linearized system
            rres = b.get_local()                        # get residual
            #print('line 466. residual rres : ', [i for i in rres]) 
            error_nw_old1 = 1 * error_nwall1            #check tolerances
            error_nw_old2 = 1 * error_nwall2            #check tolerances  
            error_nwall1 = np.ones(iiii+1)              #init errors
            error_nwall2 = np.ones(iiii+1)              #init errors
            allinactive[:] = False
            lowflux[:] = False
            converged1 = np.zeros(iiii+1, dtype='bool')  #init convergence
            converged2 = np.zeros(iiii+1, dtype='bool')  #init convergence

            # for each field
            for ai in range(iiii+1):
                # get internal forces and compute relative residuals similar to abaqus documentation for nonlinear problems   
                q0 = np.abs(fn.assemble(fint[ai]).get_local()[dofindex[ai]])
                qm = np.max(q0)
                if qm > tol_nw_forces * qm_t[ai]:
                    qmave[ai] = qm_t[ai] * (qmcounter[ai]) / (qmcounter[ai] + 1) + qm * 1 / (qmcounter[ai] + 1)
                else: 
                    allinactive[ai] = True
                    qmave[ai] = qm_t[ai]

                if qm > 0.1 * qm_t[ai]:
                        qq = np.mean(q0[q0 >= tol_nw_forces * qm_t[ai]])
                else:
                        qq = np.mean(q0)
                Dumax = np.linalg.norm(diffW.vector()[dofindex[ai]][:], ord=np.inf)
                rmax = np.linalg.norm(rres[dofindex[ai]], ord=np.inf)
                dumax = np.linalg.norm(deltaW_sol.vector()[dofindex[ai]][:], ord=np.inf)
                if qq > tol_nw_forces * qq_t[ai]:
                    qqave[ai] = qq_t[ai] * (qqcounter[ai]) / (qqcounter[ai] + 1) + qq * 1 /(qqcounter[ai] + 1)
                    if rmax < 1e-8 * qqave[ai]:  
                        converged1[ai] = True
                        converged2[ai] = True
                        continue
                    error_nwall1[ai] = rmax / qqave[ai]
                    error_nwall2[ai] = dumax / Dumax

                    if np.abs(error_nwall1[ai]) < tol_nw: 
                        converged1[ai] = True

                    if np.abs(error_nwall2[ai]) < tol_nw or Dumax < 1e-8*0.2: 
                        converged2[ai] = True
                    continue

                else: 
                    lowflux[ai] = True
                    qqave[ai] = qq_t[ai]
                    if rmax < 1e-8 * qqave[ai]:  
                        converged1[ai] = True
                        converged2[ai] = True
                        continue
                    error_nwall1[ai] = tol_nw / tol_nw_forces * rmax / qm_t[ai]
                    error_nwall2[ai] = tol_nw / 1e-3 * dumax / Dumax
                    if np.abs(error_nwall1[ai]) < tol_nw: 
                        converged1[ai] = True
                        converged2[ai] = True
                        continue
                    else:
                        if np.abs(error_nwall2[ai]) < tol_nw or Dumax < 1e-8*0.2: 
                            converged1[ai] = True
                            converged2[ai] = True
                            continue
                        else:
                            continue

            # check if converged and stop iterating if so
            print('Iter: ', iter_nw-1, 'Coupled. Time :', t,'dtime: ', dtime)
            print('error_nwall1: ', error_nwall1, 'error_nwall2: ', error_nwall2)
            print('converged1:', converged1, 'converged2:', converged2) 
            print('np.logical_or', np.logical_or(converged1, converged2))

            if np.any(np.isnan(error_nwall1)) or np.any(np.isnan(error_nwall2)): 
                converged1[:] = False
                converged2[:] = False
                iter_nw = maxiter_nw
                print('line 533 - NaN in error_nwall')
                break

            if np.all(np.logical_and(converged1, converged2)) and iter_nw > 2:
                print('line 537') 
                break
            
            # if maximum iterations reached stop
            if iter_nw == maxiter_nw:
                print('line 572 - Max iter reached') 
                break

            # Solve the linear system and apply corrections for the next iteration
            try: 
                #fn.solve(A, deltaW_sol.vector(), b, "numps")
                solver.solve(A, deltaW_sol.vector(), b)
                residual = b - A * deltaW_sol.vector()
                print("Linear solver residual norm:", residual.norm('linf'))
            except Exception as e:
                print(f"Solver failed with error: {e}") 
                converged1[:] = False
                converged2[:] = False
                iter_nw = maxiter_nw
                break

            W_sol.vector()[:] = W_sol.vector() + deltaW_sol.vector()
            diffW.vector()[:] = W_sol.vector() - W_sol_t.vector()

        # Restore to time t if not converged 
        if iter_nw == maxiter_nw and not np.all(np.logical_and(converged1,converged2)):
            t -= dtime
            inc -= 1
            dtime = dtime/2.
            qqave[:] = qq_t[:]
            qmave[:] = qm_t[:]
            print('No convergency, decrease t :', t,'delta time:', dtime)
            W_sol.vector()[:] = W_sol_t.vector()
            if dtime < dtmin: 
                dtflag = True; print('dtime < dtmin (604)')
                break
            continue

        # postprocessing if converged
        print('-----------------------------------------------------------------')
        print('Converged with iters=',iter_nw, 'time: ', t,'delta time : ', dtime)
        
        #get some macroscopic quantities
        F11_t = fn.assemble((FF_sol - fn.Identity(3))[0,0] * dx)
        F12_t = fn.assemble((FF_sol - fn.Identity(3))[0,1] * dx)
        F13_t = fn.assemble((FF_sol - fn.Identity(3))[0,2] * dx)

        F21_t = fn.assemble((FF_sol - fn.Identity(3))[1,0] * dx)
        F22_t = fn.assemble((FF_sol - fn.Identity(3))[1,1] * dx)
        F23_t = fn.assemble((FF_sol - fn.Identity(3))[1,2] * dx)

        F31_t = fn.assemble((FF_sol - fn.Identity(3))[2,0] * dx)
        F32_t = fn.assemble((FF_sol - fn.Identity(3))[2,1] * dx)
        F33_t = fn.assemble((FF_sol - fn.Identity(3))[2,2] * dx)
        '''
        H11_t = fn.assemble(fn.dot(HH_sol, fn.inv(FF_sol))[0] * dx)
        H22_t = fn.assemble(fn.dot(HH_sol, fn.inv(FF_sol))[1] * dx)
        H33_t = fn.assemble(fn.dot(HH_sol, fn.inv(FF_sol))[2] * dx)

        B11_t=fn.assemble((fn.dot(FF_sol,B(HH_sol,0,FF_sol))/fn.det(FF_sol))[0]*dx(0)+\
                       (fn.dot(FF_sol,B(HH_sol,1,FF_sol))/fn.det(FF_sol))[0]*dx(1))
        B22_t=fn.assemble((fn.dot(FF_sol,B(HH_sol,0,FF_sol))/fn.det(FF_sol))[1]*dx(0)+\
                       (fn.dot(FF_sol,B(HH_sol,1,FF_sol))/fn.det(FF_sol))[1]*dx(1))
        B33_t=fn.assemble((fn.dot(FF_sol,B(HH_sol,0,FF_sol))/fn.det(FF_sol))[2]*dx(0)+\
                       (fn.dot(FF_sol,B(HH_sol,1,FF_sol))/fn.det(FF_sol))[2]*dx(1))
        '''
        P11_t = fn.assemble(P0(FF_sol, pr_sol, 0, FVV_sol)[0,0] * dx(0) +\
                P0(FF_sol, pr_sol, 1, FVV_sol)[0,0] * dx(1))
        P12_t = fn.assemble(P0(FF_sol, pr_sol, 0, FVV_sol)[0,1] * dx(0) +\
                        P0(FF_sol, pr_sol, 1, FVV_sol)[0,1] * dx(1))
        P13_t = fn.assemble(P0(FF_sol, pr_sol, 0, FVV_sol)[0,2] * dx(0) +\
                        P0(FF_sol, pr_sol, 1, FVV_sol)[0,2] * dx(1))
        
        P21_t = fn.assemble(P0(FF_sol, pr_sol, 0, FVV_sol)[1,0] * dx(0) +\
                P0(FF_sol, pr_sol, 1, FVV_sol)[1,0] * dx(1))
        P22_t = fn.assemble(P0(FF_sol, pr_sol, 0, FVV_sol)[1,1] * dx(0) +\
                P0(FF_sol, pr_sol, 1, FVV_sol)[1,2] * dx(1))
        P23_t = fn.assemble(P0(FF_sol, pr_sol, 0, FVV_sol)[1,2] * dx(0) +\
                        P0(FF_sol, pr_sol, 1, FVV_sol)[1,2] * dx(1))
        
        P31_t = fn.assemble(P0(FF_sol, pr_sol, 0, FVV_sol)[2,0] * dx(0) +\
                P0(FF_sol, pr_sol, 1, FVV_sol)[2,0] * dx(1))
        P32_t = fn.assemble(P0(FF_sol, pr_sol, 0, FVV_sol)[2,1] * dx(0) +\
                        P0(FF_sol, pr_sol, 1, FVV_sol)[2,1] * dx(1))
        P33_t = fn.assemble(P0(FF_sol, pr_sol, 0, FVV_sol)[2,2] * dx(0) +\
                        P0(FF_sol, pr_sol, 1, FVV_sol)[2,2] * dx(1))
        
        F11=np.append(F11, F11_t); F12=np.append(F12, F12_t); F13=np.append(F13, F13_t)
        F21=np.append(F21, F21_t); F22=np.append(F22, F22_t); F23=np.append(F23, F23_t)
        F31=np.append(F31, F31_t); F32=np.append(F32, F32_t); F33=np.append(F33, F33_t)

        P11=np.append(P11, P11_t); P12=np.append(P12, P12_t); P13=np.append(P13, P13_t)
        P22=np.append(P22, P22_t); P23=np.append(P23, P23_t); P21=np.append(P21, P21_t)
        P31=np.append(P31, P31_t); P32=np.append(P32, P32_t); P33=np.append(P33, P33_t)
        
        ttt=np.append(ttt, t)

        # update output files
        post(file_results, t, mesh, domains, u_sol, lambda_tensor, F_macro, dx, P0, FF_sol, pr_sol, FVV_sol) #### ADDED JCL

        #update solution and review residuals
        W_sol_t.vector()[:] = W_sol.vector()
        diffW2.vector()[:] = diffW.vector()
        diffW.vector()[:] = 0

        qq_t[:] = qqave[:]
        qm_t[:] = qmave[:]
        qqcounter[lowflux == False] += 1
        qmcounter[allinactive == False] += 1

        # adjust dt if converged easily
        if iter_nw <= 5: # Converged easily
            dtime = dtime * 1.5
            print('Converged easily. Increasing dtime to :', dtime)
            if dtime > dtmax:  # Ensure dtime not over dtmax
                dtime = 1 * dtmax 
                print('Increase dt',dtime)
            
        # stop simulation if minimum time increment requires
        if dtflag: 
            file_results.close() 
            print('dtflag is true. breaking the simulation (line 670)')
            break

# Saving data 
#np.savetxt(foldername+'/'+output_filename+'.txt', np.array([ttt, H33, F22, P33, B33, F12, P12]).transpose())
np.savetxt(foldername+'/'+output_filename+'.txt', 
    np.array([ttt, F33, P33, F11, F12, F13, F21, F22, F23, F31, F32, F33, P11, P12, P13, P21, P22, P23, P31, P32, P33]).transpose(),
    header='Time F33 P33 F11 F12 F13 F21 F22 F23 F31 F32 F33 P11 P12 P13 P21 P22 P23 P31 P32 P33',  # Column labels
    fmt='%.6e'  # Optional: format numbers in scientific notation with 6 decimal places
)


import matplotlib.pyplot as plt
# Plot P33 vs F33 directly from arrays
plt.figure(figsize=(8, 6))
plt.plot(P33, F33, 'b-', label='Stress vs Strain')
#plt.plot(ttt[::-1]/t_ramp, F33[::-1], 'b-', label='Strain vs. normalized time')
plt.xlabel('Strain')
plt.ylabel('Stress (kPa)')
plt.title('STRESS vs. STRAIN')
plt.grid(True)
plt.legend()

# Invert both axes
plt.gca().invert_xaxis()  # Invert x-axis (F33)
plt.gca().invert_yaxis()  # Invert y-axis (P33)

plt.show()