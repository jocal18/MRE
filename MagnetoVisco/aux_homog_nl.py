#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    FEniCS complete code for magneto-visco-mechanical homogenization: 3D full-field magneto-mechanical 
    homogenization that accounts for: finite strains, viscoelasticity, incompressibility and mixed 
    stress/magnetic flux control.

    Copyright (C) 2021: Sergio Lucarini, Miguel Angel Moreno-Mateos, Kostas Danas, Daniel Garcia-Gonzalez.

    If using this code for research or industrial purposes, please cite:
    S. Lucarini, M.A. Moreno-Mateos, K. Danas, D. Garcia-Gonzalez
    Insights into the viscohyperelastic response of soft magnetorheological
    elastomers: competition of macrostructural versus microstructural players
    International Journal of Solids and Structures, 2022.
"""

import numpy as np
import fenics as fn
from fenics import * 

# funtion for reading the xdmf with subdomains
def importMesh(file):
    ## Import mesh
    mesh = fn.Mesh()
    with fn.XDMFFile(file) as infile: 
        infile.read(mesh)
    ## Import material info (physical volume)
    mvc = fn.MeshValueCollection("size_t", mesh, 3)
    with fn.XDMFFile(file) as infile: 
        infile.read(mvc, "phase")
    materials = fn.MeshFunction("size_t",mesh, mvc)
    return mesh,materials,mvc

# class definition for periodic boundary conditionss from the vetices
class PeriodicBoundary(fn.SubDomain):
    def __init__(self, vertices):
        # vertices stores the coordinates of the 8 unit cell corners
        fn.SubDomain.__init__(self)
        self.vv = vertices
        self.a1 = self.vv[1,:]-self.vv[0,:] # first vector generating periodicity
        self.a2 = self.vv[3,:]-self.vv[0,:] # second vector generating periodicity
        self.a3 = self.vv[4,:]-self.vv[0,:] # third vector generating periodicity

    def inside(self, x, on_boundary):
        # return True if on left or bottom or sleft boundary AND NOT on one of the 
        # bottom-right-sleft or bottom-left-srigth or top-left-sleft  vertices
        return bool(on_boundary and \
         ((fn.near(x[0], self.vv[0,0]) and (not (fn.near(x[1], self.vv[3,1]) or fn.near(x[2], self.vv[4,2])))) or \
          (fn.near(x[1], self.vv[0,1]) and (not (fn.near(x[0], self.vv[1,0]) or fn.near(x[2], self.vv[4,2])))) or \
          (fn.near(x[2], self.vv[0,2]) and (not (fn.near(x[0], self.vv[1,0]) or fn.near(x[1], self.vv[3,1]))))))

    def map(self, x, y):
        if fn.near(x[0], self.vv[6,0]) and fn.near(x[1], self.vv[6,1]) and fn.near(x[2], self.vv[6,2]): # if on top-right-sright corner
            y[0] = x[0] - (self.a1[0]+self.a2[0]+self.a3[0])
            y[1] = x[1] - (self.a1[1]+self.a2[1]+self.a3[1])
            y[2] = x[2] - (self.a1[2]+self.a2[2]+self.a3[2])
        elif fn.near(x[0], self.vv[2,0]) and fn.near(x[1], self.vv[2,1]): # if on top-right corner
            y[0] = x[0] - (self.a1[0]+self.a2[0])
            y[1] = x[1] - (self.a1[1]+self.a2[1])
            y[2] = x[2] - (self.a1[2]+self.a2[2])
        elif fn.near(x[0], self.vv[5,0]) and fn.near(x[2], self.vv[5,2]): # if on top-sright corner
            y[0] = x[0] - (self.a1[0]+self.a3[0])
            y[1] = x[1] - (self.a1[1]+self.a3[1])
            y[2] = x[2] - (self.a1[2]+self.a3[2])
        elif fn.near(x[1], self.vv[7,1]) and fn.near(x[2], self.vv[7,2]): # if on right-sright corner
            y[0] = x[0] - (self.a3[0]+self.a2[0])
            y[1] = x[1] - (self.a3[1]+self.a2[1])
            y[2] = x[2] - (self.a3[2]+self.a2[2])
        elif fn.near(x[0], self.vv[1,0]): # if on right boundary
            y[0] = x[0] - self.a1[0]
            y[1] = x[1] - self.a1[1]
            y[2] = x[2] - self.a1[2]
        elif fn.near(x[1], self.vv[3,1]): # if on sright boundary
            y[0] = x[0] - self.a2[0]
            y[1] = x[1] - self.a2[1]
            y[2] = x[2] - self.a2[2]
        elif fn.near(x[2], self.vv[4,2]): # should be on top boundary
            y[0] = x[0] - self.a3[0]
            y[1] = x[1] - self.a3[1]
            y[2] = x[2] - self.a3[2]
        else:
            y[0] = x[0]
            y[1] = x[1]
            y[2] = x[2]


def post(file_results, t, mesh, domains, u_sol, lambda_tensor, F_macro, dx, P0, FF_sol, pr_sol, FVV_sol):
    ''''
    Postprocessing
    Required to save results as XDMF file to be analysed with Paraview'''
    # Configure XDMF file parameters
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    file_results.parameters["rewrite_function_mesh"] = False

    # Write mesh and domains only once at t=0
    if t == 0:
        file_results.write(mesh, t)
        file_results.write(domains, t)

    # Fluctuating displacement field
    V20 = fn.VectorFunctionSpace(mesh, "CG", 2)
    u_sol1 = fn.project(u_sol, V20, solver_type='cg', preconditioner_type='hypre_amg')
    u_sol1.rename("Displacement_fluc", "u_sol1")
    file_results.write(u_sol1, t)

    # Total displacement field
    macro_disp = fn.dot(lambda_tensor + F_macro, Expression(("x[0]", "x[1]", "x[2]"), degree=1))
    u_tot = fn.project(macro_disp + u_sol, V20, solver_type='cg', preconditioner_type='hypre_amg')
    u_tot.rename("Displacement", "u_tot")
    file_results.write(u_tot, t)

    # Strain field (F - I)
    V3 = fn.TensorFunctionSpace(mesh, "DG", 1)
    e_tot = fn.project(FF_sol - Identity(3), V3, solver_type='cg', preconditioner_type='hypre_amg')
    e_tot.rename("Strain", "e_tot")
    file_results.write(e_tot, t)

    # Stress field
    def split_project(P0, FF_sol, pr_sol, FVV_sol, V3):
        u = fn.TrialFunction(V3)
        v = fn.TestFunction(V3)
        s_tot = fn.Function(V3)
        a = fn.inner(u, v) * dx
        L = fn.inner(P0(FF_sol, pr_sol, 0, FVV_sol), v) * dx(0) + fn.inner(P0(FF_sol, pr_sol, 1, FVV_sol), v) * dx(1)
        solve(a == L, s_tot, solver_parameters={"linear_solver": "cg", "preconditioner": "hypre_amg"})
        return s_tot
    
    s_tot = split_project(P0, FF_sol, pr_sol, FVV_sol, V3)
    s_tot.rename("Stress", "s_tot")
    file_results.write(s_tot, t)

    return

def adjust_xdmf(filenamexdmf):
    with open(filenamexdmf, 'r') as f:
        lines = f.readlines()
    
    # Find an existing attribute (e.g., Displacement)
    attr_lines = []
    for i in range(len(lines)):
        if lines[i].startswith('      <Attribute Name="Displacement"'):
            attr_lines = ['  ' + lines[i], '  ' + lines[i+1], '  ' + lines[i+2]]
            break
    
    # Rewrite file, inserting attribute after each time step
    with open(filenamexdmf, 'w') as f2:
        for line in lines:
            f2.write(line)
            if line.startswith('        <Time Value="') and attr_lines:
                f2.writelines(attr_lines)
    
    return



# ---------------------------------------------------
#               GRAVEYARD 
# ---------------------------------------------------
'''
# posprocessing: xdmf file
#def post(file_results,t,mesh,domains,u_sol,lambdas_tensor,Fmacro,P,dx,FF_sol):
    # Comment JCL 
    # u_sol : solution deplacement
def post2(file_results,t,mesh,domains,u_sol,lambda_tensor,F_macro,dx, P0, FF_sol, pr_sol, FVV_sol):

    #Postprocessing
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    file_results.parameters["rewrite_function_mesh"] = False
    if t==0:    # Write mesh only once 
        file_results.write(mesh)
        file_results.write(domains)
    
    # fluctuating displacement field
    V20 = fn.VectorFunctionSpace(mesh, "CG", 2)
    u_sol1=fn.project(u_sol,V20, solver_type='cg', preconditioner_type='hypre_amg')
    u_sol1.rename("Displacement_fluc", "u_sol1")
    file_results.write(u_sol1,t)

    # total displacement field
    u_tot=fn.project(fn.project(fn.dot(lambda_tensor+F_macro,\
      fn.Expression(("x[0]","x[1]","x[2]"), degree=1)),\
      V20, solver_type='cg', preconditioner_type='hypre_amg')+u_sol,\
      V20, solver_type='cg', preconditioner_type='hypre_amg')
    u_tot.rename("Displacement", "u_tot")
    file_results.write(u_tot,t=t)

    # F field
    V3 = fn.TensorFunctionSpace(mesh, "DG", 1)
    e_tot=fn.project(FF_sol-fn.Identity(3),V3,\
      solver_type='cg', preconditioner_type='hypre_amg')
    e_tot.rename("Strain", "e_tot")
    file_results.write(e_tot,t=t)

    
    # P field
    def split_project(P0,u_sol,V3):
        u = fn.TrialFunction(V3)
        v = fn.TestFunction(V3)
        s_tot = fn.Function(V3)
        fn.solve(fn.inner(u,v)*dx == fn.inner(P0(FF_sol, pr_sol, 0, FVV_sol),v)*dx(0)\
          + fn.inner(P0(FF_sol, pr_sol, 1, FVV_sol),v)*dx(1),s_tot, solver_parameters={"linear_solver": "cg"})
        return s_tot
    s_tot=split_project(P0,u_sol,V3)
    s_tot.rename("Stress", "s_tot")
    file_results.write(s_tot,t=t)
    
    return 


def adjust_xdmf3(filenamexdmf):
    f = open(filenamexdmf, 'r').readlines()
    f2 = open(filenamexdmf, 'w')
    # Initialize with empty strings or placeholders
    line1 = line2 = line3 = ''
    for i in range(len(f)):
        if f[i][0:25] == '      <Attribute Name="f"':
            line1 = '  ' + f[i]
            line2 = '  ' + f[i+1]
            line3 = '  ' + f[i+2]
        #if i > 3 and i < 15: 
        #   continue
        f2.write(f[i])
        if f[i][0:21] == '        <Time Value="':
            f2.write(line1)
            f2.write(line2)
            f2.write(line3)
    f2.close()
    return 

# funtion to correct paraview visualization modifying the xdmf
def adjust_xdmf2(filenamexdmf):
    f=open(filenamexdmf,'r').readlines()
    f2=open(filenamexdmf,'w')
    line1 = line2 = line3 = ''
    for i in range(len(f)):
        if f[i][0:25]=='      <Attribute Name=\"f\"':
            line1='  '+f[i]
            line2='  '+f[i+1]
            line3='  '+f[i+2]
        if i>3 and i<15: 
            continue
        f2.write(f[i])
        if f[i][0:21]=='        <Time Value="':
            f2.write(line1)
            f2.write(line2)
            f2.write(line3)
        i+=1
    f2.close()
    return
'''