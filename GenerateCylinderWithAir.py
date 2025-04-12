import gmsh
import numpy as np

# Define parameters
cylinder_radius = 7.5           # Cylinder radius
cylinder_height = 6.25           # Cylinder height
# Air domain dimensions (box that contains the cylinder)
air_x = 30.0                    # Length in x-direction
air_y = 30.0                    # Length in y-direction
air_z = 30.0                     # Height in z-direction

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("cylinder_with_air")

# Create the cylinder (TAG = 1)
cylinder = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, cylinder_height, cylinder_radius, tag=1)

# Create the air domain (box) centered around the cylinder
air_box = gmsh.model.occ.addBox(-air_x/2, -air_y/2, -air_z/2, air_x, air_y, air_z, tag=2)

ov, ovv = gmsh.model.occ.fragment([(3, air_box)], [(3, cylinder)], removeObject=True, removeTool=True)
gmsh.model.occ.synchronize()


# Debug: Print fragmented volumes
print("Fragment produced volumes:")
for e in ov:
    print(e)

# Synchronize the geometry
gmsh.model.occ.synchronize()

# Define physical groups for volumes
air_tag = ov[-1][1]  # Last entity is the one being fragmented. Air box in this case. 
cylinder_tag = ov[0][1]  # All but the last are inclusions

# ---------------------------------------
# Assigning TAG to physical groups
# ---------------------------------------
gmsh.model.addPhysicalGroup(3, [cylinder_tag], tag=1, name="cylinder")    # Cylinder volume
gmsh.model.addPhysicalGroup(3, [air_tag], tag=2, name="air_domain")  # Air volume

# Define the boundary (2D surfaces of the air domain)
boundary_surfaces = gmsh.model.getBoundary([(3, air_tag)], oriented=False)
boundary_tags = [surf[1] for surf in boundary_surfaces]  # Extract surface tags
#gmsh.model.addPhysicalGroup(2, boundary_tags, tag=3, name="boundary")  # Cylinder surfaces as "boundary"


# Mesh sizing
lcar1 = 3  # General mesh size
lcar2 = 1  # Corner point size (optional)
lcar3 = 1  # Inclusion boundary size

# Assign mesh size to all points
gmsh.model.mesh.setSize(gmsh.model.getEntities(dim=0), size=lcar1)

# Refine near the cylinder. Cylinder tag=1
gmsh.model.mesh.setSize(
    gmsh.model.getBoundary(dimTags=[(3, 1)], combined=False, oriented=False, recursive=True),
    lcar3
)

# Generate the mesh
gmsh.model.mesh.generate(3)

# Save mesh
output_filename = "cylinder_with_air_domain_bulky"
gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write(f"{output_filename}.msh")

# Optional: Convert to XDMF/H5 (uncomment if you have the conversion function)
from GenerateXDMFH5 import generateXDMFH5
generateXDMFH5(msh_filename=f"{output_filename}.msh")

# Visualize (optional)
gmsh.fltk.run()

#-------------------#-------------------#-------------------#-------------------
