import gmsh
import numpy as np
import meshio

# Define parameters
cylinder_radius = 1.0           # Cylinder radius
cylinder_height = 3.0           # Cylinder height
volume_fraction = 0.2           # Volume fraction of inclusions (20%)
inclusion_radius = 0.2        # Radius of each spherical inclusion

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("cylinder_with_inclusions")

# Create the cylinder (TAG = 1)
cylinder = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, cylinder_height, cylinder_radius, tag=1)

# Calculate the number of inclusions
V_cylinder = np.pi * cylinder_radius**2 * cylinder_height
V_inclusions = volume_fraction * V_cylinder
V_one = (4/3) * np.pi * inclusion_radius**3
n_inclusions = int(V_inclusions / V_one)

print(f"Number of inclusions to generate: {n_inclusions}")
k = 1.2               # Minimal relative distance factor
max_attempts = 1000   # To prevent infinite loops
inclusion_list = []
existing_positions = []  # Track positions of placed inclusions
min_distance = 2 * k * inclusion_radius  # Minimal distance between sphere centers

# Function to check if a new position is far enough from existing ones
def is_valid_position(new_pos, existing_positions, min_distance):
    for pos in existing_positions:
        distance = np.linalg.norm(new_pos - pos)
        if distance < min_distance:
            return False
    return True

# Generate inclusions
for t in range(1, n_inclusions + 1):
    attempts = 0
    while attempts < max_attempts:
        x = np.random.uniform(-1, 1)
        y = np.random.uniform(-1, 1)
        z = np.random.uniform(inclusion_radius, cylinder_height - inclusion_radius)
        if x**2 + y**2 <= (cylinder_radius - inclusion_radius)**2:
            new_pos = np.array([x, y, z])
            if is_valid_position(new_pos, existing_positions, min_distance):
                gmsh.model.occ.addSphere(x, y, z, inclusion_radius, t + 1)
                inclusion_list.append((3, t + 1))
                existing_positions.append(new_pos)
                break
        attempts += 1
    if attempts == max_attempts:
        print(f"Warning: Could not place inclusion {t} after {max_attempts} attempts.")

# Fragment the cylinder with inclusions
ov, ovv = gmsh.model.occ.fragment([(3, cylinder)], inclusion_list, removeObject=True, removeTool=True)
gmsh.model.occ.synchronize()

# Debug: Print fragmented volumes
print("Fragment produced volumes:")
for e in ov:
    print(e)

# Identify the new cylinder (matrix) and inclusions
new_cylinder_tag = ov[-1][1]  # Last entity is typically the matrix
inclusion_tags = [vol[1] for vol in ov[:-1]]  # All but the last are inclusions

# Define physical groups for volumes
gmsh.model.addPhysicalGroup(3, [new_cylinder_tag], tag=1, name="domain")  # Matrix as "domain"
gmsh.model.addPhysicalGroup(3, inclusion_tags, tag=2, name="inclusions")  # Inclusions as "inclusions"

# Define the boundary (2D surfaces of the cylinder)
boundary_surfaces = gmsh.model.getBoundary([(3, new_cylinder_tag)], oriented=False)
#boundary_surfaces = gmsh.model.getBoundary(inclusion_tags, oriented=False)
boundary_tags = [surf[1] for surf in boundary_surfaces]  # Extract surface tags
gmsh.model.addPhysicalGroup(2, boundary_tags, tag=3, name="boundary")  # Cylinder surfaces as "boundary"

# Mesh sizing
lcar1 = 0.25  # General mesh size
lcar2 = 0.25  # Corner point size (optional)
lcar3 = 0.15  # Inclusion boundary size

# Assign mesh size to all points
gmsh.model.mesh.setSize(gmsh.model.getEntities(dim=0), size=lcar1)

# Refine near inclusions
gmsh.model.mesh.setSize(
    gmsh.model.getBoundary(dimTags=[(3, i) for i in inclusion_tags], combined=False, oriented=False, recursive=True),
    lcar3
)

# Optional: Refine a point near the center (e.g., for testing)
eps = 1e-3
ov = gmsh.model.getEntitiesInBoundingBox(0.5 - eps, 0.5 - eps, 0.5 - eps, 0.5 + eps, 0.5 + eps, 0.5 + eps, 0)
gmsh.model.mesh.setSize(ov, lcar2)

# Generate the mesh
gmsh.model.mesh.generate(3)

# Save mesh
output_filename = "cylinder_with_inclusions"
gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write(f"{output_filename}.msh")

# Optional: Convert to XDMF/H5 (uncomment your conversion function if needed)
from GenerateXDMFH5 import generateXDMFH5
generateXDMFH5(msh_filename=f"{output_filename}.msh")

# Visualize (optional)
gmsh.fltk.run()

# Finalize Gmsh
gmsh.finalize()