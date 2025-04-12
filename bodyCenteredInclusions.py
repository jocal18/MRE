#import numpy as np
from gmshModel.Model import BodyCenteredCubicCell  # Replace with actual module import
import numpy as np
import gmsh

msh_filename =  'BCC_VF03_v0.msh'
# Define input parameters
cube_size = [1, 1, 1]      # Total cube dimensions (e.g., 1x1x1 units)
volume_fraction = 0.3      # Target volume fraction of inclusions

# Calculate total cube volume
cube_volume = cube_size[0] * cube_size[1] * cube_size[2]

# BCC unit cell properties
# In BCC, there are 2 inclusions per unit cell (1 center + 8 * 1/8 corners)
# Unit cell size is "distance" (a cube of side length d), and volume is d^3
# Volume of one sphere = (4/3) * pi * r^3
# Volume fraction per unit cell = (2 * sphere_volume) / (distance^3)

# Choose a reasonable distance (e.g., based on cube size or arbitrary initial guess)
distance = 1/2 # Distance between inclusions (adjustable; here, half the cube size for simplicity)
unit_cell_volume = distance ** 3

# Calculate required radius for target volume fraction
# volume_fraction = (2 * (4/3) * pi * r^3) / distance^3
# Solving for r: r^3 = (volume_fraction * distance^3) / (2 * (4/3) * pi)
sphere_volume_fraction = volume_fraction / 2  # Two spheres per unit cell
inclusion_radius = ((3 * volume_fraction * unit_cell_volume) / (8 * np.pi)) ** (1/3)

# Number of unit cells to fill cube_size
number_cells = [int(cube_size[0] / distance), int(cube_size[1] / distance), int(cube_size[2] / distance)]

print(f"Unit cell distance: {distance}")
print(f"Inclusion radius: {inclusion_radius:.4f}")
print(f"Number of unit cells: {number_cells}")
print(f"Effective volume fraction: {(2 * (4/3) * np.pi * inclusion_radius**3) / unit_cell_volume:.4f}")

# Define parameters based on volumetric fraction
initParameters = {
    "numberCells": number_cells,  # Number of BCC unit cells in each direction
    "radius": inclusion_radius,   # Radius calculated from volume fraction
    "distance": distance,         # Distance between inclusions (unit cell size)
    "inclusionType": "Sphere",    # Inclusion type
    "origin": [0, 0, 0],          # Origin at (0, 0, 0)
    "periodicityFlags": [1, 1, 1], # Periodic in all directions
    "domainGroup": "domain",       # Name for the domain group
    "inclusionGroup": "inclusions",# Name for the inclusion group
    "gmshConfigChanges": {
        "General.Terminal": 0,                            # No console output
        "Mesh.CharacteristicLengthExtendFromBoundary": 0, # Fixed mesh size
    }
}

# Create the BCC RVE
testRVE = BodyCenteredCubicCell(**initParameters)

testRVE.createGmshModel()

#  Calculate volume_fraction after placement
#n_spheres_final = round(testRVE.placementInfo[0])    # Return final number of inclusions 
#volume_fraction_final = (n_spheres_final * sphere_volume) / cube_volume
#volume_fraction_final = round(volume_fraction_final, 2)
#decimal_part = int((volume_fraction_final % 1) * 100)  # Extract the first 2 decimal digits
#decimal_part2 = int((inclusion_radius % 1) * 100)  # Extract the first 2 decimal digits
# Create new filename with the 2 decimal digits
# f"model-RI-VF-0{decimal_part}-N_{n_spheres_final}-R_{decimal_part2}.msh"

# Gmsh mesh creation
meshingParameters={                                                             # save all possible parameters in one dict to facilitate the method call
    "threads": 1,                                                               # do not activate parallel meshing by default
    "refinementOptions": {
                          "maxMeshSize": 0.1,                                   # Maximum mesh size
                          "inclusionRefinement": True,                          # flag to indicate active refinement of inclusions
                          "interInclusionRefinement": False,                     # flag to indicate active refinement of space between inclusions (inter-inclusion refinement)
                          "elementsPerCircumference": 6, #18                       # use X elements per inclusion circumference for inclusion refinement
                          "transitionElements": "auto",                          # automatically calculate number of transitioning elements (elements in which tanh function jumps from h_min to h_max) for inter-inclusion refinement
                          "algorithm3D": 7,                                      # algorithm for 3D meshing: 1 (Delaunay), 6 (frontal) or 7 (MMG3D)
                          "aspectRatio": 1,                                      # aspect ratio for inter-inclusion refinement: ratio of refinement in inclusion distance and perpendicular directions
                          "subdivisionAlgorithm": 3,                             # 1: Delaunay, 2: Frontal, 3: MMG3D
                          "subdivisionSurfaceAlgorithm": 1,                      # 1: Loop, 2: Catmull-Clark
                          "subdivisionSurfaceOptimization": 1,                   # 1: Normal, 2: Quadric
                          "subdivisionVolumeAlgorithm": 1, # 1: Loop, 2: Catmull-Clark
                          "subdivisionVolumeOptimization": 1, # 1: Normal, 2: Quadric
    }
}
testRVE.createMesh(**meshingParameters)


# Save resulting mesh to file
testRVE.saveMesh(msh_filename)

# Visualize (optional)
gmsh.fltk.run()
testRVE.visualizeMesh()

# Close Gmsh model
testRVE.close()

from GenerateXDMFH5 import generateXDMFH5
generateXDMFH5(msh_filename=msh_filename) # will generate .xdmf and h5 files. 

