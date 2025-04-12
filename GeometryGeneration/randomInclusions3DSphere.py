################################################################################
#              EXAMPLE FOR 3D RANDOMINCLUSION RVE WITH SPHERES                 #
# https://github.com/NEFM-TUDresden/gmshModel/blob/master/examples/randomInclusions3DSphere.py 
################################################################################
# This example shows the generation of an RVE with randomly placed, spherical
# inclusions. The basic procedure of the model an mesh generation are pointed
# out and the resulting mesh is visualized. For the example, only the standard
# configuration is used. However, in order to show the available options - all
# user configurations are passed as dictionaries to the individual classes and
# methods - the dictionaries containing the default values are passed. This
# means that, if they were not passed, the resulting mesh would be the same.

# Loading of the RandomInclusionRVE class
# Before the model and mesh generation can start, the required class has to be
# loaded. In this case it is the class RandomInclusionRVE
from gmshModel.Model import RandomInclusionRVE
import numpy as np


'''
Initialization of the RVE
In order to generate a mesh for RVEs with randomly placed inclusions, relevant
data have to be passed for the initialization of a new object instance. For
RVEs of the type under consideration, the following parameters are possible:

size: list/array (mandatory)
   array defining the size of the RVE in the individual directions
   -> size=[L_x, L_y, L_z]

inclusionSets: list/array (mandatory)
   array defining the relevant information (radius and amount) for the individual
   groups of spherical inclusions to be placed
   -> inclusionSets=[[r_1, n_1] [r_2, n_2], ..., [r_n, n_n]]

inclusionType: string (mandatory)
   string defining the type of inclusions within the RVE

origin: list/array (optional)
   array defining the origin of the RVE
   -> origin=[O_x, O_y, O_z]

periodicityFlags: list/array (optional)
   array with flags (0/1) whether the current coordinate direction has to be
   treated as periodic
   -> periodicityFlags=[0/1, 0/1, 0/1] 

domainGroup: string (optional)
   string defining which group the geometric objects defining the domain belong
   to (to reference this group within boolean operations)

inclusionGroup: string (optional)
   string defining which group the geometric objects defining the inclusions
   belong to (to reference this group within boolean operations)

gmshConfigChanges: dict (optional)
   dictionary for user updates of the default Gmsh configuration
'''

cube_size = [1, 1, 1]      # Cube dimension
inclusion_radius = 0.05    # Define inclusions radius. 
volume_fraction = 0.3      # = Fraction of sphere_volume / cube_volume 

cube_volume = (cube_size[0]*cube_size[1]*cube_size[2])
sphere_volume = (4 / 3) * np.pi * inclusion_radius ** 3 
n_spheres = int((volume_fraction * cube_volume) / sphere_volume)  # Approx. number based on full spheres
print(f"Target number of spheres (approx.): {n_spheres}")

initParameters={                                                                # save all possible parameters in one dict to facilitate the method call
    "inclusionSets": [inclusion_radius, n_spheres],                                                # [radius, number of inclusion]
    "inclusionType": "Sphere",                                                  # define inclusionType as "Sphere"
    "size": [cube_size[0], cube_size[1], cube_size[2]],                         # set RVE size 
    "origin": [0, 0, 0],                                                        # set RVE origin to [0,0,0]
    "periodicityFlags": [1, 1, 1],                                              # define all axis directions as periodic
    "domainGroup": "domain",                                                    # use "domain" as name for the domainGroup
    "inclusionGroup": "inclusions",                                             # use "inclusions" as name for the inclusionGroup
    "gmshConfigChanges": {"General.Terminal": 0,                                # deactivate console output by default (only activated for mesh generation)
                          "Mesh.ElementOrder": 1,                                # set element order to 1:linear 
                          "Mesh.CharacteristicLengthExtendFromBoundary": 0     # do not calculate mesh sizes from the boundary by default (since mesh sizes are specified by fields)
                          }
}
testRVE=RandomInclusionRVE(**initParameters)


# Gmsh model generation
modelingParameters={                                                            # save all possible parameters in one dict to facilitate the method call
    "placementOptions": {"maxAttempts": 10000,                                  # maximum number of attempts to place one inclusion
                         "minRelDistBnd": 1,                                  # minimum relative (to inclusion radius) distance to the domain boundaries
                         "minRelDistInc": 0.25,                                  # minimum relative (to inclusion radius) distance to other inclusions}
    }
}
testRVE.createGmshModel(**modelingParameters)

#  Calculate volume_fraction after placement
n_spheres_final = round(testRVE.placementInfo[0])    # Return final number of inclusions 
volume_fraction_final = (n_spheres_final * sphere_volume) / cube_volume
volume_fraction_final = round(volume_fraction_final, 2)
decimal_part = int((volume_fraction_final % 1) * 100)  # Extract the first 2 decimal digits
decimal_part2 = int((inclusion_radius % 1) * 100)  # Extract the first 2 decimal digits
# Create new filename with the 2 decimal digits
msh_filename = f"model-RI-VF-0{decimal_part}-N_{n_spheres_final}-R_{decimal_part2}.msh"

# Gmsh mesh creation
meshingParameters={                                                             # save all possible parameters in one dict to facilitate the method call
    "threads": 1,                                                               # do not activate parallel meshing by default
    "refinementOptions": {"maxMeshSize": "auto",                                # automatically calculate maximum mesh size with built-in method
                          "inclusionRefinement": True,                          # flag to indicate active refinement of inclusions
                          "interInclusionRefinement": True,                     # flag to indicate active refinement of space between inclusions (inter-inclusion refinement)
                          "elementsPerCircumference": 6, #18                       # use X elements per inclusion circumference for inclusion refinement
                          "elementsBetweenInclusions": 3, #5,                       # ensure 3 elements between close inclusions for inter-inclusion refinement
                          "inclusionRefinementWidth": 1,#1,                        # use a relative (to inclusion radius) refinement width of 1 for inclusion refinement
                          "transitionElements": "auto",                          # automatically calculate number of transitioning elements (elements in which tanh function jumps from h_min to h_max) for inter-inclusion refinement
                          "algorithm3D": 7,                                      # algorithm for 3D meshing: 1 (Delaunay), 6 (frontal) or 7 (MMG3D)
                          "aspectRatio": 2,                                      # aspect ratio for inter-inclusion refinement: ratio of refinement in inclusion distance and perpendicular directions
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

# Show resulting mesh
testRVE.visualizeMesh()

# Close Gmsh model
testRVE.close()

from GenerateXDMFH5 import generateXDMFH5
generateXDMFH5(msh_filename=msh_filename) # will generate .xdmf and h5 files. 

