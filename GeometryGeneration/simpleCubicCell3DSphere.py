################################################################################
#            EXAMPLE FOR 3D SIMPLE CUBIC UNIT CELL  WITH SPHERES               #
################################################################################
# This example shows the generation of a unit cell with a simple cubic distribution
# of spherical inclusions. The basic procedures of the model an mesh generation
# are pointed out and the resulting mesh is visualized. For the example, only the
# standard configuration is used. However, in order to show all available options -
# user configurations are passed as dictionaries to the individual classes and
# methods - the dictionaries containing the default values are passed. This
# means that, if they were not passed, the resulting mesh would be the same.

# Loading of the SimpleCubicUnitCell class
# Before the model and mesh generation can start, the required class has to be
# loaded. In this case it is the class SimpleCubicCell
from gmshModel.Model import SimpleCubicCell
import numpy as np
'''
Initialization of the unit cell
In order to generate a mesh for unit cells with a simple cubic distribution of
spherical inclusions, relevant data have to be passed for the initialization of
a new object instance. For unit cells of the type under consideration, the
following parameters are possible:

radius: float (mandatory)
   radius of the inclusions within the unit cell

distance: float (defining either distance or size is mandatory)
   distance of the inclusions within the unit cell
   -> if the distance is given, the cells size is calculated automatically

size: list/array (defining either distance or size is mandatory)
   array defining the size of the RVE in the individual directions
   -> size=[L_x, L_y, (L_z)]
   -> if the size is given, the inclusion distances are calculated automatically
      (this allows more flexibility and unit cells with inclusion distributions
       that are similar to the physical unit cell under consideration)

numberCells: list/array (optional)
   array defining the number of cells in the 3 spatial axis directions
   -> numberCells=[n_x, n_y, n_z]

 inclusionType: string (mandatory)
   string defining the type of inclusions within the RVE

origin: list/array (optional)
   array defining the origin of the RVE
   -> origin=[O_x, O_y, (O_z)]

periodicityFlags: list/array (optional)
   array with flags (0/1) whether the current coordinate direction has to be
   treated as periodic
   periodicityFlags=[0/1, 0/1, 0/1]

domainGroup: string (optional)
   string defining which group the geometric objects defining the domain belong
   to (to reference this group within boolean operations)

inclusionGroup: string (optional)
   string defining which group the geometric objects defining the inclusions
   belong to (to reference this group within boolean operations)

 gmshConfigChanges: dict (optional)
   dictionary for user updates of the default Gmsh configuration
'''

msh_filename = "3spheres_cube_113_VF_01.msh"
cube_size = [1, 1, 1] # Cube dimension
volume_fraction = 0.1 # = Fraction of sphere_volume / cube_volume 

cube_volume = (cube_size[0]*cube_size[1]*cube_size[2])
sphere_volume = volume_fraction * cube_volume
radius = (3 * sphere_volume / (4 * np.pi)) ** (1/3) # Good only if periodicity is null 
print(f"Calculated sphere radius: {radius}")

initParameters={                                                                # save all possible parameters in one dict to facilitate the method call
    "numberCells": [cube_size[0], cube_size[1], cube_size[2]],                  # generate 1 unit cell in every spatial direction
    "radius": radius,                                                           # set the inclusion radius to 2
    "inclusionType": "Sphere",                                                  # define inclusionType as "Sphere"
    "distance": 1,                                                              # Must stay at 1
    "origin": [0, 0, 0],                                                        # set cell origin to [0,0,0]
    "periodicityFlags": [0,0,1],                                              # define all axis directions as periodic
    "domainGroup": "domain",                                                    # use "domain" as name for the domainGroup
    "inclusionGroup": "inclusions",                                             # use "inclusions" as name for the inclusionGroup
    "gmshConfigChanges": {"General.Terminal": 0,                                # deactivate console output by default (only activated for mesh generation)
                          "Mesh.CharacteristicLengthExtendFromBoundary": 0,     # do not calculate mesh sizes from the boundary by default (since mesh sizes are specified by fields)
                          "Mesh.ElementOrder": 1                                # Mesh order. Must be 1 with legacyFEniCS
    }
}
testCell=SimpleCubicCell(**initParameters)


# Gmsh model generation
testCell.createGmshModel()

# Gmsh mesh creation

meshingParameters={                                                             # save all possible parameters in one dict to facilitate the method call
    "threads": None,                                                            # do not activate parallel meshing by default
    "refinementOptions": {"maxMeshSize":  "auto",                                # automatically calculate maximum mesh size with built-in method
                          "inclusionRefinement": True,                          # flag to indicate active refinement of inclusions
                          "interInclusionRefinement": True,                     # flag to indicate active refinement of space between inclusions (inter-inclusion refinement)
                          "elementsPerCircumference": 10,                        # use 18 elements per inclusion circumference for inclusion refinement
                          "elementsBetweenInclusions": 5,                       # ensure 3 elements between close inclusions for inter-inclusion refinement
                          "inclusionRefinementWidth": 3,                        # use a relative (to inclusion radius) refinement width of 1 for inclusion refinement                  
                          "transitionElements": "auto",                         # automatically calculate number of transitioning elements (elements in which tanh function jumps from h_min to h_max) for inter-inclusion refinement
                          "algorithm3D": 7,                                      # algorithm for 3D meshing: 1 (Delaunay), 6 (frontal) or 7 (MMG3D)
                          "aspectRatio": 0.5,                                    # aspect ratio for inter-inclusion refinement: ratio of refinement in inclusion distance and perpendicular directions
                          "subdivisionAlgorithm": 3,                             # 1: Delaunay, 2: Frontal, 3: Gambit
                          "subdivisionSurfaceAlgorithm": 1,                      # 1: Loop, 2: Catmull-Clark
                          "subdivisionSurfaceOptimization": 1,                   # 1: Normal, 2: Quadric
                          "subdivisionVolumeAlgorithm": 1, # 1: Loop, 2: Catmull-Clark
                          "subdivisionVolumeOptimization": 1, # 1: Normal, 2: Quadric
    }
}
testCell.createMesh(**meshingParameters)


# Save resulting mesh to file

testCell.saveMesh(msh_filename)

# Show resulting mesh
testCell.visualizeMesh()

# Close Gmsh model
testCell.close()
from GenerateXDMFH5 import generateXDMFH5
generateXDMFH5(msh_filename=msh_filename) # will generate .xdmf and h5 files. 
