import meshio
import numpy as np
import os


def generateXDMFH5(msh_filename): 
# Input and output filenames
      # Your .msh file
    output_filename = os.path.splitext(msh_filename)[0]

    #Read the .msh file
    mesh = meshio.read(msh_filename)

    # Inspect the mesh structure
    print("Cell types:", [cb.type for cb in mesh.cells])
    print("mesh.cell_data keys:", list(mesh.cell_data.keys()))

    # Check for physical group data
    if "gmsh:physical" not in mesh.cell_data:
        raise KeyError("No 'gmsh:physical' data found in .msh file. Physical groups may not be defined.")

    # Extract physical data
    physical_data = mesh.cell_data["gmsh:physical"]
    print("Type of physical_data:", type(physical_data))
    print("Length of physical_data:", len(physical_data))

    # Collect all 'tetra' cells and their physical tags
    tetra_cells = []
    tetra_physical_tags = []

    for i, cell_block in enumerate(mesh.cells):
        if cell_block.type == "tetra":
            tetra_cells.append(cell_block.data)
            tetra_physical_tags.append(physical_data[i])

    if not tetra_cells:
        raise ValueError("No 'tetra' cells found in the .msh file.")

    # Combine all 'tetra' cells and tags
    tetra_cells = np.vstack(tetra_cells)
    tetra_physical_tags = np.hstack(tetra_physical_tags)

    print(f"Total number of 'tetra' elements: {len(tetra_cells)}")
    print(f"Length of physical tags for 'tetra': {len(tetra_physical_tags)}")

    # Create a new mesh with only 'tetra' elements
    tetra_mesh = meshio.Mesh(
        points=mesh.points,
        cells=[("tetra", tetra_cells)],
        cell_data={"phase": [tetra_physical_tags]}
    )

    # Write to XDMF with HDF5 storage
    xdmf_file = output_filename + ".xdmf"
    h5_file = output_filename + ".h5"
    meshio.write(xdmf_file, tetra_mesh, data_format="HDF")
    print(f"Generated {xdmf_file} and {h5_file}")