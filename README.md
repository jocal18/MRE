# Magneto-Mechanical Simulations of Magnetorheological Elastomers (MREs)

This repository contains scripts, geometries, and simulation results associated with computational studies on magnetorheological elastomers (MREs). Due to file size constraints, some results and model geometries are hosted externally.

## 📁 Repository Structure

- `/model`: Contains the geometric models used in the simulations.
- `/results`: Contains representative simulation outputs (.xdmf format).
- `*.py`: Python scripts used for geometry generation, simulation, and post-processing.

---

## 🛠 Tools and Installation

### Prerequisites

This repository assumes familiarity with:

- Python programming
- Virtual environments via [Anaconda](https://www.anaconda.com/)
- Working with VS Code or equivalent IDEs (optional but recommended)

### Recommended Environment Setup

A step-by-step installation guide for **FEniCSx** is available at:  
👉 [https://me.jhu.edu/nguyenlab/doku.php?id=fenicsx](https://me.jhu.edu/nguyenlab/doku.php?id=fenicsx)  
This guide is particularly helpful for configuring FEniCSx within an Anaconda virtual environment.

> Note: **FEniCS** is the legacy version, while **FEniCSx** is the modern, actively maintained framework. The legacy version remains available and may be suitable for reproducing older research work.

---

## 📚 Useful Resources

- **Gmsh**: Open-source 3D finite element mesh generator with a built-in CAD engine and post-processor  
  - Install via pip:  
    ```bash
    pip install gmsh
    ```
  - ⚠️ The `gmsh` package from `conda-forge` may not include the Python API. Prefer `pip` installation, or alternatively:  
    ```bash
    conda install -c conda-forge gmsh python-gmsh
    ```
  - [Gmsh Installation Notes (StackOverflow)](https://stackoverflow.com/questions/70947216/python-cant-use-gmsh-after-installing-via-pip-and-conda)

- **GmshModel**: Python interface for full Gmsh functionality  
- **FEniCS/FEniCSx**: Finite Element Method framework for solving PDEs  
  - GitHub: [https://github.com/FEniCS](https://github.com/FEniCS)
- **ParaView**: Open-source visualization tool for `.xdmf` and other mesh formats

---

## 📄 Research Paper 1

**Title**: *A Computational Framework for Magnetically Hard and Soft Viscoelastic Magnetorheological Elastomers*  
**DOI**: [10.1016/j.ijsolstr.2013.08.024](https://doi.org/10.1016/j.ijsolstr.2013.08.024)  
**Code Repository**: [Zenodo Record](https://zenodo.org/records/5543516)

### Associated Files

**Main Simulation Directory**: `/MagnetoVisco/`

- `fenics_magneto_visco_final.py`: 3D full-field magneto-mechanical homogenization under finite strains, incompressibility, viscoelasticity, and mixed control conditions.
- `aux_homog_nl.py`: Post-processing and `.xdmf` file generation for ParaView visualization.

**Geometry Generation Directory**: `/GeometryGeneration/`

- `simpleCubicCell3DSphere.py`: Generates RVEs with periodic spherical inclusions (controlled via `periodicityFlags`).
- `randomInclusion3DSphere.py`: Generates RVEs with randomly positioned spheres, based on specified volume fraction and radius.
- `bodyCenteredInclusions.py`: Additional inclusion strategy for centered particles.

**Provided Geometries**:

- `s30.xdmf`: Legacy RVE with centered sphere (volume fraction = 0.3)
- `model-RI-VF-07-N_2-R_20.xdmf`: Random inclusion model with VF = 0.07, N = 2, R = 0.2
- `model-RI-VF-017-N_5-R_20.xdmf`
- `model-RI-VF-027-N_8-R_20.xdmf`

---

## 📄 Research Paper 2

**Title**: *Influence of Magnetic Boundary Conditions on the Quantitative Modelling of Magnetorheological Elastomers*  
**DOI**: [10.1016/j.mechmat.2023.104742](https://doi.org/10.1016/j.mechmat.2023.104742)  
**Code Repository**: [Zenodo Record](https://zenodo.org/records/8129310)

### Associated Files

**Main Simulation Script**:

- `fenics_magneto_hMRE_AirDomain.py`: Main FEniCSx simulation script using air domain for accurate boundary representation.

**Geometry Generation**:

- `GenerateCylinderWithAir.py`: Creates the cylinder with an embedded air domain.
- `GenerateXDMFH5.py`: Converts geometries into `.xdmf` and `.h5` formats for FEniCS compatibility.

**Provided Geometries**:

- `/model/cylinder_with_air_domain.xdmf`
- `/model/cylinder_with_air_domain_bulky.xdmf`

**Simulation Outputs**:

- `/results/AirDomain/cylinder_with_air_domain8.xdmf`
- `/results/AirDomain/cylinder_with_air_domain9.xdmf`

---

## 🧪 Citation

If you use this repository or its components in your research, please consider citing the original research papers linked above and referencing this repository.
