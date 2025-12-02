# S2P3-new - 1-D Shelf Sea Physics & Biogeochemistry

This testing version of the model (still in development) builds on the S2P3 physics core and provides a flexible interface to biogeochemical modules through FABM. It also incorporates concepts from BROM-transport to represent sediment transport processes, enabling coupled water-sediment biogeochemistry.

## Supported Platforms

The model is **developed and tested on Linux**.

CMake already contains compiler flags for GNU, Intel, and Cray Fortran compilers, so it *should* work on most UNIX-like systems:

- **Linux** (recommended, tested)
- **macOS** (via Homebrew/MacPorts) - not yet tested

---

## Required Software

To build the model you need:

### 1. Git

To clone the repository and initialise the FABM submodule.

Check by typing in the console:

```bash
git --version
```

### 2. CMake ≥ 3.23

Check by typing in the console:

```bash
cmake --version
```

### 3. Fortran & C compilers
On Linux the recommended compilers are:
- gfortran (version>=10)
- gcc

### 4. NetCDF Libraries (C + Fortran)
On Debian/Ubuntu, the required packages are typically:
- `libnetcdf-dev` and `libnetcdff-dev`

Check with:
```bash
nf-config --all
```

## Installing Dependencies on Linux

### Debian/Ubuntu
```bash 
sudo apt update
sudo apt install -y \
    git cmake gfortran gcc \
    libnetcdf-dev libnetcdff-dev
```

## Cloning the repository (including FABM)
This repository includes FABM as a git submodule under `external/fabm`.

```bash
git clone https://github.com/BeatrizArellano/model_indev model
cd model

# Fetch FABM (required)
git submodule update --init --recursive -- external/fabm
```

## Copying the FABM driver specific to the model

Create a new directory called `shelf_model` (for now) inside `external/fabm/src/drivers/` and
copy `fabm_driver.h` inside `external/driver/` into `external/fabm/src/drivers/shelf_model/`.
```bash
mkdir external/fabm/src/drivers/shelf_model
cp external/driver/fabm_driver.h external/fabm/src/drivers/shelf_model/
```

## Building the model
From the repository root, run (Ignoring the Warnings):

```bash
cmake --preset release
cmake --build --preset release
```
The executable of the model will be creasted at:
```bash
build/release/shelf_model
```

## Running the first simulation
Once the model is built, you need:

- A run directory containing the following

- A YAML configuration file (main.yaml)

- Forcing data in NetCDF format

- Optional: A YAML configuration file if including biogeochemistry (fabm.yaml)


From the simulation directory, run the shelf_model executable correcting for the path in your directory structure:

```bash
../../build/bin/shelf_model
```

## Repository Structure (Overview)

```
external/fabm/      # FABM submodule (biogeochemical models)
src/                # Model source code (physics, bio, IO, utilities)
cmake/Modules/      # Custom CMake find modules (FindNetCDF.cmake)
sims/               # Example configurations and forcing 
build/              # CMake build directory (created by you, not in git)
```

