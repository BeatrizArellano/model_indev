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

## Repository Structure (Overview)

```
external/fabm/      # FABM submodule (biogeochemical models)
src/                # Model source code (physics, bio, IO, utilities)
cmake/Modules/      # Custom CMake find modules (FindNetCDF.cmake)
sims/               # Example configurations and forcing 
build/              # CMake build directory (created by you, not in git)
```

## Running the first simulation
Once the model is built, you need:

- A run directory containing the following

- A YAML configuration file (main.yaml)

- Forcing data in NetCDF format

- Optional: A YAML configuration file if including biogeochemistry (fabm.yaml)


From the simulation directory, run the shelf_model executable adapting the path to your directory structure. This is an example considering that your simulation case is inside sims/simulationfolder/ within the repository structure. 

```bash
../../build/release/bin/shelf_model
```

## Quick guide to `main.yaml` (model configuration)

The parameters relevant to the simulation can be edited in the configuration file `main.yaml`: site, time period, grid, forcing, physics, tides, and (optionally) FABM biogeochemistry.  
You normally only need to edit a few sections for a new site / experiment.

---

### 1. YAML conventions

- **Numbers**  
  - Integers or reals are both fine.  
  - Some fields are validated (e.g. must be positive, within min/max ranges).

- **Booleans**  
  - Prefer: `true`, `false`.  
  - Also accepted: `yes/no`, `on/off`.

- **Strings**  
  - **Do not** use quotes.  
  - Example:  
    `file: path/to/file.nc`

- **Null / “unset” values**  
  - Use: `null`, `NULL`, `Null`, `none`, `None`, `~`  
  - Or simply comment a key out.  
  - `null` / missing are treated the same.  
  - **Do not** use `NaN` to mean “unset”.

- **Comments & indentation**  
  - `#` starts a comment (rest of the line is ignored).  
  - Use **spaces only** for indentation (recommended: 4 spaces).  
  - Avoid TABs – they can break parsing.

---

## Forcing

The model reads atmospheric and surface forcing from NetCDF files. A forcing file:

- Contains a time dimension with a CF-compatible time variable.

- Can contain information for different sites referenced by latitude and longitude. The model will choose the closest grid-point to the specified location. 

- Must include all required variables for physics (and optionally biogeochemistry). Variables may be in one file or separate files, according to the filename field in `main.yaml`.

### Forcing Variables:

Below are the variables the model may read, depending on your YAML configuration.
Each variable follows the pattern:
```yaml
mode: file | constant | off
name: <variable in NetCDF>
filename: off | <custom file>
```

- Surface air temperature (surf_air_temp): degrees Celsius.
- Sea level pressure (sl_pressure): hPa.
- Relative humidity (relative_humidity): percentage (%).
- Shortwave radiation (shortwave_radiation): W/m2.
- Longwave radiation (longwave_radiation): W/m2.
- 10-m zonal wind (wind_u10): m/s.
- 10-m meridional wind (wind_v10): m/s.
- Precipitation (precipitation): m/s. (Optional, but required when computing salinity)
- Evaporation (evaporation): m/s. Negative values mean evaporation. (Optional, but required when computing salinity)
- Runoff (runoff): m/s. (Optional)

The model does not convert units: preprocessing scripts must provide the expected units.

## Biogeochemistry

Biogeochemistry is configured separately using a FABM configuration file (typically fabm.yaml).
This file defines:

- Which FABM modules are active (e.g., NPZD, ERSEM, BROM components, custom models).
- The parameters for each module (rate constants, remineralisation rates, sinking speeds, light parameters, etc.).

The main model will:

1. Initialise FABM  and provide all environmental fields
2. Integrate the biogeochemical tracers consistently with the physics time-stepping.
3. Transport tracers according to vertical diffusion and their specified sinking/floating speeds. 

The biogeochemistry behaviour is entirely determined by the contents of its configuration file (fabm.yaml), allowing flexible switching between simple test models and full ecosystem models. A sample configuration file can be found inside `sims/Biogeochemistry_toy_case/fabm.yaml`.