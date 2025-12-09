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

To integrate FABM with the model, you need to provide a model-specific FABM driver.
Create a new driver directory under FABM and copy the template driver into it:

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
The executable of the model will be created at:
```bash
build/release/bin/shelf_model
```
or inside the folder for the preset you have chosen (e.g. `build/debug/bin/shelf_model`).

## Repository Structure (Overview)

```
external/fabm/      # FABM submodule (biogeochemical models)
src/                # Model source code (physics, bio, IO, utilities)
cmake/Modules/      # Custom CMake find modules (FindNetCDF.cmake)
sims/               # Example configurations and forcing 
build/              # CMake build directory (created by you, not in git)
```

## Running the first simulation

Once the model is built, you just need a **run directory** containing:

- The main model configuration file (`main.yaml`): to specify details for location, time, grid, forcing, physics, output, etc. 

- Forcing data: one or more NetCDF files providing the atmospheric and surface vriables referenced in `main.yaml`.

- Optional: A YAML configuration file (`fabm.yaml`) required only if biogeochemistry is activated.

**Note:** If your `main.yaml` specifies that output should be written into a dedicated directory, make sure this directory exists before running the model. 
For example, inside `sims/simulationfolder` you can create it with:

```bash
mkdir output
```
This assumes your `output` block in the configuration file looks like this (adapt paths and filenames as needed):

```yaml
output:
    file: output/your_output_filename   # The model will write to output/your_output_filename.nc
    overwrite: yes
    ...
```

From the simulation directory, **run the `shelf_model` executable** adapting the path to your directory structure. 
For example, if your simulation case lives under `sims/simulationfolder/` within the repository, and you compiled the model in `build/release/`, you can launch it as:

```bash
cd sims/simulationfolder/

../../build/release/bin/shelf_model
```

Alternatively, you can configure the contents of `config_launch.json` to specify the path to the executable and then run the launcher script (Unix-like environments only):

```bash
bash launch_model.sh
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

| Variable                     | YAML key            | Units    | Required                     |
|------------------------------|---------------------|----------|------------------------------|
| Surface air temperature      | `surf_air_temp`     | °C       | Yes                          |
| Sea level pressure           | `sl_pressure`       | hPa      | Yes                          |
| Relative humidity            | `relative_humidity` | %        | Yes                          |
| Shortwave radiation          | `shortwave_radiation` | W m-2  | Yes                          |
| Longwave radiation           | `longwave_radiation`  | W m-2  | Yes                          |
| 10-m zonal wind              | `wind_u10`          | m s-1    | Yes                          |
| 10-m meridional wind         | `wind_v10`          | m s-1    | Yes                          |
| Precipitation                | `precipitation`     | m s-1    | Only if salinity is computed |
| Evaporation                  | `evaporation`       | m s-1    | Only if salinity is computed |
| Runoff                       | `runoff`            | m s-1    | Optional                     |


The model does not convert units: preprocessing scripts must provide the expected units.

### Tidal elliptic parameters

The model includes barotropic tidal forcing using tidal **elliptic parameters**.  
These describe the horizontal tidal current ellipse for each tidal constituent (M2, S2, N2, K1, O1) and consist of:

- **Semi-major axis** (cm s-1): maximum tidal current amplitude.  
- **Semi-minor axis** (cm s-1): minimum amplitude; may be negative to indicate rotation direction.  
- **Inclination** (degrees): orientation of the ellipse relative to the eastward direction.  
- **Greenwich phase angle** (degrees): timing of maximum tidal current relative to a universal reference.

The model can obtain these parameters in two ways:

1. **From a text file** (`tides.read_from_file: yes`)  
   The model selects the closest grid point using the latitude/longitude in `main.yaml` (within `tol_deg`).

2. **Manually in the YAML file** (`tides.read_from_file: no`)  
   You may specify the ellipse for each constituent directly under `ellipse_params`.

Missing or `nan` values are treated as zero.


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

## Quick guide to the structure of the code

At the highest level, the model is organised around a small number of core modules that talk to each other via the `shelfseas` orchestrator. 


### Main orchestrator

- **`shelfseas`** (in `src/shelfseas.F90`)  (The name will change once a new name has been decided)
  1. Scans forcing data, main parameters, builds the vertical grid (water, and optionally sediments).
  2. Initialises the forcing system, physics, optional biogeochemistry, and output.
  3. Runs the main time-stepping loop.

Think of `shelfseas` as the place to look if you want to understand “what happens in which order”.

### Main building blocks

- **Forcing**
  - `forcing_manager`:
    - scans forcing files,
    - returns a `ForcingSnapshot` (all surface fields at the needed time step).

- **Physics**
  - `physics_main`: the core physics module:
    - `init_physics`: allocate and initialise fields.
    - `solve_physics`: one full physics time step (turbulence, mixing, tidal forcing, etc.) called from `shelfseas`. Internally, it performs a **physics sub-cycling loop** to maintain numerical stability when the main time step `dt` is too large for the mixing or momentum equations.
    - `end_physics`: clean-up.
  - Additional internal modules (e.g. turbulence, mixing, tides) are called from here.
  - `physics_types`: holds `PhysicsState` (T, S, velocities, turbulence vars, etc.) and `PhysicsEnv`.

- **Biogeochemistry (FABM)**
  - `bio_main`: the main module for plugging biogeochemistry via FABM:
    - `init_bio_fabm`:  allocates tracers, links environment fields.
    - `integrate_bio_fabm`: advances all FABM tracers for one physics time step (with sub-cycling as needed).
    - `end_bio_fabm`: cleans everything associated with biogeochemistry.
  - Only used if `biogeochemistry.enabled: yes` in `main.yaml`.

- **Output**
  - `variable_registry` – metadata for each variable (name, units, dimensions, pointers to data).
  - `output_manager` – schedules output, accumulates statistics (mean/instant), and calls the writer.
  - `output_writer` – handles NetCDF file creation, dimensions, variables, and writing to disk.

### Sequence of steps in a nutshell

1. Read configuration in main.yaml.  
2. Build grids (water + optional sediments).
3. Scan Forcing and initialise the forcing manager.
4. For each main time step:
   - Update forcing snapshot.
   - Call `solve_physics`.
   - If enabled, call `integrate_bio_fabm`.
   - Call `om_step` (from `output_manager`) to accumulate and write output.
5. Finalise physics, FABM, forcing, and output.

If you’re new to the code start in `shelfseas.F90`, then look at the modules you're interested in physics, bio, io, etc. 
