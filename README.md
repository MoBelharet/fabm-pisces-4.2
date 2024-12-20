# FABM-PISCES

This is a [FABM](https://fabm.net) port of [the PISCES model](https://doi.org/10.5194/gmd-8-2465-2015). It is based on the PISCES code that comes with the 4.2.2 version of [NEMO](https://www.nemo-ocean.eu/).

The code has been modularised to the point where it is straightforward to create configurations with any number
of phytoplankton and zooplankton types and any number of particulate organic matter classes - all by adjusting
the runtime configuration (`fabm.yaml`), no code change or recompilation needed. However, some more work would
be needed to fully support such configurations. Specifically, the zooplankton code would need to be changed to
handle a runtime-configurable number of prey types.

This new FABM-based version of PISCES (4.2) is implemented based on the previous FABM-based version (4.0) available in (https://github.com/BoldingBruggeman/fabm-pisces). In addition to updating all modules to be consistent with the original PISCES 4.2 release, this version includes new processes that were not in the previous release (fabm-pisces 4.0) such as:

* source of iron due to sea ice melt 
* iron input from hydrothermal vents
* iron source from sediment throughout the column 
* nutrient inputs from rivers


## How to build

This code must be compiled together with FABM. To do this, provide the following additional arguments to cmake [when you build FABM](https://github.com/fabm-model/fabm/wiki/Building-and-installing): `-DFABM_INSTITUTES=pisces -DFABM_PISCES_BASE=<PISCESDIR>`

Here, `<PISCESDIR>` is the directory with the FABM-PISCES code (the same directory that contains this readme file). Note that `-DFABM_INSTITUTES=pisces` will make FABM compile PISCES as the *only* available biogeochemical model. If you additionally want to have access to other biogeochemical models included with FABM, you can set `FABM_INSTITUTES` to a semi-colon separated list, e.g., `-DFABM_INSTITUTES="pisces;ersem"` (to prevent the shell from interpreting the semi-colons, you typically have to enclose this list with quotes).

For instance, to use PISCES with the latest stable release of the [General Ocean Turbulence Model (GOTM)](https://gotm.net/), do the following:

```
git clone --recurse-submodules -b v6.0 https://github.com/gotm-model/code.git gotm
git clone https://github.com/fabm-model/fabm.git
git clone https://github.com/BoldingBruggeman/fabm-pisces.git
mkdir build
cd build
cmake ../gotm -DFABM_BASE=../fabm -DFABM_INSTITUTES=pisces -DFABM_PISCES_BASE=../fabm-pisces
make install
```

This will install the GOTM executable with support for PISCES at `~/local/gotm/bin/gotm`.

## How to run a FABM-PISCES simulation

A `fabm.yaml` file with the PISCES configuration is provided under `<PISCESDIR>/testcases`. You can drop this file in the working directory of a FABM-compatible model such as GOTM to use it during simulation. Note that for GOTM, you will also need to ensure that `fabm/use` is set to `true` in `gotm.yaml`. Otherwise GOTM would run with physics only.
