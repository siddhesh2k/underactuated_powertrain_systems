# Underactuated Powertrain Systems — FCHEV SIMO Reference Tracking

Code and data accompanying the paper:

> S. Yadav and F. Assadian, "Practical Solutions to Underactuated SIMO
> Reference Tracking in Fuel-Cell Hybrid Electric Vehicles," ASME Journal of
> Dynamic Systems, Measurement, and Control (under review).

This repository reproduces the results of the paper: three control
architectures for tracking battery state-of-charge (SOC) and cumulative
hydrogen consumption with a single fuel-cell power actuator, on the UDDS drive
cycle. Results are provided for both the **linear plant** and the **nonlinear
Simulink powertrain model**.

## The problem

The fuel-cell hybrid powertrain has one actuator (fuel-cell power `Pfc`) but two
outputs to track (SOC and cumulative hydrogen mass `M_H2`), because the DC-bus
power balance `Pdem = Pfc + Pbat` makes battery power a dependent quantity. The
resulting plant is a single-input, two-output (SIMO) system whose transfer matrix
is rank-1. This repository contains the three architectures the paper develops
for tracking within that constraint, plus the benchmark.

## What's here

### Linear-plant results

| Script | Architecture | Paper figure |
|---|---|---|
| `src/Underactuated_MPC2.m` | Constrained linear MPC benchmark | Fig. 5 (`MPC_cons.png`) |
| `src/Underactuated_Optimizer15.m` | Sensor allocator + SISO H-infinity | Fig. 4 (`linear_SensorAlloc.png`) |
| `src/Underactuated_Null_allocator.m` | Joint H + G_a null-space prefilter | Fig. 3 (`linear_HGa.png`) |

Each linear script is run directly in MATLAB, loads the shared data file, runs
its closed-loop simulation on the UDDS cycle, prints the tracking metrics, and
writes its paper figure as a 300-DPI PNG.

### Nonlinear (Simulink) results

Each nonlinear architecture is a pair: a setup script that builds the controller
and workspace variables, and a Simulink model that embeds the nonlinear PEMFC,
battery (equivalent-circuit), and their low-level controllers.

| Setup script | Simulink model | Architecture | Paper figure |
|---|---|---|---|
| `src/nonlinear_HGa_setup.m` | `models/nonlinear_HGa.slx` | Joint H + G_a null-space prefilter | Fig. 6 (`nonlinear_HGa.png`) |
| `src/nonlinear_sensoralloc_setup.m` | `models/nonlinear_SensorAlloc.slx` | Sensor allocator + SISO H-infinity | Fig. 7 (`nonlinear_SensorAlloc.png`) |
| `src/nonlinear_MPC_setup.m` | `models/nonlinear_MPC.slx` | MPC benchmark | Fig. 8 (`nonlinear_MPC.png`) |

All three setup scripts first run the shared nonlinear plant parameter file
`src/PTsystem_Openloop_simulink_parameter_v03.m`, which defines the PMSM,
driveline, battery, and PEMFC parameters and builds the design objects used by
the models.

## Requirements

- MATLAB (developed on R2024a)
- Simulink
- Stateflow (the nonlinear models contain embedded charts)
- Control System Toolbox (`tf`, `ss`, `c2d`, `mixsyn`, `hinfsyn`, `augw`)
- Model Predictive Control Toolbox (for both MPC variants)
- Optimization Toolbox (`quadprog`, for the linear MPC benchmark)
- Symbolic Math Toolbox (for the null-space prefilter's Smith-McMillan step)
- Simscape / Powertrain Blockset (the nonlinear models use the equivalent-circuit
  battery block `autolibrcnetworksystem/Estimation Equivalent Circuit Battery`)

## How to run

### Linear results

From the repository root in MATLAB:

```matlab
addpath(genpath('.'));          % put src/, utils/, third_party/ on the path
cd src

Underactuated_MPC2            % MPC benchmark         -> MPC_cons.png
Underactuated_Optimizer15    % sensor allocator      -> linear_SensorAlloc.png
Underactuated_Null_allocator % joint H+G_a prefilter -> linear_HGa.png
```

The MPC and sensor-allocator scripts are self-contained. The joint H+G_a script
first runs `Hinf4_null_allocator.m` (the H-infinity synthesis and null-space
prefilter design), which depends on the helper functions in `utils/` and
`third_party/`.

### Nonlinear results

Each nonlinear result is produced in three steps. From the repository root in
MATLAB (`addpath(genpath('.'))` first), for each architecture:

1. Run the **setup script** (for example `nonlinear_HGa_setup.m`). This runs the
   parameter file, builds the controller, and populates the workspace variables
   the model's From-Workspace blocks read.
2. Run the matching **Simulink model** (for example `nonlinear_HGa.slx`). The run
   produces the logged signals in the `out` object.
3. Run the **"PAPER FIGURE from Simulink run"** section at the bottom of the same
   setup script. It reads `out`, prints the tracking metrics, and writes the
   300-DPI figure PNG.

The setup scripts save their controller variables across the parameter file's
`clear` so that the models can be run immediately after setup.

## Repository layout

```
.
├── src/          main result scripts (linear + nonlinear setup) and H-infinity synthesis
├── models/       nonlinear Simulink models
├── utils/        helper functions (all-pass generation, properness, Skogestad weight)
├── third_party/  third-party functions, original author headers retained
├── data/         shared UDDS target/demand data (see data/README.md)
└── figures/      output directory for generated PNGs
```

## Data

See `data/README.md` for the provenance of
`FCHEV_SDP_target_pdem_data_UDDS.mat`.

## License

Code in `src/`, `utils/`, and the Simulink models in `models/` is released under
the MIT License (see `LICENSE`). Files in `third_party/` retain their original
authors' licenses and headers and are redistributed under those terms.

## Citing

If you use this code, please cite the paper above.
