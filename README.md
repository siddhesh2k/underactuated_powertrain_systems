# Underactuated Powertrain Systems — FCHEV SIMO Reference Tracking

Code and data accompanying the paper:

> S. Yadav and F. Assadian, "Practical Solutions to Underactuated SIMO
> Reference Tracking in Fuel-Cell Hybrid Electric Vehicles," ASME Journal of
> Dynamic Systems, Measurement, and Control (under review).

This repository reproduces the **linear-plant results** of the paper: the three
control architectures for tracking battery state-of-charge (SOC) and cumulative
hydrogen consumption with a single fuel-cell power actuator, on the UDDS drive
cycle.

## The problem

The fuel-cell hybrid powertrain has one actuator (fuel-cell power `Pfc`) but two
outputs to track (SOC and cumulative hydrogen mass `M_H2`), because the DC-bus
power balance `Pdem = Pfc + Pbat` makes battery power a dependent quantity. The
resulting plant is a single-input, two-output (SIMO) system whose transfer matrix
is rank-1. This repository contains the three architectures the paper develops
for tracking within that constraint, plus the benchmark.

## What's here

| Script | Architecture | Paper figure |
|---|---|---|
| `src/Underactuated_MPC2.m` | Constrained linear MPC benchmark | Fig. 5 (`MPC_cons.png`) |
| `src/Underactuated_Optimizer15.m` | Sensor allocator + SISO H-infinity | Fig. 4 (`linear_SensorAlloc.png`) |
| `src/Underactuated_Null_allocator.m` | Joint H + G_a null-space prefilter | Fig. 3 (`linear_HGa.png`) |

Each script loads the shared data file, runs its closed-loop simulation on the
UDDS cycle, prints the tracking metrics (SOC RMS, M_H2 RMS, saturation, etc.)
reported in the paper's linear-results tables, and writes its paper figure as a
300-DPI PNG.

## Requirements

- MATLAB (developed on R2024a)
- Control System Toolbox (`tf`, `ss`, `c2d`, `mixsyn`, `hinfsyn`, `augw`)
- Optimization Toolbox (`quadprog`, for the MPC benchmark)
- Symbolic Math Toolbox (for the null-space prefilter's Smith-McMillan step)

## How to run

From the repository root in MATLAB:

```matlab
addpath(genpath('.'));          % put src/, utils/, third_party/ on the path
cd src

Underactuated_MPC2            % MPC benchmark        -> MPC_cons.png
Underactuated_Optimizer15    % sensor allocator     -> linear_SensorAlloc.png
Underactuated_Null_allocator % joint H+G_a prefilter -> linear_HGa.png
```

The MPC and sensor-allocator scripts are self-contained. The joint H+G_a script
first runs `Hinf4_null_allocator.m` (the H-infinity synthesis and null-space
prefilter design), which depends on the helper functions in `utils/` and
`third_party/`.

## Repository layout

```
.
├── src/          main result scripts + H-infinity synthesis
├── utils/        helper functions (all-pass generation, properness, Skogestad weight)
├── third_party/  third-party functions, original author headers retained
├── data/         shared UDDS target/demand data (see data/README.md)
└── figures/      output directory for generated PNGs
```

## Data

See `data/README.md` for the provenance of
`FCHEV_SDP_target_pdem_data_UDDS.mat`.

## License

Code in `src/` and `utils/` is released under the MIT License (see `LICENSE`).
Files in `third_party/` retain their original authors' licenses and headers and
are redistributed under those terms.

## Citing

If you use this code, please cite the paper above.
