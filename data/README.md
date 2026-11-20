# Data

## `FCHEV_SDP_target_pdem_data_UDDS.mat`

This file contains the reference trajectories and driver power demand used by all
simulations in this repository, sampled over the 1370-second UDDS drive cycle at
1-second intervals.

### Contents

| Variable | Size | Description |
|---|---|---|
| `SOC_Target` | 2 × 1370 | Target battery state-of-charge trajectory (row 2 used) |
| `H2_cons_target` | 2 × 1370 | Target cumulative hydrogen consumption, grams (row 2 used) |
| `P_dem_simulink` | 2 × 1370 | Driver power demand, watts (row 2 used) |

In each variable, row 1 is the time base and row 2 is the signal.

### Provenance

The reference trajectories and demand profile were extracted from a
high-fidelity Simulink model of the fuel-cell hybrid vehicle. That model was
developed in conjunction with Hyundai for a prior, related project. The model was
built by the author and validated against vehicle behavior together with a
Hyundai engineer. The target SOC and hydrogen-consumption trajectories were then
produced by an upstream energy-management optimization run on that model for the
UDDS cycle, and serve here as the references the downstream tracking controllers
must follow.

The data is released with permission for the purpose of reproducing the results
in the accompanying paper.
