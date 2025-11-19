# Energy-Dissipation-Limits-of-Circuit-Systems

This repository contains the source code used to generate the numerical results and figures in the paper:

> **"Energy Dissipation Limits of Circuit Systems"**  
> Y. L. Dong *et al.*, Scientific Reports, 2025.  
> (DOI: to be updated)

The main purpose of this repository is to facilitate reproducibility and further research on the thermodynamic limits of energy dissipation in circuit systems.

---

## Repository structure

```text
Energy-Dissipation-Limits-of-Circuit-Systems/
├── BSC.m                 # Variation of energy dissipation and information energy efficiency of a binary symmetric channel w.r.t. input probability and prior error probability
├── BSC_I_Q.m             # Computes the change in information energy consumption
├── BSC_eta_improved.m    # Computes the improved information energy efficiency
├── energy_optimization.m # Optimization of total energy consumption
├── proc_tran.m           # Handles the trade-off of transmission energy consumption
├── relative_entropy.m    # Computes the relative entropy between probability distributions
├── total_energy_opt_time.m  # Impact of time allocation on total energy consumption
├── tran_parallel.m       # Energy consumption of parallel transmission
└── README.md             # This file

“Requirements / How to run”：
```markdown
## Requirements
- MATLAB (R20xxa or later is recommended)
## How to run

Each `.m` file can be executed in MATLAB to reproduce the corresponding numerical results:
- `BSC.m`: sweep input probability / error probability and plot energy dissipation and information energy efficiency.
- `energy_optimization.m`: reproduce the total energy optimization results in the paper.
- `tran_parallel.m`: reproduce the energy consumption of parallel transmission, etc.
