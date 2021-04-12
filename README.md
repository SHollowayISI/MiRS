# MiRS
MATLAB project for MiRS

## Overview

This project simulates the MiRS radar scenario, performs signal and data processing, and displays results. The multilateration solver used to determine position and phase offsets is included in the "MiRS Solver" folder.

### Radar Simulation
* Full simulation is run from "FullSystem.m" or "FullSystem_Test.m"
* Settings, data, and results are passed through "scenario" object.
* Parameters are set through three Setup files:
  * SetupRadarScenario: Waveform and processing parameters
  * SetupSimulation: Simulation parameters
  * SetupTarget: Target parameters
* Use "FullSystem_Test.m" to run Monte Carlo simulation

### Position Solver
* Levenberg-Marquardt Solver can be run through "MiRS Solver/MiRS_Solver_3D.m" or "MiRS Solver/MiRS_Solver_3D_RunScripts.m"
* This file also includes phase shift calculation, gain calculation, and heatmap generation.
* Use "RunScripts.m" to run Monte Carlo simulations of solver
