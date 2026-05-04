# NACA Wind-Tunnel Simulator

MATLAB coursework project that simulates inviscid flow around NACA 4-digit aerofoils with a constant-strength doublet panel method.

## Features

- Generates NACA 4-digit aerofoil geometry with cosine-spaced panels.
- Solves for doublet strengths using the Kutta condition.
- Reports coefficient of lift.
- Produces streamline and velocity-arrow plots.
- Includes an automated NACA 2412 study against the bundled XFOIL polar data.

## Files

- `MATLAB.m` - main script.
- `panelgen.m` - NACA geometry and panel generation.
- `cdoublet.m` - velocity induced by a unit-strength doublet panel.
- `xf-naca2412-il-1000000.txt` - XFOIL comparison data for the NACA 2412 study.

## Running

Open MATLAB in this folder and run:

```matlab
MATLAB
```

Entering `2412` runs the preset comparison study with `U = 15 m/s`, panel counts of `50`, `100`, and `200`, and angles of attack from `0` to `10` degrees. Other valid NACA 4-digit codes prompt for freestream velocity, angle of attack, and panel count.

Generated plots are saved as PNG files in the project folder.
