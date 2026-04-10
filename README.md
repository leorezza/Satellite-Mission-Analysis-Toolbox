# Satellite Mission Analysis & Constellation Simulator

A MATLAB-based tool for preliminary mission analysis of satellites in Low Earth Orbit (LEO). This project allows users to simulate Sun-Synchronous Orbits (SSO), visualize Ground Tracks, and calculate visibility windows for single satellites or multi-satellite constellations.

## Features

1. **Sun-Synchronous Orbit Modeling:** Automatically calculates the required inclination for SSO based on altitude.
2. **Constellation Analysis:** Supports n-satellite constellations with customizable RAAN and Mean Anomaly.
3. **Ground Track Visualization:** Plots the satellite's path over Earth (or Moon) maps.
4. **Visibility & Access Windows:** Calculates the precise time (hours, minutes, seconds) a satellite or constellation remains in view of a specific ground target.
5. **Propagator:** Includes a restricted two-body problem numerical integrator (ODE113).
