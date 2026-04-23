# Setting up RK4 Kepler

1) Pick which initial conditions you want (either `Kepler_VEL.input` or `Kepler_ENERGY.input`):
    - VEL: Uses initial time, position (radius and angle w.r.t x-axis), velocity (radial and angular), and a final time
    - ENERGY: Uses initial time, position (radius and angle w.r.t x-axis), energy and angular momentum (per unit mass), and a final time
        - NOTE: Some combinations of initial position and E and L will cause the integrator to exit due to unphysical orbital conditions (e.g., $v_r^2 < 0$)...

2) Build executable: 
    `gcc rk4_Kepler.c -o rk4_Kepler.exe -lm`
    - Builds the executable `rk4_Kepler.exe`

3) Run the executable for either .input file, run either the command: 
    `./rk4_Kepler.exe < Kepler_ENERGY.input > Kepler_ENERGY_output.txt `
    - or 
    `./rk4_Kepler.exe < Kepler_VEL.input > Kepler_VEL_output.txt `

4) Plotting data: 
    `gnuplot Kepler_plottingv2.gnu`

    - In the set output, adjust the name of the file to match what you want plotted
    - In the plot command, adjust the name of the file to match what you want plotted

    