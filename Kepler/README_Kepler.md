# Setting up RK4 Kepler

1) Pick which initial conditions you want
    - VEL: Uses initial time, position (radius and angle w.r.t x-axis), velocity (radial and angular), and a final time
    - ENERGY: Uses initial time, position (radius and angle w.r.t x-axis), energy and angular momentum (per unit mass), and a final time
        - NOTE: Some combinations of initial position and E and L will cause the integrator to exit due to unphysical orbital conditions...

2) Run the command: gcc rk4_Kepler.c -o rk4_Kepler.exe -lm
    - This creates the executable rk4_Kepler.exe

3) To run the executable for either .input file, run either the command: 
    - ./rk4_Kepler.exe < Kepler_ENERGY.input > Kepler_ENERGY_output.txt 
    - ./rk4_Kepler.exe < Kepler_VEL.input > Kepler_VEL_output.txt 

4) Run the command: gnuplot Kepler_plottingv2.gnu
    - In the set output, adjust the name of the file to match what you want plotted
    - In the plot command, adjust the name of the file to match what you want plotted