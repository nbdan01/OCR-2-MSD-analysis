# OCR-2-MSD-analysis
Calculate diffifusion coefficient and velocity from SPT trajectories


# Simulation
## Step 1: GEnerating the trajectories
- Open the script 'Generate_trajectories.mat'
- Choose the Folder to save the trajectories in the variable Foldersave

- Copy Paste your diffusion coefficient values (D) and velocity values (v), position (X) dependent in the variable Data_diffusion, and Data_antero_velocity and Data_velocity_retro respectively.
The first column of each variable corresponds to the position (X) and the second column corresponds to the D or v values.

D and V are smoothed with an sgolay filter of degree 3 to reduce the noise.

- Copy paste your ratio between directed and diffusive transport, location dependent in the variable Percentage_transport_antero.
- And indicate the corresponding positions in the variable Pos_percentage

- Do the same for the retrograde directed transport if the ratio is different (not the case here).




