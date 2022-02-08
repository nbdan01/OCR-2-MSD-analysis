# OCR-2-MSD-analysis
Calculate diffifusion coefficient and velocity from SPT trajectories


# Simulation
## Step 1: Generate the trajectories
- Open the script `Generate_trajectories.mat`
- Choose the Folder to save the trajectories in the variable `Foldersave`

- Copy Paste your diffusion coefficient values (D) and velocity values (v), position (X) dependent in the variable `Data_diffusion`, and `Data_antero_velocity` and `Data_velocity_retro respectively`.
The first column of each variable corresponds to the position (X) and the second column corresponds to the D or v values.

D and V are smoothed with an sgolay filter of degree 3 to reduce the noise.

- Copy paste your ratio between directed and diffusive transport, location dependent in the variable `Percentage_transport_antero`.
- And indicate the corresponding positions in the variable `Pos_percentage`

- Do the same for the retrograde directed transport if the ratio is different (not the case here).

- Set the different parameters

```
% Maximal size of the in-silico cilium
lim_cilia = 7.5;
% Ratio between anterograde and retrograde IFT train frequencies
proba_antero = 1/3;
% Time between 2 consecutives data points
dt = 0.1;
Diff_on = [1 0];
% Radius of the in-silico cilium
Radius = 0.1;
% Size of the sliding window ( MSD) used to analyse experimental data 
WW =15;
% Maximal number of steps
Nstep = 5000;
% Number of trajectories per position
maxtraj = 50;
% Precision of localization
Precision_localization = 30/1000;
```
Parameters D and V are saved in the subFolder `Parameters`
Generated trajectories are saved in the subFolder `Trajectories`

## Step 2: Create the time-dependent intensity/number of molecules profiles 

- Open the script `Intensity_profiles.mat`
- Set the integration time in the variable `frame_num` to integrate/collect/gather the data points over time and make the profiles
The profiles are saved in the subFolder `Data hist` with the format `Dat_ + number of the first frame + .mat`

## Step 3: Manual sorting of the profiles
- Sort the `.mat` in the subFolder `Data hist` by Creation Date and rename the files using the format `name (1)`, `name (2)`, `name (3)`, ..., `name (50)`...

## Step 4: 


