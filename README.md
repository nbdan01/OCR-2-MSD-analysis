# OCR-2-MSD-analysis
Calculate diffifusion coefficient and velocity from SPT trajectories
Use the script `main_MSD_analysis.m`

# Simulation
## Step 1: Generate the trajectories
- Use the script `Generate_trajectories.m`
- Choose the Folder to save the trajectories in the variable `Foldersave`

- **Parameters**: Copy Paste your **diffusion coefficient** values (D) and **velocity** values (v), position (X) dependent in the variable `Data_diffusion`, and `Data_antero_velocity` and `Data_velocity_retro`respectively.
1st column: X position and 2nd column: D or v values.
- Choose your **OPTIONs**. You can either
  1) increase D by a factor `factorDCincrease` for ALL the positions
  2) increase D in the PCMC by a factor `factorDCincrease_PCMC` 
  3) set D = 0 for ALL the positions (`nulDiffusionCoefficient = 1`)

Note: D and v are smoothed with an sgolay filter of degree 3 to reduce the noise.
- Set the ratio IFT/Diffusion for each segements (dendrite, PCMC, TZ, PS, DS and Tip)

  > Standard OPTION  
  >
  > - dendrire_percentage = 1  
  > - PCMC_percentage = 0  
  > - TZ_percentage = 1  
  > - PS_percentage = 0.05  
  > - DS_percentage = 0.15  
  > - Tip_percentage = 0  
  
- Set the different **parameters** for the simulation

  > - lim_cilia = 7.5   _(Maximal size of the in-silico cilium, um)_  
  >   
  > - proba_antero = 1/3   _(Ratio between anterograde and retrograde IFT train frequencies)_  
  >  
  > - dt = 0.1   _(Time between 2 consecutives data points, s)_  
  >  
  > - Radius = 0.1   _(Radius of the in-silico cilium, um)_   
  >  
  > - WW =15   _(Size of the sliding window ( MSD) used to analyse experimental data, frame)_    
  >  
  > - maxtraj = 50   _(Number of trajectories per position)_    
  >  
  > - Nstep = 5000   _(Maximal number of steps)_  
  > 
  > - Precision_localization = 30/1000   _(Precision of localization, um)_    

Parameters D and V are saved in the subFolder `Parameters`
Generated trajectories are saved in the subFolder `Trajectories`

## Step 2: Generate the intensity profiles 

- Use the script `Intensity_profiles.m`
- Set the integration time in the variable `integration_over_time` to integrate/collect/gather the data points over time and generate the intensity profiles
Profiles are saved in the subFolder `Data hist` with the format `Dat_ + number of the first frame + .mat`

## Step 3: Manual renaming of the files
- In the subFolder `Data hist`, **Rename** the `name.mat` files with the format `name (1).mat`, `name (2).mat`, `name (3).mat`, ..., `name (50).mat`...

## Step 4: Select the steady-state
- Use the script `Steady_state_intensity_profiles.m`
- Set the **Parameters** of the profiles

> - min_X = -2  _(min position)_  
> - max_X = 8  _(max position)_  
> - binsize = 0.1 _(bin size, um)_  
> - smoothval = 5   -> 0.1*5 = 0.5 um in this case _(value used to smooth the histograms)_  

Profiles are saved in the subFolder `Movie`.
At the beginning, you will see that the profiles are correlated with the initial positions of the molecules.
- Set the steady state in the variables `time_beginning` and `time_end` as indicated in the profiles legend.
The average intensity profile during the steady state is saved in a csv file (position (um), mean I, std I)
