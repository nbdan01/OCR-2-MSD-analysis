clearvars
close all
Folderdatadata = 'D:\Modeling OCR2\Results\Trajectories proba antero retro 0p3 5000 steps\';
listmatfilesheredtat = dir([Folderdatadata,'*mat'])
%%
Foldersavedat = [Folderdatadata,'Data hist\'];
mkdir(Foldersavedat)
%%
load([Folderdatadata,listmatfilesheredtat(1).name]);
% Collect all the individual position from all the molecules within 
% a laps of 3 seconds
frame_num = round(3./dt)
for ig = 1111:frame_num:Nstep-frame_num
    Pos_time = [];
    numboftraj = 0;
    for ithad = 1: size(listmatfilesheredtat,1)
    
        load([Folderdatadata,listmatfilesheredtat(ithad).name]);
        
        [ig Nstep-frame_num ithad size(listmatfilesheredtat,1)]
	
        if size(PositionX,1)<ig+frame_num-1

        else
            numboftraj = 1+numboftraj;
            Pos_time = [Pos_time;PositionX(ig:ig+frame_num-1,1)];
        end
    end

%     pause(0.15)
    save([Foldersavedat,'Dat_' num2str(ig+15,'00%.0f'),'.mat'])
end
