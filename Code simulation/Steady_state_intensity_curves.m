clearvars
close all
% FFolderdatanoIFT = 'E:\Modeling OCR2\OCR2\OCR2 modeling\2021_05_06 Results no IFT\Trajectories proba antero retro 0p3 5000 steps 6\copy data hist\';
% FFolderdatanoIFT = 'E:\Modeling OCR2\OCR2\OCR2 modeling\2021_05_06 Results TZ only IFT\Trajectories proba antero retro 0p3 5000 steps 6\';

% FFolderdataIFT = 'E:\Modeling OCR2\OCR2\OCR2 modeling\2021_03_30 Results\Trajectories proba antero retro 0p3 5000 steps 6\Compilation\Data\';
% FFolderdataIFT = 'E:\Modeling OCR2\OCR2\OCR2 modeling\2021_03_30 Results diff PCMC x2\Trajectories proba antero retro 0p3 5000 steps 6\Data hist\';
FFolderdataIFT = 'E:\Modeling OCR2\OCR2\OCR2 modeling\2021_03_30 Results diff cst 0\Trajectories proba antero retro 0p3 5000 steps 6\Data hist\';
FFolderdatanoIFT = 'E:\Modeling OCR2\OCR2\OCR2 modeling\2021_03_30 Results no diffusion\Trajectories proba antero retro 0p3 5000 steps 6\Data hist\';
% FFolderdatanoIFT = 'E:\Modeling OCR2\OCR2\OCR2 modeling\2021_03_30 Results DC all factor 3 V2\Trajectories proba antero retro 0p3 5000 steps 6\Data hist\';

% FFolderdatanoIFT = 'E:\Modeling OCR2\OCR2\OCR2 modeling\2021_05_06 Results no IFT without automatic entry of mol\Trajectories proba antero retro 0p3 5000 steps 6\Data hist\';
% FFolderdataIFT = 'E:\Modeling OCR2\OCR2\OCR2 modeling\2021_05_06 Results without automatic entry of mol\Trajectories proba antero retro 0p3 5000 steps 6\Data hist\';

Fffoldersaveall = 'E:\Modeling OCR2\OCR2\OCR2 modeling\Compilation all smooth histo\';
%%
% FFolderdataIFT = 'D:\Modeling OCR2\Results\Trajectories proba antero retro 0p3 5000 steps\Data hist\';

%%
clearvars
close all
% Folder data
FFolderdataIFT = 'D:\Modeling OCR2\Results\Trajectories proba antero retro 0p3 5000 steps\Data hist\';

% Foldersave
Fffoldersaveall = [FFolderdataIFT, 'Intensity curves\'];
mkdir(Fffoldersaveall)
% Additional subfolder name if needed
FoldersavedatIFT = [Fffoldersaveall,'entry\IFT\'];
mkdir(FoldersavedatIFT)

FoldersavedatIFThist = [FoldersavedatIFT,'Movie Histograms\'];
mkdir(FoldersavedatIFThist)
% list mat files
llllistmatfileIFT = dir([FFolderdataIFT,'*mat']);

%%
% Parameters of the histograms
%xlim
max_X = 8; 
min_X = -2;
% bin size
binsize = 0.1;
% value used to smooth the histograms
smoothval = 5; % ==> 0.1*5 = 0.5 um in this case

DataforIFT = [];
for itall = [1:size(llllistmatfileIFT,1)]
% load data 
load([FFolderdataIFT,llllistmatfileIFT(itall).name]);
% extract filename
fffilenamehere = llllistmatfileIFT(itall).name
% extract the number of the histogram
ii2 = find(fffilenamehere=='(');
ii3 = find(fffilenamehere==')');
numberhistIFT = fffilenamehere(ii2+1:ii3-1)


Pos_timeIFT = Position_on_this_laps_time;%Pos_time;

clear Position_on_this_laps_time Pos_time

% store the histograms with time, x position, and occurrence
DataforIFT = [DataforIFT;str2num(numberhistIFT).*ones(size(smooth(nn,smoothval),1),1) smooth(nn,smoothval) smooth(xx,smoothval)./sum(xx)];

% plot the histograms smoothed and save them in a folder
% this part is used to visualize and select curves after reaching the steady-state
% to do this, drag the folder Movie histogram in ImageJ or just look manually in the
% Folder
figure(11)
hold off
[xx2,nn2] = hist(Pos_timeIFT,[min_X:binsize:max_X]);
plot(smooth(nn2,smoothval),smooth(xx2,smoothval),'color', [0.7 0.7 0.7], 'linewidth',1.5)
hold on
xlabel('Position \mum')
ylabel('Occurrence')
set(gca,'fontsize',14)
ylim([0 800])
legend(['Hist # ', numberhistIFT])
saveas(11,[FoldersavedatIFThist,'Histno_00' ,numberhistIFT,'.png'])
pause(0.15)
end


%%
% Sorting the profiles by time
Time = DataforIFT(:,1);
[~,order] = sort(Time);
Time_order = Time(order);
DataforIFT_order = DataforIFT(order,:);

%%
% Select the steady state

% time_beginning = 50;
% time_end = 90;
time_beginning =60;
time_end = 130;

DatnoX = [];
DatnoY = [];
DatX = [];
DatY = [];
% overlap all the profiles selected
% check in you are in the steady-state
figure(100)
hold off
for itime = time_beginning:time_end
    
    plot(DataforIFT_order(DataforIFT_order(:,1)==itime,2),DataforIFT_order(DataforIFT_order(:,1)==itime,3),'color',[0.7 0.7 0.7],'linewidth',1.5)
    hold on
    % Store data to calculate the mean and std
    DatX = [DatX;DataforIFT_order(DataforIFT_order(:,1)==itime,2)'];
    DatY = [DatY;DataforIFT_order(DataforIFT_order(:,1)==itime,3)'];

end
hold on
% plot the average profile
plot(mean(DatX),mean(DatY),'k','linewidth',1.5)
xlabel('Position \mum')
ylabel('Intensity a.u.')
set(gca,'fontsize',18)
saveas(100,[FoldersavedatIFT,'Avegraged_intensity_D0.png'])
saveas(100,[FoldersavedatIFT,'Avegraged_intensity_D0.fig'])
% save average profile in xls
Average_profile = [mean(DatX)' mean(DatY)' std(DatY)'];
writematrix(Average_profile,[FoldersavedatIFT,'Average_profile.csv']) 
