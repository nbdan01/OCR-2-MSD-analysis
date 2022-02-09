clearvars
close all

Folder_data='C:\Users\noemi\Documents\alpha classification analysis\Data\All_all.xlsx';
Folder_result='C:\Users\noemi\Documents\alpha classification analysis\Data\Results matlab\All_all\'

pxlSize = 0.053; % Pixel Size in um
dt = 0.1; % Time between 2 consecutive frames in second
min_number_point_in_traj = 14;
sct=0.1;
if exist(Folder_result,'dir')==0
mkdir(Folder_result)
else
end
%%
[PAR,ORT]=readdata(Folder_data); % Read xslx file and extract Parallel and orthogonal positions
checkparort(PAR,ORT); %Check if PAR and ORT have the right size

% Analysis of each molecule
for i_particule = 1 :size(PAR,2);
    %%
    close all
    
 
XC = PAR{i_particule}.*pxlSize; % X Positions in um
YC = ORT{i_particule}.*pxlSize; % Y Positions in um
TIME = [1:size(PAR{i_particule},1)]'; % Time in Frame

% Selection of trajectories with a number of positions above 15
if size(XC,1)<=min_number_point_in_traj
else

figure(1)
hold off
plot(XC,YC,'-','color',[0.8 0.8 1],'linewidth',1.5)
hold on
plot(XC,YC,'o','markerfacecolor','b','markeredgecolor',[0.8 0.8 1],'markersize',10,'linewidth',1.5)
set(gca,'fontsize',16)
xlabel('X \mum')
ylabel('Y \mum')
axis equal
saveas(1,[Folder_result,num2str(i_particule,'00%.0f'), '_Position_XY.fig'])
saveas(1,[Folder_result,num2str(i_particule,'00%.0f'), '_Position_XY.png'])
saveas(1,[Folder_result,num2str(i_particule,'00%.0f'), '_Position_XY.eps'])
%%
% Local analysis
% Initialization of the tables
out= [];
out15= [];
outalpha15= [];
outalphaT15 = [];
Tout = [];
Tout15 = [];
%%
ss = round(min_number_point_in_traj/2);
numberoflags = min_number_point_in_traj+1;

for it = ss+1:1:length(XC)-ss;
    Time_for_log_par = [];
    Dinst_for_log_par = [];
    Time_for_log_per = [];
    Dinst_for_log_per = [];
    [it,length(XC)]
    ID    = intersect(TIME', [it-ss:it+ss]);
    posit = [XC(ID),YC(ID)];                 % Positions in the window
% minimal number of consecutive points in the window to calculate the MSD
    tt = 3;                                  
    if length(ID) > tt;
        % Calculate MSD in the parallel direction
        [LagTpar,MSD_meanpar,MSD_stdpar,yfitpar,ppar] = MSD(posit(:,1),(TIME(ID)-min(TIME(ID))+1)*dt,numberoflags,dt);       
        
        % Calculate MSD in the perpendicular direction
        [LagTper,MSD_meanper,MSD_stdper,yfitper,pper] = MSD(posit(:,2),(TIME(ID)-min(TIME(ID))+1)*dt,numberoflags,dt);
        
        if size(MSD_meanpar,2)<=tt+1
        else
        
        % Fit the parallel MSD
        fittedDpar = polyfit(LagTpar(1:tt), MSD_meanpar(1:tt),1)
        D_par  =   fittedDpar(1)/2;
        Offset_par = fittedDpar(2);
        % Fit the perpendicular MSD        
        fittedDper = polyfit(LagTpar(1:tt), MSD_meanper(1:tt),1)
        D_per  =   fittedDper(1)/2;
        Offset_per = fittedDper(2);
        % Put fitted parameters in a table
        out = [out; [(XC(it)),(YC(it)),D_par D_per Offset_par Offset_per]]; 
        out15 = [out15; [posit(:,1),posit(:,2),D_par.*ones(size(posit,1),1) D_per.*ones(size(posit,1),1)] Offset_par.*ones(size(posit,1),1) Offset_per.*ones(size(posit,1),1)]; 
        % Corresponding Time in Frame
        Tout = [ Tout; it];
        Tout15 = [ Tout15;TIME(ID)];
        
% Calculate alpha per and par in this window
for iii = 1 : length(LagTper)
        D_inst_per_2   = MSD_meanper(iii)./(2.* LagTper(iii));
        Time_for_log_per = [Time_for_log_per ; iii];
        Dinst_for_log_per = [Dinst_for_log_per;D_inst_per_2];
        D_inst2par   = MSD_meanpar(iii)./(2.* LagTpar(iii));
        Time_for_log_par = [Time_for_log_par ; iii];
        Dinst_for_log_par = [Dinst_for_log_par;D_inst2par];
end

%%
    Log_of_LagTime_per = log(Time_for_log_per);
    Ratio_D_i_on_D2_inst_per = [Dinst_for_log_per]./Dinst_for_log_per(1);
    
    dLog_of_LagTime_per = Log_of_LagTime_per(2:end)-Log_of_LagTime_per(1:end-1);
    dLog_Ratio_D_i_on_D2_inst_per = log(Ratio_D_i_on_D2_inst_per(2:end))-log(Ratio_D_i_on_D2_inst_per(1:end-1));
    pente_per = 1+dLog_Ratio_D_i_on_D2_inst_per./dLog_of_LagTime_per;
    alpha_local_per = mean(pente_per(1:6));
    
    Log_of_LagTime_par = log(Time_for_log_par);
    Ratio_D_i_on_D2_inst_par = [Dinst_for_log_par]./Dinst_for_log_par(1);
   
    dLog_of_LagTime_par = Log_of_LagTime_par(2:end)-Log_of_LagTime_par(1:end-1);
    dLog_Ratio_D_i_on_D2_inst_par = log(Ratio_D_i_on_D2_inst_par(2:end))-log(Ratio_D_i_on_D2_inst_par(1:end-1));
    pente_par = 1+dLog_Ratio_D_i_on_D2_inst_par./dLog_of_LagTime_par;
    alpha_local_par = mean(pente_par(1:6));

    outalpha15 = [outalpha15; [posit(:,1),posit(:,2), alpha_local_par.*ones(size(posit,1),1) alpha_local_per.*ones(size(posit,1),1)]];
    outalphaT15 = [outalphaT15; TIME(ID) alpha_local_par.*ones(size(posit,1),1) alpha_local_per.*ones(size(posit,1),1)];
    
        end
    end
end
end
save([Folder_result,num2str(i_particule,'00%.0f'), '_Results.mat'])
end
