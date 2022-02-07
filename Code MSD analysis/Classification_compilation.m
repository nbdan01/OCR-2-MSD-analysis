clearvars
close all
 
FolderResultsG ='C:\Users\noemi\Documents\alpha classification analysis\Data\Results matlab\Results all\'
FolderSavehere = [FolderResultsG,'2020_12_17 Compilation\'];
mkdir(FolderSavehere)
listmatfiles = dir([FolderResultsG,'*mat'])
warning off
DGLOBALE= [];
VELOCITYantero=[];
VELOCITYretro=[];
ALPHAGLOBALE15=[];
RESULTSEND = [];
Pre_locGper = [];
Pre_locGpar = [];

%%
for itGlobal=1: size(listmatfiles,1)
    clear XC
    load([FolderResultsG,listmatfiles(itGlobal).name])

        close all
        % Time X Y Dpar Dper Offsetpar Offsetper
DGLOBALE = [DGLOBALE ; Tout15(:,1) out15(:,1) out15(:,2) out15(:,3) out15(:,4) out15(:,5) out15(:,6)];
% Time X Y alphapar alphaper LagTime MSD
ALPHAGLOBALE15 = [ALPHAGLOBALE15; outalphaT15(:,1) outalpha15(:,1) outalpha15(:,2) outalpha15(:,3) outalpha15(:,4) outalpha15(:,[5:end])];


%%
% Positions X Y and alpha par, per values
% Re-allocation of the alpha values
XXC = outalpha15(:,1);
YYC = outalpha15(:,2);
AAC = outalpha15(:,3);
AACp = outalpha15(:,4);
ToutGG = outalphaT15(:,1);

% Sort by time
[~,sortTC] = sort(ToutGG)
XXC = outalpha15(sortTC,1);
YYC = outalpha15(sortTC,2);
AAC = outalpha15(sortTC,3);
AACp = outalpha15(sortTC,4);
ToutGG = outalphaT15(sortTC,1);

% Eliminate nan values in positions
XXCb = XXC(isnan(XXC)==0 & isnan(YYC)==0 );
YYCb =YYC(isnan(XXC)==0 & isnan(YYC)==0 );
AACb =AAC(isnan(XXC)==0 & isnan(YYC)==0 );
AACpb =AACp(isnan(XXC)==0 & isnan(YYC)==0 );
ToutGGb =ToutGG(isnan(XXC)==0 & isnan(YYC)==0 );

% Initialization of the tables
XXCg = [];
YYCg = [];
AACg = [];
AACpg = [];
ToutGGg = [];

for itt = min(ToutGGb):max(ToutGGb)
%     find(ToutGGb==itt)
   XXCg = [XXCg;mean(XXCb(ToutGGb==itt))];
   YYCg = [YYCg;mean(YYCb(ToutGGb==itt))];
   AACg = [AACg;mean(AACb(ToutGGb==itt))];
   AACpg = [AACpg;mean(AACpb(ToutGGb==itt))]; 
   ToutGGg = [ToutGGg;mean(ToutGGb(ToutGGb==itt))];
end


%New positins and their alpha values
Xxhereconv = XXCg(isnan(XXCg)==0 & isnan(YYCg)==0 );
Yyhereconv =YYCg(isnan(XXCg)==0 & isnan(YYCg)==0 );
AAhereconv =AACg(isnan(XXCg)==0 & isnan(YYCg)==0 );
AAhereconvp =AACpg(isnan(XXCg)==0 & isnan(YYCg)==0 );
TThereconvp =ToutGGg(isnan(XXCg)==0 & isnan(YYCg)==0 );

% Position when there is transport alphapar > 1.4
threshold = 1.4;
Xxhereconvtt = Xxhereconv(AAhereconv>=threshold); % Selection de position when alphaparallel up to 1.3

% Calulate the diffusion coefficient everywhere
% 'MSD' per for three time lag
MSDper1 = (Yyhereconv(2:end)-Yyhereconv(1:end-1)).^2;
MSDper2 = (Yyhereconv(3:end)-Yyhereconv(1:end-2)).^2;
MSDper3 = (Yyhereconv(4:end)-Yyhereconv(1:end-3)).^2;
DDlper = [MSDper1(1:end-2)  MSDper2(1:end-1) MSDper3];
            
D_Holcman_fitper = [];
Pre_locper =[];

for eachdiff = 1:size(DDlper,1)
	p = polyfit([1 2 3].*dt,DDlper(eachdiff,:),1);
    Pre_locper = [Pre_locper;(p(2))];
    Pre_locGper = [Pre_locGper;(p(2))];
    D_Holcman_fitper = [D_Holcman_fitper;p(1)/2];
end

% 'MSD' par for three time lag
MSDpar1 = (Xxhereconv(2:end)-Xxhereconv(1:end-1)).^2;
MSDpar2 = (Xxhereconv(3:end)-Xxhereconv(1:end-2)).^2;
MSDpar3 = (Xxhereconv(4:end)-Xxhereconv(1:end-3)).^2;
DDlpar = [MSDpar1(1:end-2)  MSDpar2(1:end-1) MSDpar3];

D_Holcman_fitpar = [];
Pre_locpar =[];

for eachdiff = 1:size(DDlpar,1)
    p = polyfit([1 2 3].*dt,DDlpar(eachdiff,:),1);
    Pre_locpar = [Pre_locpar;(p(2))];
    Pre_locGpar = [Pre_locGpar;(p(2))];
    D_Holcman_fitpar = [D_Holcman_fitpar;p(1)/2];
end


%%
% RESULTSEND = [RESULTSEND;Xxhereconv Yyhereconv AAhereconv AAhereconvp [D_Holcman_fitpar ; D_Holcman_fitpar(end).*ones(3,1)] [D_Holcman_fitper ; D_Holcman_fitper(end).*ones(3,1)]]
RESULTSEND = [RESULTSEND;Xxhereconv Yyhereconv AAhereconv AAhereconvp [MSDpar1./(2.*dt) ; MSDpar1(end)./(2.*dt).*ones(1,1)]  [MSDper1./(2.*dt) ; MSDper1(end)./(2.*dt).*ones(1,1)] ]

%%
% Select transport and classify in anterograde and retrograde transport
if isempty(Xxhereconvtt)==0
    TThereconvtt = TThereconvp(AAhereconv>=1.4); 
    Yyhereconvtt = Yyhereconv(AAhereconv>=1.4);

    Yyyy = (Xxhereconvtt(2:end)-Xxhereconvtt(1:end-1));
    yyyyss = Yyyy(2:end)-Yyyy(1:end-1);
    iretrograde = find(Yyyy<=0);
    ianterograde = find(Yyyy>=0);
    diantero =  ianterograde(2:end)-ianterograde(1:end-1)
    diretro = iretrograde(2:end)-iretrograde(1:end-1)

    if isempty(Xxhereconvtt(ianterograde+1))==1
    else
        if size(find(diantero==1),1)>2
            Xvhere = Xxhereconvtt(ianterograde);
            Yvhere = Yyhereconvtt(ianterograde);
            Timevhere = TThereconvtt(ianterograde);
            velocityHolcman = (Xvhere(2:end)-Xvhere(1:end-1))./(Timevhere(2:end)-Timevhere(1:end-1))./dt;
            VELOCITYantero = [VELOCITYantero;Xxhereconvtt(ianterograde) Yyhereconvtt(ianterograde) [velocityHolcman ; velocityHolcman(end).*ones(2-1,1)]];
        else
        end
    end

    if isempty(Xxhereconvtt(iretrograde+1))==1
    else
        if size(find(diretro==1),1)>2
            Xvhere = Xxhereconvtt(iretrograde);
            Yvhere = Yyhereconvtt(iretrograde);
            Timevhere = TThereconvtt(iretrograde);
            velocityHolcman = (Xvhere(2:end)-Xvhere(1:end-1))./(Timevhere(2:end)-Timevhere(1:end-1))./dt;
            VELOCITYretro= [VELOCITYretro;Xxhereconvtt(iretrograde) Yyhereconvtt(iretrograde) [velocityHolcman ; velocityHolcman(end).*ones(2-1,1)] ];
        else
        end
    end
else
end

end
%%
% Plot the diffusion coefficient in function of the position by averaging
% over 60 consecutive X positions
lower_threshold_transport = 1.4
upper_threshold_diffusion = 1;

% Select Results where alpha par < 1
DHpar = RESULTSEND(RESULTSEND(:,3)<upper_threshold_diffusion,5);
Xpar = RESULTSEND(RESULTSEND(:,3)<upper_threshold_diffusion,1);

% Select results where alphapar < 1.4 and alphaper < 1 
DHper = RESULTSEND(RESULTSEND(:,4)<upper_threshold_diffusion & RESULTSEND(:,3)<lower_threshold_transport,6);
Xper = RESULTSEND(RESULTSEND(:,4)<upper_threshold_diffusion & RESULTSEND(:,3)<lower_threshold_transport,1);

[~,orderici] = sort(Xpar);
Xparo = Xpar(orderici);
DHparo = DHpar(orderici);

ssi = 60; % number of point to average and calculate the std

Diffusion_coefficient_par_mean_std = [];


    for ipar = 1+ssi:ssi: size(Xparo)
Diffusion_coefficient_par_mean_std = [Diffusion_coefficient_par_mean_std;mean(Xparo(ipar-ssi:ipar-1,1)) mean(DHparo(ipar-ssi:ipar-1,1)) std(DHparo(ipar-ssi:ipar-1,1))./sqrt(size(ipar-ssi:ipar-1,2))]
    end
[~,orderici] = sort(Xper);
Xpero = Xper(orderici);
DHpero = DHper(orderici);

Diffusion_coefficient_per_mean_std = [];


    for ipar = 1+ssi:ssi: size(Xpero)
Diffusion_coefficient_per_mean_std = [Diffusion_coefficient_per_mean_std;mean(Xpero(ipar-ssi:ipar-1,1)) mean(DHpero(ipar-ssi:ipar-1,1)) std(DHpero(ipar-ssi:ipar-1,1))./sqrt(size(ipar-ssi:ipar-1,2))]
    end
 
figure(2)
hold off

avvalue = 20; % To smooth the outputs
plot(smooth(Diffusion_coefficient_par_mean_std(1:end,1),1),smooth(Diffusion_coefficient_par_mean_std(1:end,2),avvalue),'-r','linewidth',2)
hold on

plot(smooth(Diffusion_coefficient_par_mean_std(1:end,1),1),smooth(Diffusion_coefficient_par_mean_std(1:end,2)+Diffusion_coefficient_par_mean_std(1:end,3),avvalue),'--r','linewidth',2)
plot(smooth(Diffusion_coefficient_par_mean_std(1:end,1),1),smooth(Diffusion_coefficient_par_mean_std(1:end,2)-Diffusion_coefficient_par_mean_std(1:end,3),avvalue),'--r','linewidth',2)
plot(smooth(Diffusion_coefficient_per_mean_std(1:end,1),1),smooth(Diffusion_coefficient_per_mean_std(1:end,2),avvalue),'-b','linewidth',2)
hold on

plot(smooth(Diffusion_coefficient_per_mean_std(1:end,1),1),smooth(Diffusion_coefficient_per_mean_std(1:end,2)+Diffusion_coefficient_per_mean_std(1:end,3),avvalue),'--b','linewidth',2)
plot(smooth(Diffusion_coefficient_per_mean_std(1:end,1),1),smooth(Diffusion_coefficient_per_mean_std(1:end,2)-Diffusion_coefficient_per_mean_std(1:end,3),avvalue),'--b','linewidth',2)

xlabel('Positions \mum')
ylabel('D \mum^2/s')
set(gca,'fontsize',14)
saveas(2,[FolderSavehere,'Diffusion_coefficient.fig'])
saveas(2,[FolderSavehere,'Diffusion_coefficient.png'])

% Save .mat file
save([FolderSavehere,'Results_compilation.mat'])
%%
Tx8= table(smooth(Diffusion_coefficient_par_mean_std(1:end,1),1),smooth(Diffusion_coefficient_par_mean_std(1:end,2),avvalue),smooth(Diffusion_coefficient_par_mean_std(1:end,2)+Diffusion_coefficient_par_mean_std(1:end,3),avvalue),smooth(Diffusion_coefficient_par_mean_std(1:end,2)-Diffusion_coefficient_par_mean_std(1:end,3),avvalue), 'VariableNames', { 'Position_um','D_par_um2ps', 'D_par_um2ps_plus_sem', 'D_par_um2ps_minus_sem'} )
writetable(Tx8, [FolderSavehere,'Parallel_Diffusion_coefficient_VS_Position.txt'])
% %%
Tx80= table(smooth(Diffusion_coefficient_per_mean_std(1:end,1),1),smooth(Diffusion_coefficient_per_mean_std(1:end,2),avvalue),smooth(Diffusion_coefficient_per_mean_std(1:end,2)+Diffusion_coefficient_per_mean_std(1:end,3),avvalue),smooth(Diffusion_coefficient_per_mean_std(1:end,2)-Diffusion_coefficient_per_mean_std(1:end,3),avvalue), 'VariableNames', { 'Position_um','D_per_um2ps', 'D_per_um2ps_plus_sem', 'D_per_um2ps_minus_sem'} )
writetable(Tx80, [FolderSavehere,'Perpendicular_Diffusion_coefficient_VS_Position.txt'])
%%
% Plot anterograde and retrograde velocities function of the position by averaging
% over 10 consecutive X positions
Xpar = VELOCITYantero(:,1);
Xpar_r = VELOCITYretro(:,1);
Velo_antero = VELOCITYantero(:,3);
Velo_retro = VELOCITYretro(:,3);
[~,orderici] = sort(Xpar);
Xparo = Xpar(orderici);
VHparo = Velo_antero(orderici);

ssi = 10; % Number of points to average and calculate the std

Velocity_anterograde = [];


    for ipar = 1+ssi:ssi: size(Xparo)
Velocity_anterograde = [Velocity_anterograde;mean(Xparo(ipar-ssi:ipar-1,1)) mean(VHparo(ipar-ssi:ipar-1,1)) std(VHparo(ipar-ssi:ipar-1,1))./sqrt(size(ipar-ssi:ipar-1,2))]
    end
    
[~,ordericia] = sort(Xpar_r);
Xparo_r = Xpar_r(ordericia);
VHparo_r = Velo_retro(ordericia);

Velocity_retrograde = [];


    for ipar = 1+ssi:ssi: size(Xparo_r)
Velocity_retrograde = [Velocity_retrograde;mean(Xparo_r(ipar-ssi:ipar-1,1)) mean(VHparo_r(ipar-ssi:ipar-1,1)) std(VHparo_r(ipar-ssi:ipar-1,1))./sqrt(size(ipar-ssi:ipar-1,2))]
    end
 
figure(3)
hold off

avvalue = 20;
plot(smooth(Velocity_anterograde(1:end,1),avvalue),smooth(Velocity_anterograde(1:end,2),avvalue),'-r','linewidth',2)
hold on

plot(smooth(Velocity_anterograde(1:end,1),avvalue),smooth(Velocity_anterograde(1:end,2)+Velocity_anterograde(1:end,3),avvalue),'--r','linewidth',2)
plot(smooth(Velocity_anterograde(1:end,1),avvalue),smooth(Velocity_anterograde(1:end,2)-Velocity_anterograde(1:end,3),avvalue),'--r','linewidth',2)
plot(smooth(Velocity_retrograde(1:end,1),avvalue),smooth(Velocity_retrograde(1:end,2),avvalue),'-b','linewidth',2)
hold on

plot(smooth(Velocity_retrograde(1:end,1),avvalue),smooth(Velocity_retrograde(1:end,2)+Velocity_retrograde(1:end,3),avvalue),'--b','linewidth',2)
plot(smooth(Velocity_retrograde(1:end,1),avvalue),smooth(Velocity_retrograde(1:end,2)-Velocity_retrograde(1:end,3),avvalue),'--b','linewidth',2)

xlabel('Positions \mum')
ylabel('v \mum/s')
set(gca,'fontsize',14)
saveas(3,[FolderSavehere,'Anterograde_and_retrograde_velocities.fig'])
saveas(3,[FolderSavehere,'Anterograde_and_retrograde_velocities.png'])
%%
Tx7= table(smooth(Velocity_anterograde(1:end,1),1),smooth(Velocity_anterograde(1:end,2),avvalue),smooth(Velocity_anterograde(1:end,2)+Velocity_anterograde(1:end,3),avvalue),smooth(Velocity_anterograde(1:end,2)-Velocity_anterograde(1:end,3),avvalue), 'VariableNames', { 'Position_um','V_antero_umps', 'V_antero_umps_plus_sem', 'V_antero_umps_minus_sem'} )
writetable(Tx7, [FolderSavehere,'Anterograde_velocity_VS_Position.txt'])
% %%
Tx70= table(smooth(Velocity_retrograde(1:end,1),1),smooth(Velocity_retrograde(1:end,2),avvalue),smooth(Velocity_retrograde(1:end,2)+Velocity_retrograde(1:end,3),avvalue),smooth(Velocity_retrograde(1:end,2)-Velocity_retrograde(1:end,3),avvalue), 'VariableNames', { 'Position_um','V_retro_umps', 'V_retro_umps_plus_sem', 'V_retro_umps_minus_sem'} )
writetable(Tx70, [FolderSavehere,'Retrograde_velocity_VS_Position.txt'])
%%
figure(1)
[nn,xx] = hist(RESULTSEND(:,3),[-2:0.2:2])
bar(xx,nn./sum(nn),'k')
xlabel('\alpha_{//}')
ylabel('Occurence')
set(gca,'fontsize',20)

% Fit the histogram by two gaussians
XXx2 = xx;
YYy2 = nn./sum(nn);
Paraset1 = [0.5 0.4 0.4 0.2 1.6 0.3  ]; 
curvefitoptions = optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',1E-30,'TolFun',1E-30);
[fittedParametersTotal2,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Fit_two_gaussian,Paraset1, XXx2,YYy2,[],[],curvefitoptions);   
% [Paraset1(1),fittedParametersTotal2(1)];
y_fittedval2 = Fit_two_gaussian(fittedParametersTotal2,[-10: 0.001:5]);
hold on
plot([-10: 0.001:5],fittedParametersTotal2(1).*exp(-([-10: 0.001:5]-fittedParametersTotal2(2)).^2./(2.*(fittedParametersTotal2(3)).^2)),'-c','linewidth',2)
plot([-10: 0.001:5],fittedParametersTotal2(4).*exp(-([-10: 0.001:5]-fittedParametersTotal2(5)).^2./(2.*(fittedParametersTotal2(6)).^2)),'-b','linewidth',2)
plot([-10: 0.001:5],y_fittedval2,'-r','linewidth',2)

set(gca,'fontsize',18)
ylabel('Probability')
xlim([-2 2.2])

% Confidence interval
ci = nlparci(fittedParametersTotal2,residual,'jacobian',jacobian)

saveas(1,[FolderSavehere,'alphapar_all_hist.fig'])
saveas(1,[FolderSavehere,'alphapar_all_hist.png'])


curvecyan = fittedParametersTotal2(1).*exp(-([-10: 0.001:5]-fittedParametersTotal2(2)).^2./(2.*(fittedParametersTotal2(3)).^2));
curveblue =fittedParametersTotal2(4).*exp(-([-10: 0.001:5]-fittedParametersTotal2(5)).^2./(2.*(fittedParametersTotal2(6)).^2))

Percentage_underthecurve_cyan = sum(curvecyan)./(sum(curvecyan)+sum(curveblue))
Percentage_underthecurve_blue = sum(curveblue)./(sum(curvecyan)+sum(curveblue))
Tx5= table([-10: 0.001:5]',curvecyan',curveblue', 'VariableNames', {'Xaxis','curve1','Curve2'} )
writetable(Tx5, [FolderSavehere,'fit_Alpha_par_values.txt'])
Tx11= table(RESULTSEND(:,3), 'VariableNames', {'alpha_par_values'} )
writetable(Tx11, [FolderSavehere,'Alpha_par_values.csv'])
%%
figure(4)
[nn,xx] = hist(RESULTSEND(:,4),[-2:0.2:2])
bar(xx,nn./sum(nn),'k')
xlabel('\alpha_{\perp}')
ylabel('Occurence')
set(gca,'fontsize',20)

% Fit histogram with one gaussian
XXx2 = xx;
YYy2 = nn./sum(nn);
 
Paraset1 = [0.5 0.4 0.4]; 
curvefitoptions = optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',1E-30,'TolFun',1E-30);
[fittedParametersTotal2,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Fit_one_gaussian,Paraset1, XXx2,YYy2,[],[],curvefitoptions);   
[Paraset1(1),fittedParametersTotal2(1)];
y_fittedval2 = Fit_one_gaussian(fittedParametersTotal2,[-10: 0.001:5]);

hold on
plot([-10: 0.001:5],y_fittedval2,'-r','linewidth',2)
threshold = 0.0001
x = [-10: 0.001:5];
set(gca,'fontsize',20)
ylabel('Probability')
xlim([-2 2])
% confidence interval
ci = nlparci(fittedParametersTotal2,residual,'jacobian',jacobian)

saveas(4,[FolderSavehere,'alphaper_all_hist.fig'])
saveas(4,[FolderSavehere,'alphaper_all_hist.png'])

curvecyan = fittedParametersTotal2(1).*exp(-([-10: 0.001:5]-fittedParametersTotal2(2)).^2./(2.*(fittedParametersTotal2(3)).^2));

Tx50= table([-10: 0.001:5]',curvecyan', 'VariableNames', {'Xaxis','curve1'} )
writetable(Tx50, [FolderSavehere,'fit_Alpha_per_values.txt'])
Tx10= table(RESULTSEND(:,4), 'VariableNames', {'alpha_per_values'} )
writetable(Tx10, [FolderSavehere,'Alpha_per_values.csv'])

%%
% figure(10)
% plot(ALPHAGLOBALE15(:,2),ALPHAGLOBALE15(:,3))
       

%%
% Generation of the alpha par and per MAPs - Part 1
Pixelsize   = 0.0535    
siz         = Pixelsize/53.5*25;
wx          = 0.05;
wy          = wx;

out2alphaggpar = RESULTSEND(:,[1,2,3]); 
out2alphaggpar(:,1)   =  out2alphaggpar(:,1)-min(out2alphaggpar(:,1));
out2alphaggpar(:,2)   =  out2alphaggpar(:,2)-min(out2alphaggpar(:,2));

out2alphagg = RESULTSEND(:,[1,2,4]); 
out2alphagg(:,1)   =  out2alphagg(:,1)-min(out2alphagg(:,1));
out2alphagg(:,2)   =  out2alphagg(:,2)-min(out2alphagg(:,2));

outpalphagg = [];
outpalphaggpar = [];

pp   = 150
id_hmmm = []
for it = 1:size(out2alphaggpar,1)-1
    [it ,size(out2alphaggpar,1)]
    x1 = floor(out2alphaggpar(it,1)/siz)+1 +ceil(pp*wx/siz/2);
    x2 = floor(out2alphaggpar(it+1,1)/siz)+1+ceil(pp*wx/siz/2);
    y1 = floor(out2alphaggpar(it,2)/siz)  +1+ceil(pp/2*wx/siz);
    y2 = floor(out2alphaggpar(it+1,2)/siz)+1+ceil(pp/2*wx/siz);
    if [(sqrt((out2alphaggpar(it,1)-out2alphaggpar(it+1,1)).^2 + (out2alphaggpar(it,2)-out2alphaggpar(it+1,2)).^2)*1) < 0.5]
        x = round(linspace(x1,x2,20)');
        y = round(linspace(y1,y2,20)');
        x2 = x;
        y2 = y;
        it2 = 1;
        while it2 < length(x)
            posEq = find(sum(ismember([x,y],[x(it2),y(it2)]),2)==2);
            if length(posEq) > 1
                iii = [1:length(x)];
                iii = setdiff(iii,setdiff(posEq,it2));
                x = x(iii);
                y = y(iii);
            else
                it2 = it2 +1;
            end
        end

         outpalphagg = [outpalphagg; [x,y,out2alphagg(it,3)+0*x]];

         outpalphaggpar = [outpalphaggpar; [x,y,out2alphaggpar(it,3)+0*x]];

    else
        id_hmmm = [id_hmmm,it];
        '        hmmm ... '

         outpalphagg = [outpalphagg; [x1,y1,out2alphagg(it,3)+0*x1]]; % avec alpha global

         outpalphaggpar = [outpalphaggpar; [x1,y1,out2alphaggpar(it,3)+0*x1]]; % avec alpha global

    end
end
%%
% Superresolved image of all the localizations covoluted by a gaussian 
[XX,YY]     = meshgrid([1:(ceil(max(out2alphaggpar(:,1)/siz))+ceil(pp*wx/siz))],[1:(ceil(max(out2alphaggpar(:,2)/siz))+ceil(pp*wx/siz))]);
imaSuperresolution_p       = full(sparse(floor((out2alphaggpar(:,1))/siz)+1+ceil(pp*wx/siz/2),floor((out2alphaggpar(:,2))/siz)+1+ceil(pp*wx/siz/2), 1+0*out2alphaggpar(:,2),size(XX,2),size(XX,1)));

beta        = [1,ceil(3*wx/siz),ceil(3*wx/siz),wx/siz,wy/siz,0];
[XXg,YYg]   = meshgrid(1:ceil(6*wx/siz),1:ceil(6*wx/siz));
[~,Gauss2D] = gaussian2d(beta,[YYg(:);XXg(:)],XXg);
Gauss2D     = Gauss2D/sum(Gauss2D(:));
imaSuperresolution   = conv2(imaSuperresolution_p,Gauss2D,'same');
figure(100)
imshow(imaSuperresolution',[]); 
saveas(100,[FolderSavehere,'_Superresolvedmapc.tif'])
saveas(100,[FolderSavehere,'_Superresolvedmapc.png'],'png')
%%
% Generation of the maps for alpha par and per values - Part 2
imaAlphaggpar        = full(sparse(outpalphaggpar(outpalphaggpar(:,3)>0,1),outpalphaggpar(outpalphaggpar(:,3)>0,2),outpalphaggpar(outpalphaggpar(:,3)>0,3) ,size(XX,2),size(XX,1)));
imaAlphanggpar       = full(sparse(outpalphaggpar(outpalphaggpar(:,3)>0,1),outpalphaggpar(outpalphaggpar(:,3)>0,2), 1+0*outpalphaggpar(outpalphaggpar(:,3)>0,3),size(XX,2),size(XX,1)));
imaAlphaggpar = (imaAlphaggpar);
minN = 1;
imaAlphaggpar(imaAlphanggpar>=minN)  = imaAlphaggpar(imaAlphanggpar>=minN)./imaAlphanggpar(imaAlphanggpar>=minN);


imaAlphagg        = full(sparse(outpalphagg(outpalphaggpar(:,3)>0,1),outpalphagg(outpalphaggpar(:,3)>0,2),outpalphagg(outpalphaggpar(:,3)>0,3) ,size(XX,2),size(XX,1)));
imaAlphangg       = full(sparse(outpalphagg(outpalphaggpar(:,3)>0,1),outpalphagg(outpalphaggpar(:,3)>0,2), 1+0*outpalphagg(outpalphaggpar(:,3)>0,3),size(XX,2),size(XX,1)));
imaAlphagg = (imaAlphagg);
minN = 1;
imaAlphagg(imaAlphangg>=minN)  = imaAlphagg(imaAlphangg>=minN)./imaAlphangg(imaAlphangg>=minN);
%%
% Alpha par MAP
Scale = 1./siz; % 1 um
imaAlphaggcpar = imaAlphaggpar;

imaAlphaggcpar(50:50+Scale,298:298+10)= 256
 fi = figure(40)
set(gca,'color','none')
imshow([imaAlphaggcpar'],[0 2])
Numbits = 12

maph = colormap(parula(2^Numbits));
maph(1,:)=0;
axis equal
colormap(maph);
cn= colorbar;
ylabel(cn,['\alpha_{//} '],'fontsize',24,'rotation',90)
freezeColors
whitebg('white')
set(gca,'fontsize',16)

saveas(40,[FolderSavehere,'_Alphamappar.fig'])
saveas(40,[FolderSavehere,'_Alphamappar.png'],'png')

%%
% Alpha per MAP
Scale = 1./siz;
imaAlphaggc = imaAlphagg

imaAlphaggc(50:50+Scale,298:298+10)= 256
 fi = figure(41)
set(gca,'color','none')
imshow([imaAlphaggc'],[0 2])
Numbits = 12
maph = colormap(parula(2^Numbits));
maph(1,:)=0;
axis equal

colormap(maph);
cn= colorbar;
ylabel(cn,['\alpha_{\perp} '],'fontsize',24,'rotation',90)
freezeColors
whitebg('white')
set(gca,'fontsize',16)

saveas(41,[FolderSavehere,'_Alphamapper.fig'])
saveas(41,[FolderSavehere,'_Alphamapper.png'],'png')
%%
%%
