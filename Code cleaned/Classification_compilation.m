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
% Position when this is not transport
% Xxherecondtt = Xxhereconv(AAhereconv<threshold); % Selection de position when alphaparallel up to 1.3
% Yyherecondtt = Yyhereconv(AAhereconv<threshold); % Selection de position when alphaparallel up to 1.3

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
saveas(96,[FolderSavehere,'Diffusion_coefficient.fig'])
saveas(96,[FolderSavehere,'Diffusion_coefficient.png'])

% Save .mat file
save([FolderSavehere,'Results.mat'])
%%
%  Tx8= table(smooth(Distancemeande_e5a1(1:end-1,1),1),smooth(Distancemeande_e5a1(1:end-1,2),avvalue),smooth(Distancemeande_e5a1(1:end-1,2)+Distancemeande_e5a1(1:end-1,3),avvalue),smooth(Distancemeande_e5a1(1:end-1,2)-Distancemeande_e5a1(1:end-1,3),avvalue), 'VariableNames', { 'Position_um','D_par_um2ps', 'D_par_um2ps_plus_sem', 'D_par_um2ps_minus_sem'} )
% writetable(Tx8, [FolderSavehere,'Parallel_Diffusion_coefficient_averagedon60pts_Nperp_5228_Npar_4418_VS_Position.txt'])
% %%
%  Tx8= table(smooth(Distancemeande_e5a2(1:end-1,1),1),smooth(Distancemeande_e5a2(1:end-1,2),avvalue),smooth(Distancemeande_e5a2(1:end-1,2)+Distancemeande_e5a2(1:end-1,3),avvalue),smooth(Distancemeande_e5a2(1:end-1,2)-Distancemeande_e5a2(1:end-1,3),avvalue), 'VariableNames', { 'Position_um','D_par_um2ps', 'D_par_um2ps_plus_sem', 'D_par_um2ps_minus_sem'} )
% writetable(Tx8, [FolderSavehere,'Perpendicular_Diffusion_coefficient_averagedon60pts_Nperp_5228_Npar_4418_VS_Position.txt'])
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
%%

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
% Tx5= table([-10: 0.001:5]',curvecyan',curveblue', 'VariableNames', {'Xaxis','curve1','Curve2'} )
% writetable(Tx5, [FolderSavehere,'fit_Alpha_par_all.txt'])
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

Tx5= table([-10: 0.001:5]',curvecyan', 'VariableNames', {'Xaxis','curve1'} )
writetable(Tx5, [FolderSavehere,'fit_Alpha_per_all.txt'])
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
 fi = figure(4)
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

saveas(4,[FolderSavehere,'_Alphamappar.fig'])
saveas(4,[FolderSavehere,'_Alphamappar.png'],'png')

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
% ALPHAGLOBALE(ALPHAGLOBALE(:,4)<=0.55 & ALPHAGLOBALE(:,4)>=0.45,[5:13+5])
lPS = 3.5;
lTZ = 1.5;
lDS = 6.5;
MSDPERGLOBPS =  ALPHAGLOBALE(ALPHAGLOBALE(:,4)<=0.55 & ALPHAGLOBALE(:,4)>=0.45 & ALPHAGLOBALE(:,1)>=lTZ & ALPHAGLOBALE(:,1)<=lPS,[15+13+5:15+13+5+13]);
MSDPERGLOBPS = MSDPERGLOBPS(MSDPERGLOBPS(:,end)<100,:) ;
MSDPARGLOBPS =  ALPHAGLOBALE(ALPHAGLOBALE(:,3)<=0.55 & ALPHAGLOBALE(:,3)>=0.45 & ALPHAGLOBALE(:,1)>=lTZ & ALPHAGLOBALE(:,1)<=lPS,[5:15+5-2]);
MSDPARGLOBPS = MSDPARGLOBPS(MSDPARGLOBPS(:,end)<100,:) ;

MSDPERGLOBDS =  ALPHAGLOBALE(ALPHAGLOBALE(:,4)<=0.55 & ALPHAGLOBALE(:,4)>=0.45 & ALPHAGLOBALE(:,1)>=lPS & ALPHAGLOBALE(:,1)<=lDS,[15+13+5:15+13+5+13]);
MSDPERGLOBDS = MSDPERGLOBDS(MSDPERGLOBDS(:,end)<100,:) ;
MSDPARGLOBDS =  ALPHAGLOBALE(ALPHAGLOBALE(:,3)<=0.55 & ALPHAGLOBALE(:,3)>=0.45 & ALPHAGLOBALE(:,1)>=lPS & ALPHAGLOBALE(:,1)<=lDS,[5:15+5-2]);
MSDPARGLOBDS = MSDPARGLOBDS(MSDPARGLOBDS(:,end)<100,:) ;

MSDPERGLOBt =  ALPHAGLOBALE(ALPHAGLOBALE(:,4)<=0.55 & ALPHAGLOBALE(:,4)>=0.45 & ALPHAGLOBALE(:,1)>=lDS,[15+13+5:15+13+5+13]);
MSDPERGLOBt = MSDPERGLOBt(MSDPERGLOBt(:,end)<100,:) ;
MSDPARGLOBt =  ALPHAGLOBALE(ALPHAGLOBALE(:,3)<=0.55 & ALPHAGLOBALE(:,3)>=0.45 & ALPHAGLOBALE(:,1)>=lDS,[5:15+5-2]);
MSDPARGLOBt = MSDPARGLOBt(MSDPARGLOBt(:,end)<100,:) ; 
%  mean(MSDPERGLOB)
MSDPERGLOBtz =  ALPHAGLOBALE(ALPHAGLOBALE(:,4)<=0.55 & ALPHAGLOBALE(:,4)>=0.45 & ALPHAGLOBALE(:,1)<=lTZ,[15+13+5:15+13+5+13]);
MSDPERGLOBtz = MSDPERGLOBtz(MSDPERGLOBtz(:,end)<100,:) ;
MSDPARGLOBtz =  ALPHAGLOBALE(ALPHAGLOBALE(:,3)<=0.55 & ALPHAGLOBALE(:,3)>=0.45 & ALPHAGLOBALE(:,1)<=lTZ,[5:15+5-2]);
MSDPARGLOBtz = MSDPARGLOBtz(MSDPARGLOBtz(:,end)<100,:) ; 

 %%
 %default =0.5
 alphalinmhere = +0.5 
 alphalinmheretipper = 0.36-0.5;
 alphalinmheretippar = 0.46-0.5;
  alphalinmheretzper = 0.39-0.5;
 alphalinmheretzpar = 0.65-0.5;
  alphalinmhereDSper = 0.39-0.5;
 alphalinmhereDSpar = 0.41-0.5;
   alphalinmherePSper = 0.41-0.5;
 alphalinmherePSpar = 0.47-0.5;
 alphalinm = alphalinmhere+0.5;
MSDPERGLOBPS =  ALPHAGLOBALE(ALPHAGLOBALE(:,3)<1.4 & ALPHAGLOBALE(:,4)<=0.55+alphalinmherePSper & ALPHAGLOBALE(:,4)>=0.45+alphalinmherePSper & ALPHAGLOBALE(:,1)>=lTZ & ALPHAGLOBALE(:,1)<=lPS,[15+13+5:15+13+5+13]);
MSDPERGLOBPS = MSDPERGLOBPS(MSDPERGLOBPS(:,end)<100,:) ;
MSDPARGLOBPS =  ALPHAGLOBALE(ALPHAGLOBALE(:,3)<=0.55+alphalinmherePSpar & ALPHAGLOBALE(:,3)>=0.45+alphalinmherePSpar & ALPHAGLOBALE(:,1)>=lTZ & ALPHAGLOBALE(:,1)<=lPS,[5:15+5-2]);
MSDPARGLOBPS = MSDPARGLOBPS(MSDPARGLOBPS(:,end)<100,:) ;

MSDPERGLOBDS =  ALPHAGLOBALE(ALPHAGLOBALE(:,3)<1.4 & ALPHAGLOBALE(:,4)<=0.55+alphalinmhereDSper & ALPHAGLOBALE(:,4)>=0.45+alphalinmhereDSper & ALPHAGLOBALE(:,1)>=lPS & ALPHAGLOBALE(:,1)<=lDS,[15+13+5:15+13+5+13]);
MSDPERGLOBDS = MSDPERGLOBDS(MSDPERGLOBDS(:,end)<100,:) ;
MSDPARGLOBDS =  ALPHAGLOBALE(ALPHAGLOBALE(:,3)<=0.55+alphalinmhereDSpar & ALPHAGLOBALE(:,3)>=0.45+alphalinmhereDSpar & ALPHAGLOBALE(:,1)>=lPS & ALPHAGLOBALE(:,1)<=lDS,[5:15+5-2]);
MSDPARGLOBDS = MSDPARGLOBDS(MSDPARGLOBDS(:,end)<100,:) ;

MSDPERGLOBt =  ALPHAGLOBALE(ALPHAGLOBALE(:,3)<1.4 & ALPHAGLOBALE(:,4)<=0.55+alphalinmheretipper & ALPHAGLOBALE(:,4)>=0.45+alphalinmheretipper & ALPHAGLOBALE(:,1)>=lDS,[15+13+5:15+13+5+13]);
MSDPERGLOBt = MSDPERGLOBt(MSDPERGLOBt(:,end)<100,:) ;
MSDPARGLOBt =  ALPHAGLOBALE(ALPHAGLOBALE(:,3)<=0.55+alphalinmheretippar & ALPHAGLOBALE(:,3)>=0.45+alphalinmheretippar & ALPHAGLOBALE(:,1)>=lDS,[5:15+5-2]);
MSDPARGLOBt = MSDPARGLOBt(MSDPARGLOBt(:,end)<100,:) ; 
%  mean(MSDPERGLOB)
MSDPERGLOBtz =  ALPHAGLOBALE(ALPHAGLOBALE(:,3)<1.4 & ALPHAGLOBALE(:,4)<=0.55+alphalinmheretzper & ALPHAGLOBALE(:,4)>=0.45+alphalinmheretzper & ALPHAGLOBALE(:,1)<=lTZ,[15+13+5:15+13+5+13]);
MSDPERGLOBtz = MSDPERGLOBtz(MSDPERGLOBtz(:,end)<100,:) ;
MSDPARGLOBtz =  ALPHAGLOBALE(ALPHAGLOBALE(:,3)<=0.55+alphalinmheretzpar & ALPHAGLOBALE(:,3)>=0.45+alphalinmheretzpar & ALPHAGLOBALE(:,1)<=lTZ,[5:15+5-2]);
MSDPARGLOBtz = MSDPARGLOBtz(MSDPARGLOBtz(:,end)<100,:) ; 
Time  = [1:size(mean(MSDPARGLOBtz)',1)].*dt
%  Tx2= table(Time',mean(MSDPARGLOBtz)',mean(MSDPERGLOBtz)',mean(MSDPARGLOBPS)',mean(MSDPERGLOBPS)',mean(MSDPARGLOBDS)',mean(MSDPERGLOBDS)',mean(MSDPARGLOBt)',mean(MSDPERGLOBt)', 'VariableNames', { 'Lag_Time_s','MSD_par_TZ', 'MSD_per_TZ', 'MSD_par_PS', 'MSD_per_PS', 'MSD_par_DS', 'MSD_per_DS', 'MSD_par_tip', 'MSD_per_tip'} )
% writetable(Tx2, [FolderSavehere,'MSD_confined_par_per.txt'])

 %%
 figure(1)
 hold off
 plot([1:1:14].*dt,mean(MSDPERGLOBPS),'ob','linewidth',2)
 xlabel('\tau')
 ylabel('MSD')
 hold on
 plot([1:1:14].*dt,mean(MSDPARGLOBPS),'or','linewidth',2)

 title(['PS \alpha_{//} = ' num2str(alphalinmherePSpar+0.5,'%.2f'),' \alpha_{\perp}  = ' num2str(alphalinmherePSper+0.5,'%.2f')])
 set(gca,'fontsize',20)

 XXTotal = [1:1:6].*dt; % Je prends 30% du temps total
ydataMSDTotal = mean(MSDPERGLOBPS(:,[1:6]));
Paraset1 = [mean(MSDPERGLOBPS(:,1)),mean(MSDPERGLOBPS(:,6))];

curvefitoptions = optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',1E-30,'TolFun',1E-30);
[fittedParametersTotal,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Fit_xi_MSD,Paraset1, XXTotal,ydataMSDTotal,[],[],curvefitoptions);   
[Paraset1(1),fittedParametersTotal(1)];
fit_MSD_per_PS = Fit_xi_MSD(fittedParametersTotal,[0:1:10].*dt)
curvefitoptions = optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',1E-30,'TolFun',1E-30);
[fittedParametersTotal2,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Fit_xi_MSD_2nd_order,Paraset1, XXTotal,ydataMSDTotal,[],[],curvefitoptions);   
[Paraset1(1),fittedParametersTotal(1)];
fit_MSD_per_PS2 = Fit_xi_MSD_2nd_order(fittedParametersTotal2,[0:1:10].*dt)

 hold on
 cipar = nlparci(fittedParametersTotal,residual,'jacobian',jacobian)
 dl = fittedParametersTotal(1)'-cipar(1)
 dxi = 1./(2.*sqrt(fittedParametersTotal(1))).*dl*1000
%  Time  = [0:1:10].*dt
%  Tx2= table(Time',fit_MSD_per_PS', 'VariableNames', { 'Lag_Time_s','fit_MSD_per_PS', 'MSD_per_TZ', 'MSD_par_PS', 'MSD_per_PS', 'MSD_par_DS', 'MSD_per_DS', 'MSD_par_tip', 'MSD_per_tip'} )
% writetable(Tx2, [FolderSavehere,'fit_MSDper_confined_PS.txt'])


 plot([0:1:10].*dt,fit_MSD_per_PS,'-b','linewidth',2)
  XXTotal = [1:1:6].*dt; % Je prends 30% du temps total
ydataMSDTotal = mean(MSDPARGLOBPS(:,[1:6]));
Paraset1 = [mean(MSDPARGLOBPS(:,1)),mean(MSDPARGLOBPS(:,6))];

curvefitoptions = optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',1E-30,'TolFun',1E-30);
[fittedParametersTotalp,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Fit_xi_MSD,Paraset1, XXTotal,ydataMSDTotal,[],[],curvefitoptions);   
[Paraset1(1),fittedParametersTotalp(1)];
fit_MSD_par_PS = Fit_xi_MSD(fittedParametersTotalp,[0:1:10].*dt)
curvefitoptions = optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',1E-30,'TolFun',1E-30);
[fittedParametersTotalp2,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Fit_xi_MSD_2nd_order,Paraset1, XXTotal,ydataMSDTotal,[],[],curvefitoptions);   
[Paraset1(1),fittedParametersTotalp(1)];
fit_MSD_par_PS2 = Fit_xi_MSD_2nd_order(fittedParametersTotalp2,[0:1:10].*dt)
 hold on
 plot([0:1:10].*dt,fit_MSD_par_PS,'-r','linewidth',2)
  legend(['\xi_{\perp} = ' num2str(sqrt(fittedParametersTotal(1))*1000,'%.1f'),' nm'],['\xi_{//} = ' num2str(sqrt(fittedParametersTotalp(1))*1000,'%.1f'),' nm'])
cipar = nlparci(fittedParametersTotalp,residual,'jacobian',jacobian)
 dl = fittedParametersTotalp(1)'-cipar(1)
 dxi = 1./(2.*sqrt(fittedParametersTotalp(1))).*dl*1000
%   saveas(1,[FolderSavehere,'alphaPS =', num2str(alphalinm,'%.1f'),'xi_withalphaparper2.fig'])
% saveas(1,[FolderSavehere,'alphaPS =', num2str(alphalinm,'%.1f'),'xi_withalphaparper2.png'])
[sqrt(fittedParametersTotal(1))*1000 sqrt(fittedParametersTotal2(1))*1000;sqrt(fittedParametersTotalp(1))*1000 sqrt(fittedParametersTotalp2(1))*1000]
[fittedParametersTotal(2) fittedParametersTotal(2); fittedParametersTotalp(2) fittedParametersTotalp2(2)]
   %%
 figure(2)
 hold off
 plot([1:1:14].*dt,mean(MSDPERGLOBDS),'ob','linewidth',2)
 xlabel('\tau')
 ylabel('MSD')
 hold on
 plot([1:1:14].*dt,mean(MSDPARGLOBDS),'or','linewidth',2)

 title(['DS \alpha_{//} = ' num2str(alphalinmhereDSpar+0.5,'%.2f'),' \alpha_{\perp}  = ' num2str(alphalinmhereDSper+0.5,'%.2f')])
 set(gca,'fontsize',20)
%  xlim([0 0.6])
 XXTotal = [1:1:6].*dt; % Je prends 30% du temps total
ydataMSDTotal = mean(MSDPERGLOBDS(:,[1:6]));
Paraset1 = [mean(MSDPERGLOBDS(:,1)),mean(MSDPERGLOBDS(:,6))];

curvefitoptions = optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',1E-30,'TolFun',1E-30);
[fittedParametersTotalD,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Fit_xi_MSD,Paraset1, XXTotal,ydataMSDTotal,[],[],curvefitoptions);   
% [Paraset1(1),fittedParametersTotal(1)];
fit_MSD_per_DS = Fit_xi_MSD(fittedParametersTotalD,[0:1:10].*dt)
 hold on
 cipar = nlparci(fittedParametersTotalD,residual,'jacobian',jacobian)
 dl = fittedParametersTotalD(1)'-cipar(1)
 dxi = 1./(2.*sqrt(fittedParametersTotalD(1))).*dl*1000
 plot([0:1:10].*dt,fit_MSD_per_DS,'-b','linewidth',2)
  XXTotal = [1:1:6].*dt; % Je prends 30% du temps total
ydataMSDTotal = mean(MSDPARGLOBDS(:,[1:6]));
Paraset1 = [mean(MSDPARGLOBDS(:,1)),mean(MSDPARGLOBDS(:,6))];

curvefitoptions = optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',1E-30,'TolFun',1E-30);
[fittedParametersTotalDp,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Fit_xi_MSD,Paraset1, XXTotal,ydataMSDTotal,[],[],curvefitoptions);   
% [Paraset1(1),fittedParametersTotalp(1)];
fit_MSD_par_DS = Fit_xi_MSD(fittedParametersTotalDp,[0:1:10].*dt)
 hold on
 plot([0:1:10].*dt,fit_MSD_par_DS,'-r','linewidth',2)
  legend(['\xi_{\perp} = ' num2str(sqrt(fittedParametersTotalD(1))*1000,'%.1f'),' nm'],['\xi_{//} = ' num2str(sqrt(fittedParametersTotalDp(1))*1000,'%.1f'),' nm'])
%   saveas(2,[FolderSavehere,'alphaDS =', num2str(alphalinm,'%.1f'),'xi_withalphaparper2.fig'])
% saveas(2,[FolderSavehere,'alphaDS =', num2str(alphalinm,'%.1f'),'xi_withalphaparper2.png'])
cipar = nlparci(fittedParametersTotalDp,residual,'jacobian',jacobian)
 dl = fittedParametersTotalDp(1)'-cipar(1)
 dxi = 1./(2.*sqrt(fittedParametersTotalDp(1))).*dl*1000
[fittedParametersTotalD(2); fittedParametersTotalDp(2)]

  %%
 figure(3)
 hold off
 plot([1:1:14].*dt,mean(MSDPERGLOBt),'ob','linewidth',2)
 xlabel('\tau')
 ylabel('MSD')
 hold on
 plot([1:1:14].*dt,mean(MSDPARGLOBt),'or','linewidth',2)

 title(['tip \alpha_{//} = ' num2str(alphalinmheretippar+0.5,'%.2f'),' \alpha_{\perp}  = ' num2str(alphalinmheretipper+0.5,'%.2f')])
 set(gca,'fontsize',20)
%  xlim([0 0.6])
 XXTotal = [1:1:6].*dt; % Je prends 30% du temps total
ydataMSDTotal = mean(MSDPERGLOBt(:,[1:6]));
Paraset1 = [mean(MSDPERGLOBt(:,1)),mean(MSDPERGLOBt(:,6))];

curvefitoptions = optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',1E-30,'TolFun',1E-30);
[fittedParametersTotalt,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Fit_xi_MSD,Paraset1, XXTotal,ydataMSDTotal,[],[],curvefitoptions);   
% [Paraset1(1),fittedParametersTotal(1)];
fit_MSD_per_tip = Fit_xi_MSD(fittedParametersTotalt,[0:1:10].*dt)
 hold on
 cipar = nlparci(fittedParametersTotalt,residual,'jacobian',jacobian)
 dl = fittedParametersTotalt(1)'-cipar(1)
 dxi = 1./(2.*sqrt(fittedParametersTotalt(1))).*dl*1000
[fittedParametersTotalt(2); fittedParametersTotaltp(2)]

 %%
 plot([0:1:10].*dt,fit_MSD_per_tip,'-b','linewidth',2)
  XXTotal = [1:1:6].*dt; % Je prends 30% du temps total
ydataMSDTotal = mean(MSDPARGLOBt(:,[1:6]));
Paraset1 = [mean(MSDPARGLOBt(:,1)),mean(MSDPARGLOBt(:,6))];

curvefitoptions = optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',1E-30,'TolFun',1E-30);
[fittedParametersTotaltp,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Fit_xi_MSD,Paraset1, XXTotal,ydataMSDTotal,[],[],curvefitoptions);   
% [Paraset1(1),fittedParametersTotalp(1)];
fit_MSD_par_tip = Fit_xi_MSD(fittedParametersTotaltp,[0:1:10].*dt)
 hold on
 plot([0:1:10].*dt,fit_MSD_par_tip,'-r','linewidth',2)
  legend(['\xi_{\perp} = ' num2str(sqrt(fittedParametersTotalt(1))*1000,'%.1f'),' nm'],['\xi_{//} = ' num2str(sqrt(fittedParametersTotaltp(1))*1000,'%.1f'),' nm'])
%   saveas(3,[FolderSavehere,'alphatip =', num2str(alphalinm,'%.1f'),'xi_withalphaparper2.fig'])
% saveas(3,[FolderSavehere,'alphatip =', num2str(alphalinm,'%.1f'),'xi_withalphaparper2.png'])
cipar = nlparci(fittedParametersTotaltp,residual,'jacobian',jacobian)
 dl = fittedParametersTotaltp(1)'-cipar(1)
 dxi = 1./(2.*sqrt(fittedParametersTotaltp(1))).*dl*1000
  %%
 figure(4)
 hold off
 plot([1:1:14].*dt,mean(MSDPERGLOBtz),'ob','linewidth',2)
 xlabel('\tau')
 ylabel('MSD')
 hold on
 plot([1:1:14].*dt,mean(MSDPARGLOBtz),'or','linewidth',2)
 title(['TZ \alpha_{//} = ' num2str(alphalinmheretzpar+0.5,'%.2f'),' \alpha_{\perp}  = ' num2str(alphalinmheretzper+0.5,'%.2f')])
 set(gca,'fontsize',20)
%  xlim([0 0.6])
 XXTotal = [1:1:6].*dt; % Je prends 30% du temps total
ydataMSDTotal = mean(MSDPERGLOBtz(:,[1:6]));
Paraset1 = [mean(MSDPERGLOBtz(:,1)),mean(MSDPERGLOBtz(:,6))];

curvefitoptions = optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',1E-30,'TolFun',1E-30);
[fittedParametersTotaltz,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Fit_xi_MSD,Paraset1, XXTotal,ydataMSDTotal,[],[],curvefitoptions);   
% [Paraset1(1),fittedParametersTotal(1)];
fit_MSD_per_TZ = Fit_xi_MSD(fittedParametersTotaltz,[0:1:10].*dt)
 hold on
 ci = nlparci(fittedParametersTotaltz,residual,'jacobian',jacobian)
 dl = fittedParametersTotaltz(1)'-ci(1)
 dxi = 1./(2.*sqrt(fittedParametersTotaltz(1))).*dl*1000
 plot([0:1:10].*dt,fit_MSD_per_TZ,'-b','linewidth',2)
  XXTotal = [1:1:6].*dt; % Je prends 30% du temps total
ydataMSDTotal = mean(MSDPARGLOBtz(:,[1:6]));
Paraset1 = [mean(MSDPARGLOBtz(:,1)),mean(MSDPARGLOBtz(:,6))];

curvefitoptions = optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000,'TolX',1E-30,'TolFun',1E-30);
[fittedParametersTotaltpz,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@Fit_xi_MSD,Paraset1, XXTotal,ydataMSDTotal,[],[],curvefitoptions);   
% [Paraset1(1),fittedParametersTotalp(1)];
fit_MSD_par_TZ = Fit_xi_MSD(fittedParametersTotaltpz,[0:1:10].*dt)
 hold on
 plot([0:1:10].*dt,fit_MSD_par_TZ,'-r','linewidth',2)
  legend(['\xi_{\perp} = ' num2str(sqrt(fittedParametersTotaltz(1))*1000,'%.1f'),' nm'],['\xi_{//} = ' num2str(sqrt(fittedParametersTotaltpz(1))*1000,'%.1f'),' nm'])
%   saveas(4,[FolderSavehere,'alphaTZ =', num2str(alphalinm,'%.1f'),'xi_withalphaparper2.fig'])
% saveas(4,[FolderSavehere,'alphaTZ =', num2str(alphalinm,'%.1f'),'xi_withalphaparper2.png'])
[fittedParametersTotaltz(2); fittedParametersTotaltpz(2)]

   Time  = [0:1:10].*dt
 Tx2= table(Time',fit_MSD_par_TZ',fit_MSD_per_TZ',fit_MSD_par_PS',fit_MSD_per_PS',fit_MSD_par_DS',fit_MSD_per_DS',fit_MSD_par_tip',fit_MSD_per_tip', 'VariableNames', { 'Lag_Time_s','fit_MSD_par_TZ', 'fit_MSD_per_TZ', 'fit_MSD_par_PS', 'fit_MSD_per_PS', 'fit_MSD_par_DS', 'fit_MSD_per_DS', 'fit_MSD_par_tip', 'fit_MSD_per_tip'} )
writetable(Tx2, [FolderSavehere,'fit_MSDperpar_confined.txt'])

fittedParametersTotaltpz-cipar
cipar = nlparci(fittedParametersTotaltpz,residual,'jacobian',jacobian)
 dl = fittedParametersTotaltpz(1)'-cipar(1)
 dxi = 1./(2.*sqrt(fittedParametersTotaltpz(1))).*dl*1000

