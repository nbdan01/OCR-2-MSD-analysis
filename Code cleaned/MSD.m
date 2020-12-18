function [LagT,MSD_mean,MSD_std,yfit,p] = MSD(pos,T,numberoflags,time_image,LagT)

if (max(size(pos)))^2 > 1%1*2500^2   
    indice = round(T/time_image);
    X = zeros(max(indice),1)+NaN; % J'ai changé length par max
    X(indice) = pos(:,1);    
    T = time_image*[1:length(X)];
    
    if size(X,1) ~= size(T,1)
        T = T';        
    end
    if size(pos,2) >= 2
        Y = zeros(max(indice),1)+NaN;
        Y(indice) = pos(:,2);
    end
    if size(pos,2) >= 3
        Z = zeros(max(indice),1)+NaN;
        Z(indice) = pos(:,3);
    end
    
    if max(size(numberoflags)) > 1
        LagT     = numberoflags;
        MSD_mean = 0*LagT;
        MSD_std  = 0*LagT;
    
        for lagit = 1:length(numberoflags)       
            lag = find(abs(T-LagT(lagit)) < eps);
%             LagT(lag) = lag*time_image;
            dX = X(lag+1:end) - X(1:end-lag)  ;  
            if size(pos,2) >= 2
                dY = Y(lag+1:end) - Y(1:end-lag); 
            else
                dY = 0*dX; 
            end
            if size(pos,2) >= 3
                dZ = Z(lag+1:end) - Z(1:end-lag); 
            else
                dZ = 0*dX; 
            end        
            d2   = (dX.^2+dY.^2+dZ.^2);   
            iiii = find(~isnan(d2));     
            MSD_mean(lagit) = [mean(d2(~isnan(d2)))]; 
            MSD_std(lagit)  = [std(d2(~isnan(d2)))/sqrt(length(d2(~isnan(d2))))];         
        end    
    else
        MSD_mean = zeros(1,numberoflags);
        MSD_std  = zeros(1,numberoflags);
        LagT     = zeros(1,numberoflags);
    
        for lagit = 1:numberoflags 
            LagT(lagit) = lagit*time_image;
            
            lag = find(abs(T-LagT(lagit)) < eps);
            
            dX = X(lag+1:end) - X(1:end-lag)  ;  
            if size(pos,2) >= 2
                dY = Y(lag+1:end) - Y(1:end-lag); 
            else
                dY = 0*dX; 
            end
            if size(pos,2) >= 3
                dZ = Z(lag+1:end) - Z(1:end-lag); 
            else
                dZ = 0*dX; 
            end        
            d2   = (dX.^2+dY.^2+dZ.^2);      
            iiii = find(~isnan(d2));
            MSD_mean(lag) = [mean(d2(~isnan(d2)))]; 
            MSD_std(lag)  = [std(d2(~isnan(d2)))/sqrt(length(d2(~isnan(d2))))];         
        end    
    end
    ii=find(abs(MSD_mean)>eps);
    LagT     = LagT(ii);
    MSD_mean = MSD_mean(ii);
    MSD_std  = MSD_std(ii);
    
    p = polyfit(LagT,MSD_mean,1);
    yfit = polyval(p,LagT);
    
%     MSD_mean2 = MSD_mean;
%     MSD_std2  = MSD_std;
%     LagT2     = LagT;
else
    X = pos(:,1);
    if size(X,1) ~= size(T,1)
        T = T';
    end
    dt  = time_image*(repmat(T',[length(T),1])-repmat(T,[1,length(T)]));
    dX  = repmat(X',[length(X),1])-repmat(X,[1,length(X)]);
    if size(pos,2) == 2
        Y = pos(:,2);
        dY  = repmat(Y',[length(Y),1])-repmat(Y,[1,length(Y)]);
        d2   = (dX.^2+dY.^2);
    elseif size(pos,2) == 3        
        Y = pos(:,2);  
        Z = pos(:,3);
        dY  = repmat(Y',[length(Y),1])-repmat(Y,[1,length(Y)]);
        dZ  = repmat(Z',[length(Z),1])-repmat(Z,[1,length(Z)]);
        d2   = (dX.^2+dY.^2+dZ.^2);
    else        
        d2   = dX.^2;
    end
    
    MSD_mean = zeros(1,numberoflags);
    MSD_std  = zeros(1,numberoflags);
    LagT     = zeros(1,numberoflags);
%     
%     eps = 1E-3;
    for lag = 1:numberoflags
        ii       = find(abs(dt-lag*time_image)<eps);
        LagT(lag)     = [lag*time_image];
        MSD_mean(lag) = [ mean(d2(ii))]; 
        MSD_std(lag)  = [std(d2(ii))/sqrt(length(ii))]; 
    end
    p = polyfit(LagT,MSD_mean,1);
    yfit = polyval(p,LagT);
end