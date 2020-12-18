function [yfitted,yfitted2D]= gaussian2d(beta,X,ima)

amp     = beta(1);
x0      = beta(2);
y0      = beta(3);
wx      = beta(4);
wy      = beta(5);
Offset  = beta(6);

xx = X(1:length(X)/2);
yy = X(length(X)/2+1 : end);

yfitted = amp*exp(-log(2)*(xx-x0).^2 / (wx^2)).*exp(-log(2)*(yy-y0).^2 / (wy^2))+Offset;
% yfitted2D = 0*ima;
S = sparse(yy,xx,yfitted,size(ima,2),size(ima,1));
yfitted2D =  full(S);