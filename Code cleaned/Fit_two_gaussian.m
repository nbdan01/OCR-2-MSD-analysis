function y = Fit_two_gaussian(param,x)

 a_1 = param(1);
 x_1 = param(2);
 s_1 = param(3);
 a_2 = param(4);
 x_2 = param(5);
 s_2 = param(6);


  y = a_1.*exp(-(x-x_1).^2./(2*s_1.^2))+(a_2).*exp(-(x-x_2).^2./(2*s_2.^2));
%     y = a_1./(1+((x-x_1)./s_1).^2)+a_2./(1+((x-x_2)./s_2).^2)
end