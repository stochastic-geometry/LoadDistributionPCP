% Author: Chiranjib Saha
% This script contains the exact expression of the variance of load for
% Thomas Cluster Process. See Eq. (8). 
% Returns Normalized Variance 
% Date Created: 09/24/2019
% Input: sigma = cluster size std, l_b = intensity of base station process, l_p =
% intensity of parent PPP 

function NormalizedVariance = VarianceThomas(sigma,l_b,l_p)
%m = 10;
 

% Second moment density
%rho = @(r) l_p^2*m^2 + l_p*m^2/(4*pi*sigma^2)*exp(-r.^2/(4*sigma^2));
d = @(x1,x2,theta) sqrt(x1.^2+x2.^2-2*cos(theta).*x1.*x2);


% I = @(x1,x2,theta)  exp(-area_union(x1,x2,d(x1,x2,theta))-d(x1,x2,theta).^2/(4*sigma^2*l_b)).*x1.*x2
% theta = 0.1;
% integral2(@(x1,x2)exp(-area_union(x1,x2,d(x1,x2,theta))).*...
%     exp(-d(x1,x2,theta).^2/(4*sigma^2*l_b)).*x1.*x2,...
%      0,100,0,100)
NormalizedVariance_term2 =  2/(4*pi*sigma^2)*integral3(@(theta,x1,x2)exp(-area_union(x1,x2,d(x1,x2,theta))).*...
    exp(-d(x1,x2,theta).^2/(4*sigma^2*l_b)).*x1.*x2,...
     0,pi,0,inf,0,@(theta,x1) x1,'reltol',1e-6,'abstol',1e-6);
%  1/(4*pi*sigma^2)*integral3(@(theta,x1,x2)exp(-area_union(x1,x2,d(x1,x2,theta))).*...
%     exp(-d(x1,x2,theta).^2/(4*sigma^2*l_b)).*x1.*x2,...
%      0,pi,0,inf,0,inf,'reltol',1e-4,'abstol',1e-4)
%  
 
 
 
 NormalizedVariance_term2 = NormalizedVariance_term2 * 2 * 2*pi /l_p; % Multiply constants for integration
 
 
%  NormalizedVariance_PPP = 4*pi* integral3(@(x1,x2,theta)exp(-area_union(x1,x2,d(x1,x2,theta))).*x1.*x2,...
%      0,inf,0,inf,0,pi,'reltol',1e-6,'abstol',1e-6);
 
NormalizedVariance_PPP = 8*pi* integral3(@(theta,x1,x2)exp(-area_union(x1,x2,d(x1,x2,theta))).*x1.*x2,...
    0,pi,0,inf,0,@(theta,x1)x1,'reltol',1e-6,'abstol',1e-6);
 
 
 
 
NormalizedVariance =  NormalizedVariance_PPP + NormalizedVariance_term2 -1;
end
% xlim = @(x)x;
%  2* integral3(@(theta,x1,x2)exp(-area_union(x1,x2,d(x1,x2,theta))).*x1.*x2,...
%      0,pi,0,inf,0,xlim,'reltol',1e-4,'abstol',1e-4)
 
 
