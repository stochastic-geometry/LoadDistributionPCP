% Author: Chiranjib Saha
% This script contains the exact expression of the variance of typical cell load for
% Matern Cluster Process. See Eq. (8). 
% Returns Normalized Variance 
% Date Created: 09/25/2019
% Input: R = cluster radius, l_b = intensity of base station process, l_p =
% intensity of parent PPP 
function NormalizedVariance = VarianceMatern(R,l_b,l_p)
%m = 10;
 % Second moment density
%rho = @(r) l_p^2*m^2 + l_p*m^2/(4*pi*sigma^2)*exp(-r.^2/(4*sigma^2));
d = @(x1,x2,theta) sqrt(x1.^2+x2.^2-2*cos(theta).*x1.*x2);

NormalizedVariance_term2 =  1/(pi*R^2)^2*integral3(@(x,theta,r)...
    exp(-l_b*area_union(x,sqrt(x.^2+r.^2+2*x.*r.*cos(theta)),r)).*...
       area_intersection(R,R,r).*x.*r,...
          0,inf,0,pi,0,2*R,'reltol',1e-6,'abstol',1e-6);
 
 
 
 NormalizedVariance_term2 = NormalizedVariance_term2 * 2 * 2*pi /l_p; % Multiply constants for integration
 
 
%  NormalizedVariance_PPP = 4*pi* integral3(@(x1,x2,theta)exp(-area_union(x1,x2,d(x1,x2,theta))).*x1.*x2,...
%      0,inf,0,inf,0,pi,'reltol',1e-6,'abstol',1e-6);
 
NormalizedVariance_PPP = 8*pi* integral3(@(theta,x1,x2)exp(-area_union(x1,x2,d(x1,x2,theta))).*x1.*x2,...
    0,pi,0,inf,0,@(theta,x1)x1,'reltol',1e-6,'abstol',1e-6);
 
NormalizedVariance =  NormalizedVariance_PPP + NormalizedVariance_term2;
end
 
 
