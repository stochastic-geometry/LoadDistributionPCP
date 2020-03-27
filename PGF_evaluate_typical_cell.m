function PGF = PGF_evaluate_typical_cell(theta,clustersize,m,l_p,l_b,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Chiranjib Saha
% Computes the PGF of the typical cell load
% Input:       theta = variable of the PGF, 
%        clustersize = cluster size which is sigma for TCP and R for MCP
%                 m  = number of users per cluster
%                l_p = intensity of parent PPP 
%                l_b = intensity of base station process
%               type = 'Thomas' for TCP, 'Matern' for MCP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_bar = m;
%% Nakagami 
m = 3.5;%4.5;
w = 1;%4.5/3.5;
f_nakagami = @(x) ((2*m^m)/(gamma(m)*w^m))*x.^(2*m-1).*exp(-((m/w).*x.^2));
cell_radius_PDF = @(r) sqrt(pi*l_b).*f_nakagami(sqrt(pi*l_b)*r); 
if strcmp(type,'Thomas')
 sigma = clustersize;
 f = @(r) integral( @(v) (1-exp(-m_bar.* (1-theta).*(1-marcumq(v./sigma,r./sigma)))).*v,0,20*sigma,'reltol',1e-3,'abstol',1e-3,'arrayvalued',true);
 PGF =integral(@(r)...
     cell_radius_PDF(r) .* ...
       exp(-2 *pi* l_p .* f(r)),...
          0,450,'reltol',1e-3,'abstol',1e-3);
elseif strcmp(type,'Matern')
   R = clustersize; 
   G = @(r,v) 1/R^2*((min(r,max(R-v,0))).^2 + 2/pi*integral(@(u) u.*real(acos((u.^2+v.^2-R^2)./(2*u.*v))), min(r,abs(R-v)),min(r,R+v)));
   G_arr = @(r,v) arrayfun(@(v) G(r,v),v);
   f = @(r) integral(@(v) (1-exp(-m_bar.* (1-theta).*G_arr(r,v))).*v,0,20*R,'reltol',1e-3,'abstol',1e-3);
   f_arr = @(r) arrayfun(@(r) f(r),r);
   PGF =integral(@(r)...
     cell_radius_PDF(r) .* ...
       exp(-2 *pi* l_p .* f_arr(r)),...
          0,20*R,'reltol',1e-3,'abstol',1e-3);
end

