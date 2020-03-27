function [ordered_pdf N]= pgf_inversion(N,clustersize,m,l_p,l_b,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Chiranjib Saha
% Computes the PD F of the typical cell load
% Input:       theta = variable of the PGF, 
%        clustersize = cluster size which is sigma for TCP and R for MCP
%                 m  = number of users per cluster
%                l_p = intensity of parent PPP 
%                l_b = intensity of base station process
%               type = 'Thomas' for TCP, 'Matern' for MCP
%% specifying lambda
%lambda = 4;
%P = @(z) exp(lambda*(z-1)); % Z-transform
%%%%%%%%%%%%%%%%%%%%%
R = 1; % Radius of region of covergence
for i=0:N-1
  i
  X(i+1) =  PGF_evaluate_typical_cell(R.*exp(j*2*pi*i/N),clustersize,m,l_p,l_b,type);  %P(R.*exp(j*2*pi*[0:N-1]./N));
end
pdf = ifft(X,N);
sum(real(pdf))
ordered_pdf = [real(pdf(1)) real(pdf(end:-1:2))];
end