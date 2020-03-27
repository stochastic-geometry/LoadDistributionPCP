function A = area_intersection(r1,r2,d)
% Author: Chiranjib Saha
% Computes the area of intersection of two circles of radius r1 and r2 and
% the distance between the centers d. 
   %(d+r1+r2).*(d+r1-r2).*(d-r1+r2).*(-d+r1+r2)
   t = sqrt((d+r1+r2).*(d+r1-r2).*(d-r1+r2).*(-d+r1+r2));
   A = r1.^2.*atan2(real(t),d.^2+r1.^2-r2.^2)+r2.^2.*atan2(real(t),d.^2-r1.^2+r2.^2)-t/2;
   % Corner cases
   A(d>r1+r2) = 0;
   A_min = pi*min(r1,r2).^2; 
   A(d<abs(r1-r2)) = A_min(d<abs(r1-r2)); 
end