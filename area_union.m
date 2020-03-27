function A = area_union(r1,r2,d) 
% Author: Chiranjib Saha
% Computes the area of union of two circles of radius r1 and r2 and
% the distance between the centers d. 
A = pi*(r1.^2+r2.^2)- area_intersection(r1,r2,d);
end