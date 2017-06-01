function [h2] = h2(r,t,Bo,rmax)

m = sqrt(Bo);

f3 = @(s)  (besselj(0, m*s, 0)*s)*(3*(1-t)/(1 - t + s^2/2)^2);
f4 = @(s)  (bessely(0, m*s, 0)*s)*(3*(1-t)/(1 - t + s^2/2)^2);
q3 = quad (f3, 0, r);
q4 = quad (f4, r, rmax);
h2 = pi/2*bessely(0, m*r, 0)*q3 + pi/2*besselj(0, m*r, 0)*q4;

