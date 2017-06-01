% This file construct the central difference operator matrix for a fourth order
% operator in the cylindrical coordinate
% operator = 1/r*(d(rA)/dr), A = dB/dr, B = 1/r*(d(rC)/dr), C = dh/dr
% boundary conditions given at h(J), C(0), B(J), A(0)


% define CFL = dt/dr*O(q0/r0)
 clear; clc; close all;

ca = 0.02;
fr = 0.02;
bo = ca/fr;
ma = 0.01;
tstop = 2.5;
R1 = 10;
T1 = 1;
dr = 0.5; dr2 = dr*dr; dr4 = dr2*dr2;
dt = 0.01;
po = 1;
pobar = -1;

Jmax = R1/dr;
Nmax = T1/dt;


% construct derivative matrices
DF = sparse(Jmax, Jmax);
DA = sparse(Jmax, Jmax);
DP = sparse(Jmax, Jmax);
DN = sparse(Jmax, Jmax);
J = Jmax;
for i = 2:Jmax
	

		%/ backward-gradient
		DA(i, i) =  1 + 0.5/(i-1);
		DA(i,i-1) = -1 + 0.5/(i-1);

		% forward
		DF(i, i) = -1;
		if (i < Jmax)
			DF(i,i+1) = 1;
		end
	
		% positive convection speed
		DP(i,i) = 1;
		DP(i,i-1) = 1/i-1;

		% negative convection speed
		DN(i,i)= -1;
		if (i < Jmax)
			DN(i,i+1) = 1+1/i;
		end
		

end
	DA(1,1) =  4; 
	DP(1,1) =  2; DP(1,2) = -2; 
	DN(1,1) = -2; DN(1,2) =  2;
	DF(1,1) = -1; DF(1,2) =  1;

L = sparse(Jmax,Jmax);
for i = 3:Jmax-2
	L(i, i-2) = (i-1/2)*(i-3/2)/(i)/(i-1);
	L(i, i-1) = -4*(1-0.5/i);
	L(i, i  ) = 4 + (i+1/2)*(i+1/2)/(i)/(i+1) + (i-1/2)*(i-1/2)/(i)/(i-1);
	L(i, i+1) = -4*(1+0.5/i);
	L(i, i+2) = (i+1/2)*(i+3/2)/(i)/(i+1);
end

i = 1;
L(i, i  ) = 5.5;
L(i, i+1) = -8;
L(i, i+2) = 5/2;

i = 2;
L(i, i-1) = -4;
L(i, i  ) = 6.125;
L(i, i+1) = -4;
L(i, i+2) = 15/8;

i = Jmax-1;
	L(i, i-2) = (i-1/2)*(i-3/2)/(i)/(i-1);
	L(i, i-1) = -4*(1-0.5/i);
	L(i, i  ) = 4 + (i+1/2)*(i+1/2)/(i)/(i+1) + (i-1/2)*(i-1/2)/(i)/(i-1);
	L(i, i+1) = -4*(1+0.5/i);

i = Jmax;
	L(i, i-2) = (i-1/2)*(i-3/2)/(i)/(i-1);
	L(i, i-1) = -4*(1-0.5/i);
	L(i, i  ) = 4 + (i+1/2)*(i+1/2)/(i)/(i+1) + ...
			(i-1/2)*(i-1/2)/(i)/(i-1) - (i+1/2)*(i+3/2)/i/(i+1);




	
	
% initialize variables:
%rvec = 0:1:(Jmax-1);
%rvec = rvec';
%h0 = 1+ dr2*rvec.*rvec/2;
%q0 = zeros(size(rvec));
%
%alpha = dt/2/dr;
%I = eye(Jmax);
%Z = zeros(Jmax, Jmax);
%D4 = L/dr4-bo*bo*I;
%
%DD = DA*DF/dr2;
%b2 = ones(size(rvec))*po*2;
%plot(DD\b2)


DD = DA*DF;
DD(1,1) = -2; DD(1,2) = 2;
DDDD = DD*DD-L;
%figure; hold on;
%
%for i = 1:1:Nmax 
%time = i*dt;
%h1J = 1+R1*R1/2-time;
%
%b2 = bo*pobar*ones(size(rvec));
%j = Jmax-1;
%b2(j) = b2(j) - h1J/dr4*(1+1/2/j)*(j+1.5)/(j+1);
%j = Jmax;
%%b2(j) = b2(j) - ((2*po-bo*h1J)/dr2+2/dr4*h1J)*(j+1.5)/(j+1) + 4/dr4*h1J;
%b2(j) = b2(j) - ((2*po-bo*h1J)/dr2+2/dr4*h1J)*(j+1.5)/(j+1) + 4/dr4*h1J;
%
%
%hnew = D4\b2;
%
%
%%plot(hnew(1:floor(end*0.7)));
%plot(hnew);
%pause(0.05);
%
%end
