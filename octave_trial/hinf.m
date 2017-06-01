close all; clear; clc

ca = 0.02;
fr = 0.02;
bo = 1;
t = 0.9;

R1 = 10;
rspan = 1;
rmin = R1-rspan;
rmax = R1+rspan;
rmin = 0.0;
rmax = R1;
Jmax = 10001;
dr = (rmax - rmin)/(Jmax-1);

rvec = rmin:dr:rmax;
rvec = rvec';


% construct derivative matrices
DF = sparse(Jmax, Jmax);
DA = sparse(Jmax, Jmax);
DP = sparse(Jmax, Jmax);
DN = sparse(Jmax, Jmax);
DD = sparse(Jmax, Jmax);
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
DA = DA/dr;
DP = DP/dr;
DN = DN/dr;
DF = DF/dr;

DD = DA*DF;

L33 = DD;

for jj = 1:Jmax
L33(jj,jj) = L33(jj,jj) + bo;
end
pvec = (1-t)./(1+rvec.^2/2-t).^2;

h2fd = L33\(3*pvec);


Jint = 100;
h2int = zeros(1,Jint);
rint = linspace(0,R1,Jint);
rmax = 100;
for rr = 1:Jint
	r = rint(rr);
	h2int(rr) = h2(r,t,bo,rmax);
end




figure;hold on;
plot(rvec,h2fd,'r--','linewidth',4);
plot(rint,h2int,'k-','linewidth',4);
legend('FINITE DIFFERENCE','INTEGRATION')




