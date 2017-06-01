close all; clear; 

Ca = 0.01;
Bo = 0.1;

rmin = 0;
rmax = 10;
dr = 0.1
Jmax = (rmax - rmin)/dr+1;
dt = 0.000000001
time = 0.0;

rvec = dr/2:dr:rmax+dr/2;
rvec = 0:dr:rmax;
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

%h1vec = zeros(size(rvec));
%for k = 1:length(rvec)
%	r = rvec(k);
%	[h1 h2] = h1h2(r,time,Bo,rmax*5);
%	h1vec(k) = h1;
%end

% implicit
Z = sparse(Jmax,Jmax);
L11 = eye(Jmax)/dt;
L12 = Z;
L13 = DA*diag(1+rvec.*rvec/2-time);
L11p = L11;
L12p = Z;
L13p = Z;
L21 = DD; 
L22 = -Bo*eye(Jmax);
L23 = 6*Ca*DA;
L31 = -Bo*eye(Jmax);
L32 = DD;
L33 = Z;

%% CN
%Z = sparse(Jmax,Jmax);
%L11 = eye(Jmax)/dt;
%L12 = Z;
%L13 = DA*diag(1+rvec.*rvec/2-time)/2;
%L11p = L11;
%L12p = -L12;
%L13p = -L13;
%L21 = DD; 
%L22 = -Bo*eye(Jmax);
%L23 = 6*Ca*DA;
%L31 = -Bo*eye(Jmax);
%L32 = DD;
%L33 = Z;



TRPLHS = [L11 L12 L13; L21 L22 L23; L31 L32 L33];
TRPRHS = [L11p L12p L13p; Z Z Z; Z Z Z];
Lall = inv(TRPLHS)*TRPRHS;

maxeig = max(abs(real(eig(Lall))))
mineig = min(abs(real(eig(Lall))))
