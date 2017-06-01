clear; clc;%close all;
tvec = [0.1 0.6 0.99];
Bo = 1.0;
Ca = 0.1;
rmax = 10;
rvec = 0.0:0.05:rmax;

nt = length(tvec);
nr = length(rvec);

h1vec = zeros(nt,nr);
h2vec = zeros(nt,nr);
h1comp = zeros(nt,nr);
h2comp = zeros(nt,nr);

% compute the coefficient for h2
m = sqrt(Bo);
for tcount = 1:nt
	t = tvec(tcount);
	for k = 1:nr
		r = rvec(k);
		[h1 h2] = h1h2(r,t,Bo,rmax);
		h1vec(tcount,k) = h1;
		h2vec(tcount,k) = h2;
	end
	h1comp(tcount,:) = Ca*(3/2*log(rvec.^2/2./(1-t+rvec.^2/2))+3*besselk(0,m*rvec,0));
	h2comp(tcount,:) = -1-rvec.^2/2+t+Ca*(-3/2*log(rvec.^2/2./(1-t+rvec.^2/2))+3*pi/2*bessely(0,m*rvec,0));
end
h1tot = Ca*h1vec;
temp = [1-tvec(1)+rvec.^2/2;1-tvec(2)+rvec.^2/2;1-tvec(3)+rvec.^2/2];
h2tot = (-temp + Ca*h2vec);


figure; 

subplot(1,3,1); hold on;
plot(rvec, h1tot (1,:), 'b- ', 'linewidth',9);
plot(rvec, h2tot (1,:), 'k- ', 'linewidth',9);
plot(rvec, h1comp(1,:), 'ro ', 'linewidth',9);
plot(rvec, h2comp(1,:), 'g+ ', 'linewidth',9);
axis([0 2 -2 1.0])
set(gca,"fontsize",28,"linewidth",3)
xlabel('r','fontsize',32);
ylabel('h','fontsize',32);
title('t = 0.1','fontsize',32);
h = legend('h1 VOP','h2 VOP','h1 MAE','h2 MAE',"location",'north');
set(h,'fontsize',18,'color','none')
legend boxoff


subplot(1,3,2); hold on;
plot(rvec, h1tot (2,:), 'b- ', 'linewidth',9);
plot(rvec, h2tot (2,:), 'k- ', 'linewidth',9);
plot(rvec, h1comp(2,:), 'ro ', 'linewidth',9);
plot(rvec, h2comp(2,:), 'g+ ', 'linewidth',9);
xlabel('r','fontsize',32);
ylabel('h','fontsize',32);
title('t = 0.6','fontsize',32);
axis([0 2 -2 1.0])
set(gca,"fontsize",28,"linewidth",3)


subplot(1,3,3); hold on;
plot(rvec, h1tot (3,:), 'b- ', 'linewidth',9);
plot(rvec, h2tot (3,:), 'k- ', 'linewidth',9);
plot(rvec, h1comp(3,:), 'ro ', 'linewidth',9);
plot(rvec, h2comp(3,:), 'g+ ', 'linewidth',9);
xlabel('r','fontsize',32);
ylabel('h','fontsize',32);
title('t = 0.99','fontsize',32);
axis([0 2 -2 1.0])
set(gca,"fontsize",28,"linewidth",3)


%hold on
%plot(rvec, h2vec(1,:), 'g--', 'linewidth',4);
%plot(rvec, h2vec(2,:), 'm--', 'linewidth',4);
%plot(rvec, h2vec(3,:), 'c--', 'linewidth',4);

