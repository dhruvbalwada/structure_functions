clear all
close all

load S2.mat

%%
figure(25) , hold all 
semilogx(dist_axis./1000, (s2dd./s2rr).^0.5,'linewidth',3)
axis([10^0 10^3 0 3])
set(gca,'fontsize',16)
xlabel('r (km)')
ylabel(' Divergent/ Rotational')