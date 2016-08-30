%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Script to calculate the signature function %%%%%%
%%%%%%                        2D                  %%%%%%
%%%%%% V(r) = -(r^2)/4 d/dr(1/r d <delv^2>/dr)    %%%%%%
%%%%%%                                            %%%%%%
%%%%%% Dhruv Balwada ; July 14 2016               %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
load S2.mat 
% s2ll = nanmoving_average(s2ll,16);
% s2tt = nanmoving_average(s2tt,16);
% s2rr = nanmoving_average(s2rr,16);
% s2dd = nanmoving_average(s2dd,16);

%% polynomial fit 
% 
% [pll,~,mu] = polyfit(dist_axis(15:42), s2ll(15:42)', 4);
r = [dist_axis(14):1:dist_axis(40)];
% s2ll_fit = polyval(pll,r,[],mu);
[prr,~,mu] = polyfit(dist_axis(14:40), s2rr(14:40), 5);
s2rr_fit = polyval(prr,r,[],mu);
[pdd,~,mu] = polyfit(dist_axis(14:40), s2dd(14:40), 5);
s2dd_fit = polyval(pdd,r,[],mu);

%
close
figure,
loglog(dist_axis/1000, s2dd)
hold all
loglog(r/1000, s2dd_fit)
loglog(dist_axis/1000, s2rr)
loglog(r/1000, s2rr_fit)
axis([10^-2 10^3 10^-6 10^0])
%
% dslldr = nan(length(dist_axis),1) ; 
% dsttdr = nan(length(dist_axis),1) ; 
% dsrrdr = nan(length(dist_axis),1) ; 
% dsdddr = nan(length(dist_axis),1) ; 
% 
% r      = nan(length(dist_axis),1) ; 
% 
% for ii = 2 : length(dist_axis)-1
%     dr = dist_axis(ii+1) - dist_axis(ii);
%     dslldr(ii) = (s2ll(ii+1) - s2ll(ii))/dr;
%     dsttdr(ii) = (s2tt(ii+1) - s2tt(ii))/dr;
%     dsrrdr(ii) = (s2rr(ii+1) - s2rr(ii))/dr;
%     dsdddr(ii) = (s2dd(ii+1) - s2dd(ii))/dr;
%     r(ii) = 0.5*(dist_axis(ii+1) + dist_axis(ii)); 
% end

dslldr = nan(length(r),1) ; 
dsttdr = nan(length(r),1) ; 
dsrrdr = nan(length(r),1) ; 
dsdddr = nan(length(r),1) ; 

r2 = nan*r;
for ii = 2 : length(r)-1
    dr = r(ii+1) - r(ii);
%     dslldr(ii) = (s2ll(ii+1) - s2ll(ii))/dr;
%     dsttdr(ii) = (s2tt(ii+1) - s2tt(ii))/dr;
    dsrrdr(ii) = (s2rr_fit(ii+1) - s2rr_fit(ii))/dr;
    dsdddr(ii) = (s2dd_fit(ii+1) - s2dd_fit(ii))/dr;
    r2(ii) = 0.5*(r(ii+1) + r(ii)); 
end

Vll   = nan(length(r),1) ; 
Vtt   = nan(length(r),1) ; 
Vdd   = nan(length(r),1) ; 
Vrr   = nan(length(r),1) ; 

for ii = 2 : length(r)-1
    dr = r2(ii+1) - r2(ii); 
    
%     Vll(ii) = -r(ii)^2/4*(1/r2(ii+1)*dslldr(ii+1) ...
%         - 1/r2(ii)*dslldr(ii))/dr;
%     Vtt(ii) = -r(ii)^2/4*(1/r2(ii+1)*dsttdr(ii+1) ...
%         - 1/r2(ii)*dsttdr(ii))/dr;
    Vdd(ii) = -r(ii)^2/4*(1/r2(ii+1)*dsdddr(ii+1) ...
        - 1/r2(ii)*dsdddr(ii))/dr;
    Vrr(ii) = -r(ii)^2/4*(1/r2(ii+1)*dsrrdr(ii+1) ...
        - 1/r2(ii)*dsrrdr(ii))/dr;
end

%


figure
loglog(r/1000, Vdd.*r')
hold all
loglog(r/1000, Vrr.*r')
axis([10^-2 10^3 10^-6 10^0])

%%
% close all
% figure
% loglog(dist_axis/1000, Vll.*dist_axis')
% hold all
% loglog(dist_axis/1000, Vtt.*dist_axis')
% loglog(dist_axis/1000, Vdd.*dist_axis')
% loglog(dist_axis/1000, Vrr.*dist_axis')
% axis([10^-2 1000 10^-6 10^-1])
% 
% figure
% loglog(dist_axis/1000, s2ll)
% hold all
% loglog(dist_axis/1000, s2tt)
% loglog(dist_axis/1000, s2dd)
% loglog(dist_axis/1000, s2rr)
% axis([10^-2 1000 10^-5 10^0])