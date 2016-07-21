%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Script to calculate the signature function %%%%%%
%%%%%%                        2D                  %%%%%%
%%%%%% V(r) = -(r^2)/4 d/dr(1/r d <delv^2>/dr)    %%%%%%
%%%%%%                                            %%%%%%
%%%%%% Dhruv Balwada ; July 14 2016               %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
load S2.mat 
% s2ll = nanmoving_average(s2ll,6);
% s2tt = nanmoving_average(s2tt,6);
% s2rr = nanmoving_average(s2rr,6);
% s2dd = nanmoving_average(s2dd,6);

dslldr = nan(length(dist_axis),1) ; 
dsttdr = nan(length(dist_axis),1) ; 
dsrrdr = nan(length(dist_axis),1) ; 
dsdddr = nan(length(dist_axis),1) ; 

r      = nan(length(dist_axis),1) ; 

for ii = 2 : length(dist_axis)-1
    dr = dist_axis(ii+1) - dist_axis(ii);
    dslldr(ii) = (s2ll(ii+1) - s2ll(ii))/dr;
    dsttdr(ii) = (s2tt(ii+1) - s2tt(ii))/dr;
    dsrrdr(ii) = (s2rr(ii+1) - s2rr(ii))/dr;
    dsdddr(ii) = (s2dd(ii+1) - s2dd(ii))/dr;
    r(ii) = 0.5*(dist_axis(ii+1) + dist_axis(ii)); 
end


Vll   = nan(length(dist_axis),1) ; 
Vtt   = nan(length(dist_axis),1) ; 
Vdd   = nan(length(dist_axis),1) ; 
Vrr   = nan(length(dist_axis),1) ; 

for ii = 2 : length(dist_axis)-1
    dr = r(ii+1) - r(ii); 
    
    Vll(ii) = -dist_axis(ii)^2/4*(1/r(ii+1)*dslldr(ii+1) ...
        - 1/r(ii)*dslldr(ii))/dr;
    Vtt(ii) = -dist_axis(ii)^2/4*(1/r(ii+1)*dsttdr(ii+1) ...
        - 1/r(ii)*dsttdr(ii))/dr;
    Vdd(ii) = -dist_axis(ii)^2/4*(1/r(ii+1)*dsdddr(ii+1) ...
        - 1/r(ii)*dsdddr(ii))/dr;
    Vrr(ii) = -dist_axis(ii)^2/4*(1/r(ii+1)*dsrrdr(ii+1) ...
        - 1/r(ii)*dsrrdr(ii))/dr;
end

%% 
close all
figure
loglog(dist_axis/1000, Vll.*dist_axis')
hold all
loglog(dist_axis/1000, Vtt.*dist_axis')
loglog(dist_axis/1000, Vdd.*dist_axis')
loglog(dist_axis/1000, Vrr.*dist_axis')
axis([10^-2 1000 10^-6 10^-1])

figure
loglog(dist_axis/1000, s2ll)
hold all
loglog(dist_axis/1000, s2tt)
loglog(dist_axis/1000, s2dd)
loglog(dist_axis/1000, s2rr)
axis([10^-2 1000 10^-5 10^0])