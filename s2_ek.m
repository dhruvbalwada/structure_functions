% Find intersection point

%% Spectra
k = 2*pi./[0:100:10^6];

Lx = 10*10^3;

a = 1;
c = a*0.24*(pi/Lx)^(4/3);

figure
loglog(k, 0.12*a*k.^(-5/3))
hold all
loglog(k, c/2*k.^(-3))

%% Structure function

r = [0:100:10^6];

b = 0.0001;
figure
loglog(r, a*r.^(2/3))
hold all
loglog(r,  c*(r.^2).*log(r))

%% Using Hankel transforms
clear
k = 2*pi./10.^[6:-0.001:4];

K_0 = min(k);
K_1 = max(k);

dk = diff(k);
dk = [dk(1) dk];

alpha1 = [5/3];
alpha2 = [2.9];
L_0 = 2*pi/K_0;
L_1 = 2*pi/K_1;

a_0 = floor(log(L_0)/log(1.2));
a_1 = ceil(log(L_1)/log(1.2));

r = 1.2.^[a_1:a_0+1];


E1 = 1;
E2 = 20;

Ek1 = E1*(k/K_0).^(-alpha1);
Ek2 = E2*(k/K_0).^(-alpha2);

% Evaluate with J_1
% for ii = 1:length(r)
%     S1(ii) = 2*sum(dk.*Ek1.*(1- 2*besselj(1,k*r(ii))./k/r(ii)));
%     S2(ii) = 2*sum(dk.*Ek2.*(1- 2*besselj(1,k*r(ii))./k/r(ii)));
% end
% S1th_com = E1*K_0^alpha1*(2*pi)^(3-alpha1)/4/(3-alpha1)*r.^(alpha1-1) - ...
%     E1*r.^2*K_0^3/4/(3-alpha1) + 2*E1*K_0^alpha1*K_1^(1-alpha1)/(1-alpha1) - ...
%     2*E1*K_0^alpha1*(2*pi)^(1-alpha1)/(1-alpha1);
% S2th_com = E2*K_0^alpha2*(2*pi)^(3-alpha2)/4/(3-alpha2)*r.^(alpha2-1) - ...
%     E2*r.^2*K_0^3/4/(3-alpha2) + 2*E2*K_0^alpha2*K_1^(1-alpha2)/(1-alpha2) - ...
%     2*E2*K_0^alpha2*(2*pi)^(1-alpha2)/(1-alpha2);


% Evaluate with J_0 (needs fixing)
for ii = 1:length(r)
    S1(ii) = 2*sum(dk.*Ek1.*(1- besselj(0,k*r(ii))));
    S2(ii) = 2*sum(dk.*Ek2.*(1- besselj(0,k*r(ii))));
end
% limits of 2pi/r
S1th_com =  E1*K_0^alpha1*(2*pi)^(3-alpha1)/4/(3-alpha1)*r.^(alpha1-1) - ...
            E1*r.^2*K_0^3/4/(3-alpha1) + ...
            E1*K_0^alpha1*K_1^(1-alpha1)/(1-alpha1) - ...
            E1*K_0^alpha1*(2*pi)^(1-alpha1)/(1-alpha1);

S2th_com =  E2*K_0^alpha2*(2*pi)^(3-alpha2)/4/(3-alpha2)*r.^(alpha2-1) - ...
            E2*r.^2*K_0^3/4/(3-alpha2) + ...
            E2*K_0^alpha2*K_1^(1-alpha2)/(1-alpha2) - ...
            E2*K_0^alpha2*(2*pi)^(1-alpha2)/(1-alpha2);

% S1_th = 2*E1*(K_0^alpha1)*r.^(alpha1-1)/(alpha1-1)/(2*pi)^(alpha1-1);
% S2_th = E2*(K_0^3)*r.^2/(alpha2-3)/4;
% S2_th = 2*E2*(K_0^alpha2)*r.^(alpha2-1)/(alpha2-1)/(2*pi)^(alpha2-1);

% limits of 1/r.
% S1th_com = E1*K_0^alpha1/4/(3-alpha1)*r.^(alpha1-1) - ...
%     E1*r.^2*K_0^3/4/(3-alpha1) + E1*K_0^alpha1*K_1^(1-alpha1)/(1-alpha1) - ...
%     E1*K_0^alpha1/(1-alpha1);
% S2th_com = E2*K_0^alpha2*(2*pi)^(3-alpha2)/4/(3-alpha2)*r.^(alpha2-1) - ...
%     E2*r.^2*K_0^3/4/(3-alpha2) + E2*K_0^alpha2*K_1^(1-alpha2)/(1-alpha2) - ...
%     E2*K_0^alpha2*(2*pi)^(1-alpha2)/(1-alpha2);

 close all
%
figure
subplot(1,2,1)
loglog(r, S1,'.-')
hold all
loglog(r, S2,'.-')
loglog(r, 2*S1th_com,'s-')
loglog(r, 2*(S2th_com),'s-')
legend('S_{5/3}', ['S_' num2str(alpha2)], 'S^{th}_{5/3}', ['S^{th}_' num2str(alpha2)] ,'location','best')
% axis([min(r) max(r) 10^-7 10^-4])
% set(gca,'fontsize',16)


%
subplot(1,2,2)
%figure
loglog(k,Ek1,'.-')
hold all
loglog(k,Ek2,'.-')
axis([K_0 K_1 10^-5 10^2])
legend('E_{5/3}', ['E_{' num2str(alpha2) '}'])
% set(gca,'fontsize',16)

%%
k = 2*pi./10.^[8:-0.1:4];

K_0 = min(k);
K_1 = max(k);
syms r
a1 = 5/3;
a2mat = [1.5:0.1:2.9];
% 
% E1 = 1;
% E2 = 10000;
% E1divE2 = E1/E2;
E1divE2min = (K_1/K_0)^(1./(a1 - 2.9)); 
powmin = floor(log10(E1divE2min));
E1divE2mat = 10.^[powmin:0.05:0];

a1 = 5/3;
a2mat = [1.5:0.1:2.9];
clear solL

% eqn = (E1divE2*K_0^3/(a1 - 3) - K_0^3/4/(a2-3))*r^2 + (E1divE2*K_0^a1*(2*pi)^(1-a1)*(pi^2/(3-a1) - 2/(1-a1)))*r^(a1-1) ...
%     - (K_0^a2*(2*pi)^(1-a2)*(pi^2/(3-a2) - 2/(1-a2)))*r^(a2-1) + 2*E1divE2*K_0^a1*K_1^(1-a1)/(1-a1) - 2*K_0^a2*K_1^(1-a2)/(1-a2) == 0;
for jj = 1:length(E1divE2mat)
    for ii = 1:length(a2mat)
        a2 = a2mat(ii);
        E1divE2 = E1divE2mat(jj);
%         eqn = (E1divE2*K_0^3/(a1 - 3) - K_0^3/4/(a2-3))*r^2 + (E1divE2*K_0^a1*(2*pi)^(1-a1)*(pi^2/(3-a1) - 2/(1-a1)))*r^(a1-1) ...
%             - (K_0^a2*(2*pi)^(1-a2)*(pi^2/(3-a2) - 2/(1-a2)))*r^(a2-1) + 2*E1divE2*K_0^a1*K_1^(1-a1)/(1-a1) - 2*K_0^a2*K_1^(1-a2)/(1-a2) == 0;
        eqn = (E1divE2*K_0^3/(a1 - 3) - K_0^3/4/(a2-3))*r^2 + (E1divE2*K_0^a1*(2*pi)^(1-a1)*(pi^2/(3-a1) - 1/(1-a1)))*r^(a1-1) ...
            - (K_0^a2*(2*pi)^(1-a2)*(pi^2/(3-a2) - 1/(1-a2)))*r^(a2-1) + 1*E1divE2*K_0^a1*K_1^(1-a1)/(1-a1) - 1*K_0^a2*K_1^(1-a2)/(1-a2) == 0;
%         
        
        solLk =  K_0*(E1divE2).^(1./(a1-a2));
        
        if solLk < K_1 & solLk>K_0
            solx = vpasolve(eqn,r,pi./solLk);
            if ~isempty(solx)
                solL(1,ii,jj) = pi./solLk;
                solL(2,ii,jj) = solx;
            else
                solL(1,ii,jj) = NaN;
                solL(2,ii,jj) = NaN;
            end
        else
            solL(1,ii,jj) = NaN;
            solL(2,ii,jj) = NaN;
        end
    end
end


ratio = squeeze(solL(1,:,:)./solL(2,:,:));
Lk = squeeze(solL(1,:,:))/(pi./K_0); 
% pcolor(E1divE2mat, a2mat, ratio)
% figure, pcolor(E1divE2mat, a2mat, ratio)

%%
figure, 
% set(gca,'fontsize',16)
contourf(log10(E1divE2mat), a2mat, ratio,'edgecolor','none')
hold all 
[h,c] = contour(log10(E1divE2mat), a2mat, ratio,[1 2],'edgecolor','r','linewidth',2)
clabel(h,c)
[h, c] = contour(log10(E1divE2mat), a2mat, log10(Lk));
clabel(h,c)
caxis([0.5 3])
h= colorbar;
ylabel(h,'\pi / K_{i} / r_{i}')
xlabel('log_{10}(E_\alpha/E_\beta)') 
ylabel('\beta')
title(['K_1/K_0 = ' num2str(K_1/K_0) '; \alpha = 5/3']) 
set(gca,'fontsize',16)
%% How the intersection point of the two spectra varies
clear
E1divE2 = 10.^[-10:0];
a1mina2 = -[0.2:0.1:2.9];
K_0 = 10^-10;

for i =1:length(a1mina2)
    Kintersect(:,i) = K_0*(E1divE2).^(1./a1mina2(i));
end

close all
figure
pcolor(log10(E1divE2), a1mina2, log10(pi./(Kintersect')))
xlabel('E1/E2')
ylabel('\alpha 1 - \alpha 2')
caxis([-1 5])