%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 Check  Istropy                                      %%
%%%                 Dhruv Balwada                                       %%
%%%                  July 14 2016                                       %%
%%% Purpose: Test different measures to see if the data set is isotropic%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

load glad_traj.mat

%% Horizontal structure of the mean flow and energy

dx = 0.5;

X_grid = [nanmin(nanmin(X)):dx:nanmax(nanmax(X))];
Y_grid = [nanmin(nanmin(Y)):dx:nanmax(nanmax(Y))];
X_axis = 0.5*(X_grid(1:end-1) + X_grid(2:end));
Y_axis = 0.5*(Y_grid(1:end-1) + Y_grid(2:end));

Umean = nan(length(X_axis), length(Y_axis));
Vmean = nan(length(X_axis), length(Y_axis));
Ueke  = nan(length(X_axis), length(Y_axis));
Veke  = nan(length(X_axis), length(Y_axis));
nids  = Ueke;

for ii = 1 : length(X_axis)
    for jj = 1 : length(Y_axis)
        id = find(X>=X_grid(ii) & X<X_grid(ii+1) ...
            & Y>=Y_grid(jj) & Y<Y_grid(jj+1));
        
        nids(ii,jj) = length(id);
        
        Updf = U(id);
        Vpdf = V(id);
        
        Umean(ii,jj) = nanmean(Updf);
        Vmean(ii,jj) = nanmean(Vpdf);
        
        Ueke (ii,jj) = nanvar(Updf);
        Veke (ii,jj) = nanvar(Vpdf);
        
        % variance elipses
        clear XYd temp
        XYd(:,1) = Updf - Umean(ii,jj);
        XYd(:,2) = Vpdf - Vmean(ii,jj);
        
        idt = find(~isnan(XYd(:,1)));
        
        for h=1:2
            temp(:,1)=XYd(idt,1);
            temp(:,2)=XYd(idt,2);
        end
        
        XYd=temp;
        % find covariance
        Cv=1/length(XYd)*XYd'*XYd;
        
        %--- Find eigenvectors (Ev) and eigenvalues (El): ---%
        % find eigen vectors and values
        [Ev El]=eig(Cv);
        e1  =Ev(:,1)'; e2  =Ev(:,2)';
        klmb(1)=abs(El(1,1));  klmb(2)=abs(El(2,2));
        % take absolutes because sometimes we get very small negative values
        
        %--- Major axis (Applied multivariate statistical analysis, p. 463): ---%
        if (klmb(1) > klmb(2))
            kmj=klmb(1)^0.5*e1; % scaling using ||e1||=sqrt(lambda)
            kmn=klmb(2)^0.5*e2;
        else
            kmj=klmb(2)^0.5*e2;
            kmn=klmb(1)^0.5*e1;
        end
        
        % Plot ellipse:
        tt=(0:0.1:2*pi+0.1);
        aa=sqrt(kmj*kmj');
        bb=sqrt(kmn*kmn');
        A=[aa*cos(tt);bb*sin(tt)];
        alf=atan2(kmj(2),kmj(1));
        Rot=[cos(-alf), sin(-alf); -sin(-alf), cos(-alf)];
        fac_el=0.5;
        A2=Rot*A*fac_el;             % in units of velocity scaled up
        Ell(ii,jj,1:2,:)=[A2(1,:)/cosd(Y_axis(jj))+X_axis(ii);A2(2,:)+Y_axis(jj)];
    end
end

%% Velocity contours 
figure
m_proj('Miller','lon',[-95 -82.5], 'lat',[22 31])
hold all
m_contourf(X_axis, Y_axis, Umean', [-10 -1:0.2:1])
m_grid('linestyle','none','xaxisloc','top','yaxisloc','right')
m_coast('patch',[0 0 0]);
caxis([-1 1])

figure
m_proj('Miller','lon',[-95 -82.5], 'lat',[22 31])
hold all
m_contourf(X_axis, Y_axis, Vmean', [-10 -1:0.2:1])
m_grid('linestyle','none','xaxisloc','top','yaxisloc','right')
m_coast('patch',[0 0 0]);
caxis([-1 1])

%% Velocity quivers
[Xg, Yg] = meshgrid(X_axis, Y_axis); 
figure
colormap(brewermap(10,'GnBu'))
set(gca,'fontsize',16)
m_proj('Miller','lon',[-95 -82.5], 'lat',[22 31])
hold all
h=m_pcolor(X_axis, Y_axis, 0.5*(Umean.^2 + Vmean.^2)');
set(h,'edgecolor','none')
m_quiver(Xg+dx/2, Yg+dx/2, Umean',Vmean',2.5,'color','r','linewidth',2)
m_grid('linestyle','none','xaxisloc','top','yaxisloc','left')
m_coast('patch',[0 0 0]);
caxis([0 0.25])
h=colorbar;
ylabel(h,'MKE (m^2/s^2)')

%% Variance plots 

figure
colormap(brewermap(10,'GnBu'))
set(gca,'fontsize',16)
m_proj('Miller','lon',[-95 -82.5], 'lat',[22 31])
hold all
h= m_pcolor(X_axis, Y_axis, 0.5*(Ueke' + Veke'));
set(h,'edgecolor','none')
m_grid('linestyle','none','xaxisloc','top','yaxisloc','left')
m_coast('patch',[0 0 0]);
caxis([0 0.25])

for i = 1:size(Ell,1)
    for j = 1:size(Ell,2)
        m_plot(squeeze(Ell(i,j,1,:))+dx/2,squeeze(Ell(i,j,2,:))+dx/2,'r','linewidth',2)
    end
end
h= colorbar; 
ylabel(h, 'EKE (m^2/s^2)')
%% Number of samples 

nids(nids == 0) = NaN;
figure
set(gca,'fontsize',16)
colormap('cool')
m_proj('Miller','lon',[-95 -82.5], 'lat',[22 31])
hold all
h= m_pcolor(X_axis, Y_axis, nids');
set(h,'edgecolor','none')
caxis([0 5*10^4])
m_grid('linestyle','none','xaxisloc','top','yaxisloc','left')
m_coast('patch',[0 0 0]);
h=colorbar;
ylabel(h,'Number of samples','fontsize',16)
% caxis([-1 1])


