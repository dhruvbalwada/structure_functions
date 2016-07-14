clear all 
close all 

load('glad_traj.mat')

topo.z = ncread('../bathymetry/GEBCO_GRIDONE_2D.nc', 'elevation');
topo.lon = ncread('../bathymetry/GEBCO_GRIDONE_2D.nc', 'lon');
topo.lat = ncread('../bathymetry/GEBCO_GRIDONE_2D.nc', 'lat');

%% 

idx =find(topo.lon > -99 & topo.lon<-80);
idy =find(topo.lat > 19 & topo.lat<32);
%%

figure 
set(gca,'fontsize',16,'fontname','times')
colormap(gray)
m_proj('Miller','lon',[-98 -80.5], 'lat',[19 31])
hold all
m_contourf(topo.lon(idx), topo.lat(idy), topo.z(idx,idy)', [-5000:1000:0],'edgecolor','none')
m_contourf(topo.lon(idx), topo.lat(idy), topo.z(idx,idy)', [0 0],'color','k')

m_plot(X,Y,'linewidth',1.5)

for i =1:size(X,2)
    ids = find(~isnan(X(:,i)),1);
    h=m_plot(X(ids,i), Y(ids,i),'.','color','k','markersize',10);
%     set(h,'markerfacecolor','k')    
end
% m_etopo2('patch',[0 0 0])
% m_tbase('contourf')
% m_coast('patch',[0 0 0])
m_grid('linestyle','none','xaxisloc','top','yaxisloc','right')
