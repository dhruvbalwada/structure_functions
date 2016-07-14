% new code to try and do the helmholtz decomposition on the observational
% structure functions.

clear all
close all

traj = load ('glad_traj.mat');
for i =1:size(traj.X,2) 
    id = find(~isnan(traj.X(:,i)));
    ndays(i) = length(id);
end
id =find(ndays~=0); 
traj.X = traj.X(:,id);
traj.Y = traj.Y(:,id);
traj.U = traj.U(:,id);
traj.V = traj.V(:,id);

% remove the means
meanU = nanmean(traj.U(:));
meanV = nanmean(traj.V(:));

traj.U = traj.U - meanU; 
traj.V = traj.V - meanV; 

int_over = 4; % what length of gaps should be filled in using an interpolant
pressure_diff = 100; % what depth must a float drop to remove that section of the trajectory
% calculate time series of pair separation
sep = calculate_seperation_timeseries(traj);

%%

gamma = 1.4;

dist_bin(1) = 0.1; % in m
dist_bin = gamma.^[0:100]*dist_bin(1);
id = find(dist_bin>1000*10^3,1);
dist_bin = dist_bin(1:id);
dist_bin(2:end+1) = dist_bin(1:end);
dist_bin(1) = 0;
dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));
%     separationr = sqrt(dispersion(d).disp);
%%

clear Rll Rtt
Rll = zeros(length(dist_axis),1);
Rtt = zeros(length(dist_axis),1);
npairs = zeros(length(dist_axis),1);
%
% loop for different distance classes
for i =1:length(dist_axis)
    ull = []; utt = [];
    % loop over different pairs
    for j = 1:length(sep)
        ull_temp = []; utt_temp = [];
        % find id of the pairs in a particular geographical regime (lon
        % and pressure and within a certain distance from each other)
        
        id = find(sep(j).dist<dist_bin(i+1) & sep(j).dist>=dist_bin(i));
        
%         id = find(sep(j).dist<dist_bin(i+1) & sep(j).dist>=dist_bin(i) ...
%             & abs(sep(j).T1-sep(j).T2)<diff_temp & sep(j).X1<-70 & sep(j).X2<-70 & ...
%             sep(j).P1>plevel(1) & sep(j).P1<plevel(2));

%         ids2keep = find(id~=1117);
        %  id = id(ids2keep);
        
        % loop over the different pairs that lie in the range
        for k =1:length(id)
            % components
            rx(k) = (sep(j).X1(id(k)) - sep(j).X2(id(k)))*cosd(0.5*(sep(j).Y1(id(k))+sep(j).Y2(id(k))));
            ry(k) = (sep(j).Y1(id(k)) - sep(j).Y2(id(k)));
            magr(k) = sqrt(rx(k).^2+ry(k).^2);
            % normalize to unit vectors
            rx(k) = rx(k)/magr(k); ry(k) = ry(k)/magr(k);
            
           
            % convert to longitudnal and
            u1l = sep(j).U1(id(k))*rx(k) + sep(j).V1(id(k))*ry(k);
            
            u2l = sep(j).U2(id(k))*rx(k) + sep(j).V2(id(k))*ry(k);
            
            u1t = sep(j).V1(id(k))*rx(k) - sep(j).U1(id(k))*ry(k);
            u2t = sep(j).V2(id(k))*rx(k) - sep(j).U2(id(k))*ry(k);
            
            ull_temp(k) = u1l*u2l;
            
            utt_temp(k) = u1t*u2t;
            

           npairs(i) = npairs(i)+1; 

        end
        if ~isempty(ull_temp) 
            ull = [ull; ull_temp'];
            utt = [utt; utt_temp'];
        end
    end
    
    struct_pairs(i).ull = ull;
    struct_pairs(i).utt = utt;
    
    if ~isempty(ull) & length(ull)>20
         Rll(i) = nanmean(ull);
         Rtt(i) = nanmean(utt);

    else
        Rll(i) = NaN;
        Rtt(i) = NaN;

    end
    disp(i)
end


%% do the integrals to calculate the decomposition to rotational and divergent part 
clear Rrr Rdd
mid_dist_axis = 0.5*(dist_axis(1:end-1)+dist_axis(2:end));
mid_diff_du = 0.5*((Rtt(1:end-1)-Rll(1:end-1))+(Rtt(2:end)-Rll(2:end))); 

diff_dist = diff(dist_axis);

for i =1:length(dist_axis)-1
    Rrr(i) = Rtt(i) - nansum(1./mid_dist_axis(i:end).*mid_diff_du(i:end)'.*diff_dist(i:end)); 
    Rdd(i) = Rll(i) + nansum(1./mid_dist_axis(i:end).*mid_diff_du(i:end)'.*diff_dist(i:end)); 
end

Rrr(end+1) = 0;
Rdd(end+1) = 0;
%% plot the components 


% close all
figure 
loglog(dist_axis/1000, Rll,'-','linewidth',2)
 hold all 
 loglog(dist_axis/1000, Rtt,'-','linewidth',2)
 hold all
loglog(dist_axis/1000, Rtt+Rll, '-','linewidth',2)
hold all 
loglog(dist_axis/1000, Rrr+Rdd, '+','linewidth',2)
loglog(dist_axis/1000, Rrr, '-','linewidth',2)
loglog(dist_axis/1000, abs(Rdd), '-','linewidth',2)
id = find (Rdd<=0);
loglog(dist_axis(id)/1000, abs(Rdd(id)), 'o','linewidth',2)
% 
% ao = 10^-5;
% a1 = 10^-4;
% xax = [0.001:10:10^3];
% yax2 = ao*(xax.^2);
% yax1 = a1*(xax.^1);
% yax43 = 2*a1*(xax.^(2/3));
% loglog(xax, yax2,'--','color',[0.5 0.5 0.5],'linewidth',2)
% loglog(xax, yax1,'-.','color',[0.5 0.5 0.5],'linewidth',2)
% loglog(xax, yax43,'-','color',[0.5 0.5 0.5],'linewidth',2)
 axis([10^-2 10^3 10^-4 6*10^-1])
 xlabel('r (km)')
 ylabel('R (r) (m^2/s^2)')
legend('R_{ll}','R_{tt}','R_{tt} + R_{ll}','R_{rr}+R_{dd}', 'R_{rr}','R_{dd}') 
grid

%% convert to other kind

Rorr = Rrr(1); 
ido = find(Rrr == Rorr); 
Rodd = Rdd(ido);

Roll = max(Rll); 
Rott = max(Rtt); 

Dll = 2*Roll - 2*Rll; 
Dtt = 2*Rott - 2*Rtt; 

Drr = 2*Rorr - 2*Rrr;

Ddd = 2*Rodd - 2*Rdd;

figure 
loglog(dist_axis/1000, Drr), hold all
loglog(dist_axis/1000, Ddd)
figure,
loglog(dist_axis/1000, Dll,'linewidth',2), hold all
loglog(dist_axis/1000, Dtt,'linewidth',2)
axis([10^-2 10^3 10^-4 6*10^-1])
legend('S_{ll}','S_{tt}')
 xlabel('r (km)')
 ylabel('S(r) (m^2/s^2)')