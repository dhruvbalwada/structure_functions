% new code to try and do the helmholtz decomposition on the observational
% structure functions.

clear all
close all

traj = load ('glad_traj.mat');

%% filtering 
f = 2*2*pi/24/3600*sind(27.89); 
dt = 15*60; 
filt_flag=0;

for i =1:size(traj.X,2)
    if length(find(~isnan(traj.U(:,i))))>0
         ids(i) = find(~isnan(traj.U(:,i)),1) ;
        ide(i) = find(~isnan(traj.U(:,i)),1, 'last') ;
        tot_nonans(i) = length(find(~isnan(traj.U(:,i))));
    else
        ids(i) = NaN;
        ide(i) = NaN;
        tot_nonans(i) = 0;
    end
end

if filt_flag==1
[b, a] = butter(2, [f-0.20*f f+0.22*f]*(dt/pi),'stop');

U_filt = NaN*traj.U;
V_filt = NaN*traj.V;

for i =1:size(traj.X,2)
    if tot_nonans(i)>0
        U_filt(ids(i):ide(i),i) = filtfilt(b,a, traj.U(ids(i):ide(i),i));
        V_filt(ids(i):ide(i),i) = filtfilt(b,a, traj.V(ids(i):ide(i),i));
    end
end
disp('filt things')
traj.U=U_filt;
traj.V=V_filt;
end
%%
for i =1:size(traj.X,2) 
    id = find(~isnan(traj.X(:,i)));
    ndays(i) = length(id);
end
id =find(ndays~=0); 

id = id(end/2:end);

%%
traj.X = traj.X(:,id);
traj.Y = traj.Y(:,id);
traj.U = traj.U(:,id);
traj.V = traj.V(:,id);


int_over = 4; % what length of gaps should be filled in using an interpolant
% calculate time series of pair separation
sep = calculate_seperation_timeseries(traj);

%%
% plevel = [1200 1600];

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

clear s2ll s2tt
s2ll = zeros(length(dist_axis),1);
s2tt = zeros(length(dist_axis),1);
npairs = zeros(length(dist_axis),1);
%
% loop for different distance classes
for i =1:length(dist_axis)
    dull = []; dutt = [];
    % loop over different pairs
    for j = 1:length(sep)
        dull_temp = []; dutt_temp = [];
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
            
            % components of velocity differences
            dux(k) = (sep(j).U1(id(k))-sep(j).U2(id(k)));
            duy(k) = (sep(j).V1(id(k))-sep(j).V2(id(k)));
            
            % convert to longitudnal and
            dull_temp(k) = dux(k)*rx(k) + duy(k)*ry(k);
            dutt_temp(k) = duy(k)*rx(k) - dux(k)*ry(k);
%             dull_temp = dux(k)*rx(k) + duy(k)*ry(k);
%             dutt_temp = duy(k)*rx(k) - dux(k)*ry(k);
            
%             if ~isnan(dull_temp) & ~isnan(dutt_temp)
%             s2ll(i) = s2ll(i)+dull_temp;
%             s2tt(i) = s2tt(i)+dutt_temp;
%             npairs(i) = npairs(i)+1; 
%             end
        end
        if ~isempty(dull_temp) 
            dull = [dull; dull_temp'];
            dutt = [dutt; dutt_temp'];
        end
    end
    
    struct_pairs(i).dull = dull;
    struct_pairs(i).dutt = dutt;
    
    if ~isempty(dull) & length(dull)>20
         s2ll(i) = nanmean(dull.^2);
         s2tt(i) = nanmean(dutt.^2);
%         s3(i) = nanmean(vel_pairs.^3);
%         s2err(i) = nanstd(vel_pairs.^2)/sqrt(length(vel_pairs));
%         s3err(i) = nanstd(vel_pairs.^3)/sqrt(length(vel_pairs));
    else
        s2ll(i) = NaN;
        s2tt(i) = NaN;
%         s3(i) = NaN;
%         s2err(i) = NaN;
%         s3err(i) = NaN;
    end
    disp(i)
end

%%
for i =1:length(struct_pairs)
    s3lll(i) = nanmean(struct_pairs(i).dull.^3);
    s3ltt(i) = nanmean(struct_pairs(i).dutt.^2.*struct_pairs(i).dull);
end

%%
for i =1:length(struct_pairs)
    s4l(i) = nanmean(struct_pairs(i).dull.^4)/(s2ll(i).^2);
    s4t(i) = nanmean(struct_pairs(i).dutt.^4)/(s2tt(i).^2);
end

%%
figure 

loglog(dist_axis/1000, s4l,'linewidth',2), hold all
loglog(dist_axis/1000, s4t,'linewidth',2)
axis([10^-2 10^3 1 100])
ylabel('Flatness')
xlabel('r(km)')
legend('F_l', 'F_t')
grid
%%
s3 = s3lll+s3ltt;

figure
loglog(dist_axis/1000, abs(s3./dist_axis),'+-','linewidth', 2), hold all
loglog(dist_axis/1000, -(s3./dist_axis),'o','linewidth', 2)
loglog(dist_axis/1000, abs(s3lll./dist_axis),'+-','linewidth', 2), hold all
loglog(dist_axis/1000, -(s3lll./dist_axis),'o','linewidth', 2)
loglog(dist_axis/1000, abs(s3ltt./dist_axis),'+-','linewidth', 2), hold all
loglog(dist_axis/1000, -(s3ltt./dist_axis),'o','linewidth', 2)
legend('|S3_{lll} + S3{ltt}|','-(S3_{lll} + S3{ltt})', '|S3_{lll}|', '-S3_{lll}','|S3{ltt}|','-(S3{ltt})')
grid
axis([10^-2 10^3 10^-9 10^-6])
xlabel('r (km)')
ylabel('(S3)/r')

%%
figure
loglog(dist_axis/1000, abs(s3),'+-','linewidth', 2), hold all
loglog(dist_axis/1000, -(s3),'o','linewidth', 2)
% loglog(dist_axis/1000, abs(s3lll),'+-','linewidth', 2), hold all
% loglog(dist_axis/1000, -(s3lll),'o','linewidth', 2)
% loglog(dist_axis/1000, abs(s3ltt),'+-','linewidth', 2), hold all
% loglog(dist_axis/1000, -(s3ltt),'o','linewidth', 2)
% legend('|S3_{lll} + S3{ltt}|','-(S3_{lll} + S3{ltt})', '|S3_{lll}|', '-S3_{lll}','|S3{ltt}|','-(S3{ltt})')

ao = 10^-8;
a1 = 10^-4;
xax = [0.001:10:10^3];
yax2 = ao*(xax.^3);
yax1 = a1*(xax.^1);
% yax43 = 2*a1*(xax.^(2/3));
loglog(xax, yax2,'--','color',[0.5 0.5 0.5],'linewidth',2)
loglog(xax, yax1,'-.','color',[0.5 0.5 0.5],'linewidth',2)
grid

axis([10^-2 10^3 10^-8 1])
xlabel('r (km)')
ylabel('(S3)')

%%
for i =1:length(struct_pairs)
    npairs(i) = length(struct_pairs(i).dull);
end
%%
figure,
loglog(dist_axis/1000,npairs,'s-','linewidth',2)
axis([10^-2 10^3 1000 10^8]),grid
xlabel('r (km)')
ylabel('Number of pairs')
%% do the integrals to calculate the decomposition to rotational and divergent part 
clear s2rr s2dd
mid_dist_axis = 0.5*(dist_axis(1:end-1)+dist_axis(2:end));
mid_diff_du = 0.5*((s2tt(1:end-1)-s2ll(1:end-1))+(s2tt(2:end)-s2ll(2:end))); 

for i =2:length(dist_axis)
    s2rr(i) = s2tt(i) + nansum(1./mid_dist_axis(1:i-1).*mid_diff_du(1:i-1)'.*diff(dist_axis(1:i))); 
    s2dd(i) = s2ll(i) - nansum(1./mid_dist_axis(1:i-1).*mid_diff_du(1:i-1)'.*diff(dist_axis(1:i))); 
end

%% plot the components 
% close all
figure 
loglog(dist_axis/1000, s2ll,'linewidth',2)
 hold all 
 loglog(dist_axis/1000, s2tt,'linewidth',2)
 hold all
loglog(dist_axis/1000, s2tt+s2ll, '-','linewidth',2)
hold all 
loglog(dist_axis/1000, s2rr+s2dd, '+','linewidth',2)
loglog(dist_axis/1000, s2rr, '-','linewidth',2)
loglog(dist_axis/1000, s2dd, '-','linewidth',2)

ao = 10^-5;
a1 = 10^-4;
xax = [0.001:10:10^3];
yax2 = ao*(xax.^2);
yax1 = a1*(xax.^1);
yax43 = 2*a1*(xax.^(2/3));
loglog(xax, yax2,'--','color',[0.5 0.5 0.5],'linewidth',2)
loglog(xax, yax1,'-.','color',[0.5 0.5 0.5],'linewidth',2)
loglog(xax, yax43,'-','color',[0.5 0.5 0.5],'linewidth',2)
 axis([10^-2 10^3 5*10^-5 6*10^-1])
 xlabel('r (km)')
 ylabel('<\delta u^2>')
 legend('S2_{ll}','S2_{tt}','S2{tt} + S2_{ll}','S2_{rr}+S2_{dd}', 'S2_{rr}','S2_{dd}') 
grid

%% 

figure 
semilogx(dist_axis/1000, s2tt./s2ll,'o-','linewidth',2)
xlabel('r(km)')
ylabel('S_{tt}/S_{ll}')
grid
axis([10^-2 10^3 0.5 2.5]) 