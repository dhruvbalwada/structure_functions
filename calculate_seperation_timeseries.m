function [sep] = calculate_seperation_timeseries(traj)

nflts = size(traj.X,2);
ndays = size(traj.X,1);

n=1;
npairs = factorial(nflts)/factorial(nflts-2)/factorial(2);
% sep(npairs).X = nan*ones(ndays,1);
% sep(npairs).Y = nan*ones(ndays,1);
% sep(npairs).dist = nan*ones(ndays,1);
% sep(npairs).P1 = nan*ones(ndays,1);
% sep(npairs).P2 = nan*ones(ndays,1);
% sep(npairs).dP = nan*ones(ndays,1);
% sep(npairs).T1 = nan*ones(ndays,1);
% sep(npairs).T2 = nan*ones(ndays,1);
% sep(npairs).names = [];

for i=1:nflts-1
    for j = i+1:nflts
        
        sep(n).X = abs(traj.X(:,i) - traj.X(:,j)).*cosd(0.5*(traj.Y(:,i)+traj.Y(:,j)))*111321;
        sep(n).Y = abs(traj.Y(:,i) - traj.Y(:,j))*111321;
        sep(n).dist = sqrt(sep(n).X.^2 + sep(n).Y.^2);
        sep(n).X1 = traj.X(:,i);
        sep(n).X2 = traj.X(:,j);
        sep(n).Y1 = traj.Y(:,i);
        sep(n).Y2 = traj.Y(:,j);
        sep(n).U1 = traj.U(:,i);
        sep(n).U2 = traj.U(:,j);
        sep(n).V1 = traj.V(:,i);
        sep(n).V2 = traj.V(:,j);
%         sep(n).P1 = traj.Pi(:,i);
%         sep(n).P2 = traj.Pi(:,j);
%         sep(n).dP = abs(sep(n).P1 - sep(n).P2);
%         sep(n).T1 = traj.Ti(:,i);
%         sep(n).T2 = traj.Ti(:,j);
%         sep(n).names(1) = str2num(traj.name(:,i)');
%         sep(n).names(2) = str2num(traj.name(:,j)');
        n = n+1;
    end
end
