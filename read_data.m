clear all

fileid = fopen('GLAD_15min_filtered.dat');


C = textscan(fileid, '%s %s %s %f %f %f %f %f %f');


%%
float_num_last =1;
iter = 1;

float_num = zeros(length(C{1}),1);
float_datenum = zeros(length(C{1}),1);

for i =1:length(C{1})
    float_name = char(C{1}{i});
    
    float_num(i) = str2num(float_name(end-2:end));
    
    date = char(C{2}{i});
    time = char(C{3}{i});
     
    float_day = str2num(date(end-1:end));
    float_month = str2num(date(6:7));
    float_year = str2num(date(1:4));
    
    float_hour = str2num(time(1:2));
    float_min = str2num(time(4:5));
    float_sec = 0;
    
    float_datenum(i) = datenum(float_year, float_month, float_day, float_hour, float_min, float_sec);
%     
%     if float_num ~=float_num_last
%         iter = 1;
%     end
%     float_data(float_num).lat(iter) = C{4}(i);
%     float_data(float_num).lon(iter) = C{5}(i);
%     float_data(float_num).U(iter) = C{7}(i);
%     float_data(float_num).V(iter) = C{8}(i);
%     
%     float_num_last = float_num;
%     iter=iter+1;
    
end

%%
min_time = min(float_datenum); 
max_time = max(float_datenum); 
numfloats = max(float_num); 
dt = (float_datenum(2)- float_datenum(1));
ntimesteps = floor((max_time-min_time)/dt+1); 

%%

X = nan(ntimesteps,numfloats);
Y = nan(ntimesteps,numfloats);
U = nan(ntimesteps,numfloats);
V = nan(ntimesteps,numfloats);


for i =1:length(C{1})
    
    ii = floor((float_datenum(i)-min_time)/dt+1);
    jj = float_num(i); 
    
    X(ii,jj) = C{5}(i);
    Y(ii,jj) = C{4}(i);
    U(ii,jj) = C{7}(i);
    V(ii,jj) = C{8}(i);
%     
%     if float_num ~=float_num_last
%         iter = 1;
%     end
%     float_data(float_num).lat(iter) = C{4}(i);
%     float_data(float_num).lon(iter) = C{5}(i);
%     float_data(float_num).U(iter) = C{7}(i);
%     float_data(float_num).V(iter) = C{8}(i);
%     
%     float_num_last = float_num;
%     iter=iter+1;
    
end

%%
save glad_traj.mat X Y U V