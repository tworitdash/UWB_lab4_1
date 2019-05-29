%% UWB lab 4 - indoor tracking with 4 PulsOn 410 units
clear all;
close all;
%Positions of the 4 units in meters
P101 = [-1.6,1.6];           
P102 = [1.5,1.6];        
P103 = [1.9,-1.7];
P104 = [-1.2,-2.3];

% % Positions of the 4 units in meters
% P101 = [-1.35,2];           
% P102 = [-1.32,-1.5];        
% P103 = [1.73,2];        
% P104 = [1.90,-1.54];    


figure;
hold on; grid
plot(P101(1),P101(2),'ok','markerfacecolor','k');text(P101(1)+.1,P101(2),'101')
plot(P102(1),P102(2),'ok','markerfacecolor','k');text(P102(1)+.1,P102(2),'102')
plot(P103(1),P103(2),'ok','markerfacecolor','k');text(P103(1)+.1,P103(2),'103')
plot(P104(1),P104(2),'ok','markerfacecolor','k');text(P104(1)+.1,P104(2),'104')
    
plot(-1,1,'xr');text(-1.2,1.2,'A')
plot(0,1,'xr');text(.2,1.2,'B')
plot(1,1,'xr');text(1.2,1.2,'C')
plot(-1,0,'xr');text(-1.2,0.2,'D')
plot(0,0,'xr');text(.2,0.2,'E')
plot(1,0,'xr');text(1.2,0.2,'F')
plot(-1,-1,'xr');text(-1.2,-1.2,'G')
plot(0,-1,'xr');text(.2,-1.2,'H')
plot(1,-1,'xr');text(1.2,-1.2,'I')

axis image
title('GEOMETRY');xlabel('X (m)');ylabel('Y (m)')


%% ------------------------ Time zero calibration -------------------------
% !!!!!!!!!! NOT ESSENTIAL !!!!!!!!!!!!!
% Sphere placed at (0,0)%
% Loading data for t0 cal
%cd('Q:\UWB LAB Tracking\P');
% Load data from unit 101
[config,control,scans] = readMrmRetLog('sphere101006.csv');
Nscans101 = length(scans);
tstmp101 = [scans.T];
data101 = [scans.scn];
data101 = reshape(data101,[],Nscans101);
% Load data from unit 102
[config,control,scans] = readMrmRetLog('sphere102006.csv');
Nscans102 = length(scans);
tstmp102 = [scans.T];
data102 = [scans.scn];
data102 = reshape(data102,[],Nscans102);
% Load data from unit 103
[config,control,scans] = readMrmRetLog('sphere103006.csv');
Nscans103 = length(scans);
tstmp103 = [scans.T];
data103 = [scans.scn];
data103 = reshape(data103,[],Nscans103);
% Load data from unit 103
[config,control,scans] = readMrmRetLog('sphere102006.csv');
Nscans104 = length(scans);
tstmp104 = [scans.T];
data104 = [scans.scn];
data104 = reshape(data104,[],Nscans104);

% Load Sphere backgrounds

[config,control,scans] = readMrmRetLog('BGsphere101007.csv');
BGNscans101 = length(scans);
BGtstmp101 = [scans.T];
BGdata101 = [scans.scn];
BGdata101 = reshape(BGdata101,[],BGNscans101);
% Load data from unit 102
[config,control,scans] = readMrmRetLog('BGsphere102007.csv');
BGNscans102 = length(scans);
BGtstmp102 = [scans.T];
BGdata102 = [scans.scn];
BGdata102 = reshape(BGdata102,[],BGNscans102);
% Load data from unit 103
[config,control,scans] = readMrmRetLog('BGsphere103007.csv');
BGNscans103 = length(scans);
BGtstmp103 = [scans.T];
BGdata103 = [scans.scn];
BGdata103 = reshape(BGdata103,[],BGNscans103);
% Load data from unit 103
[config,control,scans] = readMrmRetLog('BGsphere104007.csv');
BGNscans104 = length(scans);
BGtstmp104 = [scans.T];
BGdata104 = [scans.scn];
BGdata104 = reshape(BGdata104,[],BGNscans104);

% fast time / range uncalibrated
t = linspace(scans(1).Tstrt,scans(1).Tstp,scans(1).Nscn)/1000;   % ns
rng = 3e8*(t-t(1))/2e9; % m

% 101
figure;imagesc(abs(data103));title('101')
%BG103 = repmat(mean(data103(:,50:110),2),1,Nscans103);
BG101 = repmat(mean(BGdata101,2),1,Nscans101);
data101bg = data101 - BG101;
figure;imagesc(abs(data101bg(200:400,:)));title('101')

rngSph101 = rng(230); % range of sphere in unit 101 (m)
rngSph101Real = norm(P101); % real range from unit 101 (m)
rngOffset101 = rngSph101Real - rngSph101; % range offset for unit 103

% 102
figure;imagesc(abs(data102));title('102')
%BG102 = repmat(mean(data102(:,50:110),2),1,Nscans102);
BG102 = repmat(mean(BGdata102,2),1,Nscans102);
data102bg = data102 - BG102;
figure;imagesc(abs(data102bg));title('102')

rngSph102 = rng(317); % range of sphere in unit 102 (m)
rngSph102Real = norm(P102); % real range from unit 102 (m)
rngOffset102 = rngSph102Real - rngSph102; % range offset for unit 102

% 103
figure;imagesc(abs(data103));title('103')
%BG103 = repmat(mean(data103(:,50:110),2),1,Nscans103);
BG103 = repmat(mean(BGdata103,2),1,Nscans103);
data103bg = data103 - BG103;
figure;imagesc(abs(data103bg));title('103')

rngSph103 = rng(395); % range of sphere in unit 103 (m)
rngSph103Real = norm(P103); % real range from unit 103 (m)
rngOffset103 = rngSph103Real - rngSph103; % range offset for unit 103

% 104
figure;imagesc(abs(data104));title('104')
%BG103 = repmat(mean(data103(:,50:110),2),1,Nscans103);
BG104 = repmat(mean(BGdata104,2),1,Nscans104);
data104bg = data104 - BG104;
figure;imagesc(abs(data104bg));title('104')

rngSph104 = rng(345); % range of sphere in unit 103 (m)
rngSph104Real = norm(P104); % real range from unit 103 (m)
rngOffset104 = rngSph104Real - rngSph104; % range offset for unit 103


%% Loading Data for tracking
%cd('D:\EDUCATION\UWB - Q4\Tracking lab\Data');
% Load data from unit 101
[config,control,scans] = readMrmRetLog('A101008.csv');
Nscans101 = length(scans);
tstmp101 = [scans.T];
data101 = [scans.scn];
data101 = reshape(data101,[],Nscans101);
% Load data from unit 102
[config,control,scans] = readMrmRetLog('A102008.csv');
Nscans102 = length(scans);
tstmp102 = [scans.T];
data102 = [scans.scn];
data102 = reshape(data102,[],Nscans102);
% Load data from unit 103
[config,control,scans] = readMrmRetLog('A103008.csv');
Nscans103 = length(scans);
tstmp103 = [scans.T];
data103 = [scans.scn];
data103 = reshape(data103,[],Nscans103);
% Load data from unit 104
[config,control,scans] = readMrmRetLog('A103008.csv');
Nscans104 = length(scans);
tstmp104 = [scans.T];
data104 = [scans.scn];
data104 = reshape(data104,[],Nscans104);

%% Data alignment
TStampStart = max([tstmp101(1),tstmp102(1),tstmp103(1),tstmp104(1)]);
TStampStop = min([tstmp101(end),tstmp102(end),tstmp103(end),tstmp104(end)]);

t101start_i = find(tstmp101>=TStampStart,1,'first');
t101stop_i = find(tstmp101<=TStampStop,1,'last');

t102start_i = find(tstmp102>=TStampStart,1,'first');
t102stop_i = find(tstmp102<=TStampStop,1,'last');

t103start_i = find(tstmp103>=TStampStart,1,'first');
t103stop_i = find(tstmp103<=TStampStop,1,'last');

t104start_i = find(tstmp104>=TStampStart,1,'first');
t104stop_i = find(tstmp104<=TStampStop,1,'last');

data101_algnd = data101(:,t101start_i:t101stop_i);
data102_algnd = data102(:,t102start_i:t102stop_i);
data103_algnd = data103(:,t103start_i:t103stop_i);
data104_algnd = data104(:,t104start_i:t104stop_i);

[Nsamp,Nscans101] = size(data101_algnd);
[Nsamp,Nscans102] = size(data102_algnd);
[Nsamp,Nscans103] = size(data103_algnd);
[Nsamp,Nscans104] = size(data104_algnd);

figure;imagesc(abs(data101_algnd));title('101')
figure;imagesc(abs(data102_algnd));title('102')
figure;imagesc(abs(data103_algnd));title('103')
figure;imagesc(abs(data104_algnd));title('104')

%% Background removal
BGidx = 110:130;  % select the slow time indexes containing the background data
BG101 = mean(data101_algnd(:,BGidx),2);   
BG102 = mean(data102_algnd(:,BGidx),2);
BG103 = mean(data103_algnd(:,BGidx),2);
BG104 = mean(data104_algnd(:,BGidx),2);

data101_algnd_bg = data101_algnd - repmat(BG101,1,Nscans101); % remove background
data102_algnd_bg = data102_algnd - repmat(BG102,1,Nscans102);
data103_algnd_bg = data103_algnd - repmat(BG103,1,Nscans103);
data104_algnd_bg = data104_algnd - repmat(BG104,1,Nscans104);

figure;imagesc(abs(data101_algnd_bg));title('101')
figure;imagesc(abs(data102_algnd_bg));title('102')
figure;imagesc(abs(data103_algnd_bg));title('103')
figure;imagesc(abs(data104_algnd_bg));title('104')

Motionidx = 25:95;  % select the slow time indexes containing the trajectory data
path101 = data101_algnd_bg(:,Motionidx); 
path102 = data102_algnd_bg(:,Motionidx);
path103 = data103_algnd_bg(:,Motionidx);
path104 = data104_algnd_bg(:,Motionidx);

% normalisation
path101 = path101/max(abs(path101(:))); 
path102 = path102/max(abs(path102(:))); 
path103 = path103/max(abs(path103(:)));
path104 = path104/max(abs(path104(:)));

figure;imagesc(abs(path101));title('101')
figure;imagesc(abs(path102));title('102')
figure;imagesc(abs(path103));title('103')
figure;imagesc(abs(path104));title('104');
%% Motion filter FIR4 if backgound subtraction unsuccessful (needed for some trajectories)
Motionidx = 20:75;
w = [1 -.6 -.3 -.1];

path101 = filter(w,1,data101_algnd,[],2);
path102 = filter(w,1,data102_algnd,[],2);
path103 = filter(w,1,data103_algnd,[],2);
path104 = filter(w,1,data104_algnd,[],2);

% normalisation
path101 = path101(:,Motionidx)/max(max(abs(path101(:,Motionidx)))); 
path102 = path102(:,Motionidx)/max(max(abs(path102(:,Motionidx)))); 
path103 = path103(:,Motionidx)/max(max(abs(path103(:,Motionidx))));
path104 = path104(:,Motionidx)/max(max(abs(path104(:,Motionidx))));

figure;imagesc(abs(path101));title('101')
figure;imagesc(abs(path102));title('102')
figure;imagesc(abs(path103));title('103')
figure;imagesc(abs(path104));title('104');

%% Range estimation
% fast time / range uncalibrated
t = linspace(scans(1).Tstrt,scans(1).Tstp,scans(1).Nscn)/1000;   % ns
rng = 3e8*(t-t(1))/2e9; % m
[Nrng,Nscans] = size(path101);
t_slow = 0:.5:.5*(Nscans-1);

% 101
close all;
estimated_range=[];
for i1=1:Nscans % Long time
    i2 = 80; % Avoiding picking antenna coupling
    while abs(path101(i2,i1)) < 0.04 && i2 < Nrng-1    % threshold 0.06
        i2 = i2 + 1;
    end
    estimated_range(i1) = rng(i2);
end
figure;plot(estimated_range)
median_range101 = medfilt1(estimated_range,3)+rngOffset101;
figure;imagesc(t_slow,rng,abs(path101));title('101')
hold
plot(t_slow,median_range101,'w')  % readjust the threshold and coupling index based on this plot

% 102
estimated_range=[];
for i1=1:Nscans % Long time
    i2 = 1; % Avoiding picking antenna coupling
    while abs(path102(i2,i1)) < 0.09 && i2 < Nrng-2
        i2 = i2 + 1;
    end
    estimated_range(i1) = rng(i2);
end
figure;plot(estimated_range)
median_range102 = medfilt1(estimated_range,3);%+rngOffset102;
figure;imagesc(t_slow,rng,abs(path102));title('102')
hold
plot(t_slow,median_range102,'w')

% 103
estimated_range=[];
for i1=1:Nscans % Long time
    i2 = 100; % Avoiding picking antenna coupling
    while abs(path103(i2,i1)) < 0.07 && i2 < Nrng-1
        i2 = i2 + 1;
    end
    estimated_range(i1) = rng(i2);
end
figure;plot(estimated_range)
median_range103 = medfilt1(estimated_range,3);%+rngOffset103;
figure;imagesc(t_slow,rng,abs(path103));title('103')
hold
plot(t_slow,median_range103,'w')

% 104
estimated_range=[];
for i1=1:Nscans % Long time
    i2 = 100; % Avoiding picking antenna coupling
    while abs(path104(i2,i1)) < 0.04 && i2 < Nrng-1
        i2 = i2 + 1;
    end
    estimated_range(i1) = rng(i2);
end
figure;plot(estimated_range)
median_range104 = medfilt1(estimated_range,3);%+rngOffset104;
figure;imagesc(t_slow,rng,abs(path104));title('104')
hold
plot(t_slow,median_range104,'w')

%% Target localisation

for scani = 1:Nscans
    [a, b] = circcirc(P101(1),P101(2),median_range101(scani),P102(1),P102(2),median_range102(scani));
    x12(scani,:) = a;
    y12(scani,:) = b;
end
figure;plot(x12(:,1),y12(:,2));title('101 - 102');axis([-1.5 1.5 -1.5 1.5])   

for scani = 1:Nscans
    [a, b] = circcirc(P101(1),P101(2),median_range101(scani),P103(1),P103(2),median_range103(scani));
    x13(scani,:) = a;
    y13(scani,:) = b;
end
figure;plot(x13(:,1),y13(:,2));title('101 - 103');axis([-1.5 1.5 -1.5 1.5])    

for scani = 1:Nscans
    [a, b] = circcirc(P101(1),P101(2),median_range101(scani),P104(1),P104(2),median_range104(scani));
    x14(scani,:) = a;
    y14(scani,:) = b;
end
figure;plot(x14(:,1),y14(:,2));title('101 - 104');axis([-1.5 1.5 -1.5 1.5])  

for scani = 1:Nscans
    [a, b] = circcirc(P102(1),P102(2),median_range102(scani),P103(1),P103(2),median_range103(scani));
    x23(scani,:) = a;
    y23(scani,:) = b;
end
figure;plot(x23(:,2),y23(:,2));title('102 - 103');axis([-1.5 1.5 -1.5 1.5])  

for scani = 1:Nscans
    [a, b] = circcirc(P102(1),P102(2),median_range102(scani),P104(1),P104(2),median_range104(scani));
    x24(scani,:) = a;
    y24(scani,:) = b;
end
figure;plot(x24(:,1),y24(:,1));title('102 - 104');axis([-1.5 1.5 -1.5 1.5])  

for scani = 1:Nscans
    [a, b] = circcirc(P103(1),P103(2),median_range103(scani),P104(1),P104(2),median_range104(scani));
    x34(scani,:) = a;
    y34(scani,:) = b;
end
figure;plot(x34(:,2),y34(:,2));title('103 - 104');axis([-1.5 1.5 -1.5 1.5])  


%% FINAL PLOT
% x = nanmean([x14(:,1),x24(:,2),x34(:,2),x12(:,1)],2);
% y = nanmean([y14(:,2),y24(:,1),y34(:,2),y12(:,2)],2);
% 
x = nanmean([x24(:,1)],2);
y = nanmean([y24(:,1)],2);


figure;plot(smooth(x,3),smooth(y,3),'linewidth',2);axis([-2.5 2.5 -3 3]);grid
hold on
plot(P101(1),P101(2),'ok','markerfacecolor','k');text(P101(1)+.1,P101(2),'101')
plot(P102(1),P102(2),'ok','markerfacecolor','k');text(P102(1)+.1,P102(2),'102')
plot(P103(1),P103(2),'ok','markerfacecolor','k');text(P103(1)+.1,P103(2),'103')
plot(P104(1),P104(2),'ok','markerfacecolor','k');text(P104(1)+.1,P104(2),'104')
    
plot(-1,1,'xr');text(-1.2,1.2,'A')
plot(0,1,'xr');text(.2,1.2,'B')
plot(1,1,'xr');text(1.2,1.2,'C')
plot(-1,0,'xr');text(-1.2,0.2,'D')
plot(0,0,'xr');text(.2,0.2,'E')
plot(1,0,'xr');text(1.2,0.2,'F')
plot(-1,-1,'xr');text(-1.2,-1.2,'G')
plot(0,-1,'xr');text(.2,-1.2,'H')
plot(1,-1,'xr');text(1.2,-1.2,'I')

axis image
title('AEI trajectory');xlabel('X (m)');ylabel('Y (m)')

%% Least Square method:
Theta_x = zeros(1, Nscans);
Theta_y = zeros(1, Nscans);

    x1 = P101(1);
    x2 = P102(1);
    x3 = P103(1);
    x4 = P104(1);
    
    y1 = P101(2);
    y2 = P102(2);
    y3 = P103(2);
    y4 = P104(2);

for i = 1:Nscans
    
    r1 = median_range101(i);
    r2 = median_range102(i);
    r3 = median_range103(i);
    r4 = median_range104(i);
    
    
    
    X = [r1^2 - r4^2 - x1^2 + x4^2 - y1^2 + y4^2; r2^2 - r4^2 - x2^2 + x4^2 - y2^2 + y4^2; ... 
        r3^2 - r4^2 - x3^2 + x4^2 - y3^2 + y4^2];
    
    H = 2.* [(x4 - x1) (y4 - y1);  (x4 - x2) (y4 - y2); (x4 - x3) (y4 - y3)];
    
    Theta = inv(H'*H) * H' * X;
    Theta_x(i) = Theta(1);
    Theta_y(i) = Theta(2);
    
    
    

    
    
    
end

    figure;plot(smooth(Theta_x,3),smooth(Theta_y,3),'linewidth',2);axis([-2.5 2.5 -3 3]);grid
    hold on
    plot(P101(1),P101(2),'ok','markerfacecolor','k');text(P101(1)+.1,P101(2),'101')
    plot(P102(1),P102(2),'ok','markerfacecolor','k');text(P102(1)+.1,P102(2),'102')
    plot(P103(1),P103(2),'ok','markerfacecolor','k');text(P103(1)+.1,P103(2),'103')
    plot(P104(1),P104(2),'ok','markerfacecolor','k');text(P104(1)+.1,P104(2),'104')

    plot(-1,1,'xr');text(-1.2,1.2,'A')
    plot(0,1,'xr');text(.2,1.2,'B')
    plot(1,1,'xr');text(1.2,1.2,'C')
    plot(-1,0,'xr');text(-1.2,0.2,'D')
    plot(0,0,'xr');text(.2,0.2,'E')
    plot(1,0,'xr');text(1.2,0.2,'F')
    plot(-1,-1,'xr');text(-1.2,-1.2,'G')
    plot(0,-1,'xr');text(.2,-1.2,'H')
    plot(1,-1,'xr');text(1.2,-1.2,'I')

    axis image
    title('AEI trajectory');xlabel('X (m)');ylabel('Y (m)')

    