%% UWB lab 4 - indoor tracking with 4 PulsOn 410 units
% Units positions in meters
P101 = [-1.5,-1.6];     
P102 = [-1.4,1.5];
P103 = [1.7,1.9];
P104 = [2.3,-1.3];

% Floor markers position in meters
A = [-1,1];
B = [0,1];
C = [1,1];
D = [-1,0];
E = [0,0];
F = [1,0];
G = [-1,-1];
H = [0,-1];
I = [1,-1];

MarkPos = [A;B;C;D;E;F;G;H;I];


%% ------------------------ Time zero calibration -------------------------
% Sphere placed at (0,0)%% Loading data for t0 cal
cd('D:\UWB LABS\Tracking gr2');
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
% Load data from unit 104
[config,control,scans] = readMrmRetLog('sphere104006.csv');
Nscans104 = length(scans);
tstmp104 = [scans.T];
data104 = [scans.scn];
data104 = reshape(data104,[],Nscans104);

% SPHERE BCKGRND
% Load data from unit 101
[config,control,scans] = readMrmRetLog('BGsphere101007.csv');
bgNscans101 = length(scans);
bgtstmp101 = [scans.T];
bgdata101 = [scans.scn];
bgdata101 = reshape(bgdata101,[],bgNscans101);
% Load data from unit 102
[config,control,scans] = readMrmRetLog('BGsphere102007.csv');
bgNscans102 = length(scans);
bgtstmp102 = [scans.T];
bgdata102 = [scans.scn];
bgdata102 = reshape(bgdata102,[],bgNscans102);

% Load data from unit 103
[config,control,scans] = readMrmRetLog('BGsphere103007.csv');
bgNscans103 = length(scans);
bgtstmp103 = [scans.T];
bgdata103 = [scans.scn];
bgdata103 = reshape(bgdata103,[],bgNscans103);
% Load data from unit 104
[config,control,scans] = readMrmRetLog('BGsphere104007.csv');
bgNscans104 = length(scans);
bgtstmp104 = [scans.T];
bgdata104 = [scans.scn];
bgdata104 = reshape(bgdata104,[],bgNscans104);

% fast time / range uncalibrated
t = linspace(scans(1).Tstrt,scans(1).Tstp,scans(1).Nscn)/1000;   % ns
rng = 3e8*(t-t(1))/2e9; % m

% 101
figure;imagesc(abs(data101));title('101')
figure;imagesc(bgdata101);
BG101 = repmat(mean(bgdata101(:,:),2),1,Nscans101);
data101bg = data101 - BG101;
figure;imagesc(abs(data101bg));title('101')

rngSph101 = rng(153); % range of sphere in unit 105 (m)
rngSph101Real = norm(P101); % real range from unit 105 (m)
rngOffset101 = rngSph101Real - rngSph101; % range offset for unit 101

% 102
figure;imagesc(abs(data102));title('102')
BG102 = repmat(mean(bgdata102(:,:),2),1,Nscans102);
data102bg = data102 - BG102;
figure;imagesc(abs(data102bg));title('102')

rngSph102 = rng(150); % range of sphere in unit 106 (m)
rngSph102Real = norm(P102); % real range from unit 106 (m)
rngOffset102 = rngSph102Real - rngSph102; % range offset for unit 106

% 103
figure;imagesc(abs(data103));title('103')
BG103 = repmat(mean(bgdata103,2),1,Nscans103);
data103bg = data103 - BG103;
figure;imagesc(abs(data103bg));title('103')

rngSph103 = rng(126); % range of sphere in unit 103 (m)
rngSph103Real = norm(P103); % real range from unit 103 (m)
rngOffset103 = rngSph103Real - rngSph103; % range offset for unit 103

% 104
figure;imagesc(abs(data104));title('104')
figure;imagesc(bgdata104)
BG104 = repmat(mean(bgdata104(:,:),2),1,Nscans104);
data104bg = data104 - BG104;
figure;imagesc(abs(data104bg));title('104')

rngSph104 = rng(200); % range of sphere in unit 104 (m)
rngSph104Real = norm(P104); % real range from unit 104 (m)
rngOffset104 = rngSph104Real - rngSph104; % range offset for unit 104


%% Loading Data for tracking
cd('D:\UWB LABS\Tracking gr4');
% Load data from unit 101
[config,control,scans] = readMrmRetLog('B101011.csv');
Nscans101 = length(scans);
tstmp101 = [scans.T];
data101 = [scans.scn];
data101 = reshape(data101,[],Nscans101);
% Load data from unit 102
[config,control,scans] = readMrmRetLog('B102011.csv');
Nscans102 = length(scans);
tstmp102 = [scans.T];
data102 = [scans.scn];
data102 = reshape(data102,[],Nscans102);
% Load data from unit 103
[config,control,scans] = readMrmRetLog('B103011.csv');
Nscans103 = length(scans);
tstmp103 = [scans.T];
data103 = [scans.scn];
data103 = reshape(data103,[],Nscans103);
% Load data from unit 104
[config,control,scans] = readMrmRetLog('B104011.csv');
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

[Nsamp,Nscans] = size(data101_algnd);

figure;imagesc(abs(data101_algnd));title('101')
figure;imagesc(abs(data102_algnd));title('102')
figure;imagesc(abs(data103_algnd));title('103')
figure;imagesc(abs(data104_algnd));title('104')


%% Background removal
BG101 = mean(data101_algnd(:,230:240),2); % select the slow time indexes containing the background data (here 170:230)
BG102 = mean(data102_algnd(:,230:240),2);
BG103 = mean(data103_algnd(:,230:240),2);   
BG104 = mean(data104_algnd(:,230:240),2);


data101_algnd_bg = data101_algnd - repmat(BG101,1,size(data101_algnd,2));
data102_algnd_bg = data102_algnd - repmat(BG102,1,size(data102_algnd,2));  
data103_algnd_bg = data103_algnd - repmat(BG103,1,size(data103_algnd,2)); % remove background
data104_algnd_bg = data104_algnd - repmat(BG104,1,size(data104_algnd,2));


figure;imagesc(abs(data101_algnd_bg));title('101')
figure;imagesc(abs(data102_algnd_bg));title('102')
figure;imagesc(abs(data103_algnd_bg));title('103')
figure;imagesc(abs(data104_algnd_bg));title('104')


path101 = data101_algnd_bg(:,10:190);
path102 = data102_algnd_bg(:,10:190);
path103 = data103_algnd_bg(:,10:190); % select the slow time indexes containing the trajectory data (here 15:115)
path104 = data104_algnd_bg(:,10:190);


path101 = path101/max(abs(path101(:)));
path102 = path102/max(abs(path102(:)));
path103 = path103/max(abs(path103(:))); % normalisation
path104 = path104/max(abs(path104(:))); 


figure;imagesc(abs(path101));title('101')
figure;imagesc(abs(path102));title('102');
figure;imagesc(abs(path103));title('103')
figure;imagesc(abs(path104));title('104')


%% Bandpass filter
b = [0.058918593549 0.003704122993 -0.130605206968 0 0.130605206968 -0.003704122993 -0.058918593549];
a = [1 0.339893240317 1.247471159638 0.315004577848 0.752494992039 0.094346011045 0.145214408359];

path105BPF = filter(b,a,path105);

%% Motion filter FIR4 if backgound subtraction unsuccessful
w = [1 -.6 -.3 -.1];

path104mti = filter(w,1,path104,[],2);
path105 = filter(w,1,path105BPF,[],2);

%% Range estimation
[Nrng,Nscans] = size(path103);
t_slow = 0:.5:.5*(Nscans-1);

% 101
for i1=1:Nscans % Long time
    i2 = 1; % Avoiding picking antenna coupling
    while abs(path101(i2,i1)) < 0.05 && i2 < Nrng-1
        i2 = i2 + 1;
    end
    estimated_range(i1) = rng(i2);
end
figure;plot(estimated_range)
median_range101 = medfilt1(estimated_range,6)+rngOffset101;
figure;imagesc(t_slow,rng,abs(path101));title('101')
hold
plot(t_slow,median_range101,'w')
plot(t_slow,estimated_range,'y')


% 102
for i1=1:Nscans % Long time
    i2 = 1; % Avoiding picking antenna coupling
    while abs(path102(i2,i1)) < 0.04 && i2 < Nrng-1
        i2 = i2 + 1;
    end
    estimated_range(i1) = rng(i2);
end
figure;plot(estimated_range)
median_range102 = medfilt1(estimated_range,6)+rngOffset102;
figure;imagesc(t_slow,rng,abs(path102));title('102')
hold
plot(t_slow,median_range102,'w')

% 103
for i1=1:Nscans % Long time
    i2 = 1; % Avoiding picking antenna coupling
    while abs(path103(i2,i1)) < 0.05 && i2 < Nrng-1    % threshold 0.125
        i2 = i2 + 1;
    end
    estimated_range(i1) = rng(i2);
end
figure;plot(estimated_range)
median_range103 = medfilt1(estimated_range,6)+rngOffset103; % !!!!! test median filter 6 instead of 3
figure;imagesc(t_slow,rng,abs(path103));title('103')
hold
plot(t_slow,median_range103,'y')  % readjust the threshold and coupling index based on this plot

% 104
for i1=1:Nscans % Long time
    i2 = 40; % Avoiding picking antenna coupling
    while abs(path104(i2,i1)) < 0.03 && i2 < Nrng-1
        i2 = i2 + 1;
    end
    estimated_range(i1) = rng(i2);
end
figure;plot(estimated_range)
median_range104 = medfilt1(estimated_range,6)+rngOffset104;
figure;imagesc(t_slow,rng,abs(path104));title('104')
hold
plot(t_slow,median_range104,'w')


%% Target localisation
% 101 - 102
for scani = 1:Nscans
    [a, b] = circcirc(P101(1),P101(2),median_range101(scani),P102(1),P102(2),median_range102(scani));
    x12(scani,:) = a;
    y12(scani,:) = b;

end
figure;plot(x12(:,2),y12(:,2),'o-');title('101 - 102')   

% 101 - 103
for scani = 1:Nscans
    [a, b] = circcirc(P101(1),P101(2),median_range101(scani),P103(1),P103(2),median_range103(scani));
    x13(scani,:) = a;
    y13(scani,:) = b;
end
figure;plot(x13(:,2),y13(:,2));title('101 - 103')  

% 101 - 104
for scani = 1:Nscans
    [a, b] = circcirc(P101(1),P101(2),median_range101(scani),P104(1),P104(2),median_range104(scani));
    x14(scani,:) = a;
    y14(scani,:) = b;
end
figure;plot(x14(:,1),y14(:,1));title('101 - 104')

% 102 - 103
for scani = 1:Nscans
    [a, b] = circcirc(P102(1),P102(2),median_range102(scani),P103(1),P103(2),median_range103(scani));
    x23(scani,:) = a;
    y23(scani,:) = b;
end
figure;plot(x23(:,2),y23(:,2));title('102 - 103') % select coordinates inside the area delimited by the units

% 102 - 104
for scani = 1:Nscans
    [a, b] = circcirc(P102(1),P102(2),median_range102(scani),P104(1),P104(2),median_range104(scani));
    x24(scani,:) = a;
    y24(scani,:) = b;
end
figure;plot(x24(:,2),y24(:,1));title('102 - 104')

% 103 - 104
for scani = 1:Nscans
    [a, b] = circcirc(P103(1),P103(2),median_range103(scani),P104(1),P104(2),median_range104(scani));
    x34(scani,:) = a;
    y34(scani,:) = b;
end
figure;plot(x34(:,2),y34(:,2));title('103 - 104')

%%
x = nanmean([x14(:,1)],2);
y = nanmean([y14(:,1)],2);

figure;plot(smooth(x,3),smooth(y,3),'x-','linewidth',2);axis([-3 3 -3 3]);grid
axis square
hold on
plot(P101(1),P101(2),'ok','markerfacecolor','k');text(P101(1)+.1,P101(2),'101')
plot(P102(1),P102(2),'ok','markerfacecolor','k');text(P102(1)+.1,P102(2),'102')
plot(P103(1),P103(2),'ok','markerfacecolor','k');text(P103(1)+.1,P103(2),'103')
plot(P104(1),P104(2),'ok','markerfacecolor','k');text(P104(1)+.1,P104(2),'104')
    
text(MarkPos(1,1),MarkPos(1,2),'A')
text(MarkPos(2,1),MarkPos(2,2),'B')
text(MarkPos(3,1),MarkPos(3,2),'C')
text(MarkPos(4,1),MarkPos(4,2),'D')
text(MarkPos(5,1),MarkPos(5,2),'E')
text(MarkPos(6,1),MarkPos(6,2),'F')
text(MarkPos(7,1),MarkPos(7,2),'G')
text(MarkPos(8,1),MarkPos(8,2),'H')
text(MarkPos(9,1),MarkPos(9,2),'I')

set(gca,'PlotBoxAspectRatio',[1 1 1])
title('BCFEB trajectory');xlabel('X (m)');ylabel('Y (m)')

plot([E(1),F(1),H(1),D(1),A(1),B(1),C(1)],[E(2),F(2),H(2),D(2),A(2),B(2),C(2)],'color',[0,.5,0],'linewidth',2)


%% Check speed threshold

% 105
estimated_range = zeros(Nscans,1);

i2 = 220; % Avoiding picking antenna coupling
while abs(path105(i2,1)) < 0.17 && i2 < Nrng-1
        i2 = i2 + 1;
end
estimated_range(1) = rng(i2);

for i1=2:Nscans % Slow time
    i2 = 220; % Avoiding picking antenna coupling
    while i2 < Nrng-1
        i2 = i2 + 1;
        if abs(path105(i2,i1)) > 0.19
            rngtmp = rng(i2);
            if abs(rngtmp-estimated_range(i1-1)) < .2
                estimated_range(i1) = rng(i2);
                break
            end
        end
        estimated_range(i1) = estimated_range(i1-1);
    end
 end

median_range105 = medfilt1(estimated_range,6)+rngOffset105;
figure;imagesc(t_slow,rng,abs(path105));title('105')
hold
plot(t_slow,median_range105,'w')
plot(t_slow,estimated_range,'y')

    