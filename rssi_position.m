clc,clear,close all

beaconNumber = 12;
beaconName = {12,1};
unitRSSI = 4.42;
passLoss = 10;

%refs = [0,0; 1,1; 2,0; 3,1; 4,0; 5,1; 6,0; 7,1; 8,0; 8,0; 8,2; 7,4];
refs = [1.34,4.76; 6.20,4.76; 6.20,1.54; 1.34,1.54; 2.50,3.16; 4.60,4.13; 1000,1000; 1000,1000; 1000,1000; 1000,1000; 1000,1000; 1000,1000];
% for i = 1:beaconNumber
%    refs(i,1) = refs(i,1) * 5;
%    refs(i,2) = refs(i,2) * 2.1;
% end

for i = 1:beaconNumber;
    beaconName{i} = ['[', char(64 + i), ']'];
end
defaultRSSI = -1000;
StepCount = [];
RSSI = [];
ACC = {};
directVector = [];
acc = zeros(3,1);

cd('honor8');
%fData = load('20170625114525122_xia.txt');
ffid = fopen('20170705161619409_Zhao.txt','r');
%tline = fgetl(ffid);

while feof(ffid) == 0
    rssi = ones(beaconNumber,1) * defaultRSSI;
    S = regexp(fgetl(ffid), ',', 'split');
    [temp,len] = size(S);
    flag = str2double(S{3});
    if flag == 10   %加速度
        % 主成分分析
        a = [str2double(S{4});str2double(S{5});str2double(S{6})];
        acc = [acc,a];
    elseif flag == 19   %计步器、iBeacon
        StepCount = [StepCount, str2double(S{4})];        
        for i = 1:beaconNumber;
            ss = strfind(S, beaconName{i});
            for m = 6:len;
                s = ss{1,m};            
                if s > 0
                    str = strtok(S{m},beaconName{i});
                    rssi(i) = str2num(str);
                    break;
                end
            end
        end
        RSSI = [RSSI,rssi];
        % 计算方向矢量
        sigma = cov(acc');  % 协方差矩阵
        [v,d]=eig(sigma);   % 特征向量与特征值 
        eigenValue = max(max(d));   % 最大特征值
        eigenVector = v(:,find(max(d) == eigenValue));  % 最大特征值对应特征向量
        dv = [eigenVector(1);eigenVector(2)];  % 最大特征值对应特征向量在x,y平面上的投影
        directVector = [directVector,dv];
        % ACC重新开始记录
        ACC = [ACC, acc];
        acc = [];
    end    
end

%% 计算距离
[m, n] = size(RSSI);
dists = RSSI;
for i = 1:m
    for j = 1:n
        rssi = abs(RSSI(i,j));
        dist =10^((rssi - unitRSSI)/(10 * passLoss));
        %dist =10^((abs(RSSI(i,j)) - unitRSSI)/(10 * 2));
        dists(i,j) = dist;
    end
end


%% 最小二乘法
coors = [];
for col = 1:n
    A = [];
    b = [];
   for row = 1:m
      if dists(row,col) < 1000
          A = [A; 1, -2 * refs(row,1), -2 * refs(row,2)];
          b = [b; dists(row,col)^2 - refs(row,1)^2 - refs(row,2)^2];
      end
   end
   X = (A' * A)^-1 * A' * b;
   coors = [coors,X];
end

cd('..');

figure(1)
hold on
acc = ACC{1,20};
scatter3(acc(1,:),acc(2,:),acc(3,:),'ro');
grid on

acc = acc';
center = mean(acc);
sigma = cov(acc);
[v,d]=eig(sigma);
eigenValue = max(max(d));
eigenVector = v(:,find(max(d) == eigenValue));
dv = [eigenVector(1),eigenVector(2)];
nor = norm(v,1);

%[coef,score,latent,t2] = pca(acc);
plot3([0 d(1,1) * v(1,1)],[0 d(1,1) * v(2,1)],[0 d(1,1) * v(3,1)]);
plot3([0 d(2,2) * v(1,2)],[0 d(2,2) * v(2,2)],[0 d(2,2) * v(3,2)]);
plot3([0 d(3,3) * v(1,3)],[0 d(3,3) * v(2,3)],[0 d(3,3) * v(3,3)]);

plot3([0 eigenValue * dv(1)],[0 eigenValue * dv(2)],[0 0]);
xlabel('x');
ylabel('y');
zlabel('z');
hold off
% a=100;
% windowSize = a;
% y=filter(ones(1,windowSize)/windowSize,1,dists(1,:));
% hold on
% plot(y);
% plot(dists(1,:));
% hold off

% cc = find(coors(2,:) < 10 & coors(2,:) > 1);
% 
% scatter(coors(2,cc), coors(3,cc));

%hist(coors(2,cc))

% stdRSSI = mean(RSSI(13,:));
% 
% hist(RSSI(13,:));
% 
% dist13 = mean(dists(13,:));