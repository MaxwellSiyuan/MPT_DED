clc;clear;close all

mpc = loadcase('node14case');
define_constants;

WindPower=[8.70,9.50 ,10.27,16.66 ,7.22 ,4.91,14.66 ,26.56 ,...
    20.88 ,17.85 ,12.80 ,18.65 ,14.35 ,10.35 ,8.26,5.71 ,8.44,...
    16.87 ,20.75 ,23.17 ,28.15 ,21.31 ,10.07 ,7.58]*5;
SolarPower=[zeros(7,1);1;3;6;9;12;13;12.9;10.5;9;5;3;1;zeros(5,1)]'*8;
Load=[140;150;155;160;165;170;175;190;240;220;240;250;240;220;200; ...
    180;170;185;200;240;225;190;160;145]'*1.5;

Nd = 14;
Nl = 13;
LineLimit = 140;

hour = 1:24;
Tnode = 0.25:0.25:24;
WindPower = interp1 (hour,WindPower,Tnode,'pchip');
SolarPower = interp1 (hour,SolarPower,Tnode,'pchip');
Load = interp1 (hour,Load,Tnode,'pchip');

LoadPower = zeros(Nd,96);
LoadHeavy = Load *3/5 / 6;
LoadLight = Load *2/5 / 8;
for i=1:14
    if i>=4 && i<=9
        LoadPower(i,:) = LoadHeavy + unifrnd(-5,2,1,96);
    else
        LoadPower(i,:) = LoadLight + unifrnd(-3,1,1,96);
    end
%     plot(Tnode,LoadPower(i,:))
%     hold on
end

LoadQ = LoadPower ./ unifrnd(4,8,14,96);

LoadSum = sum(LoadPower);

% DG1 (MT) node4
% DG2 (FC) node7
% DG3 (DE) node8
DGout = zeros(3,96); 
K = [0.0001 0.0417 0.023
     0.0001 0.0440 0.021
     0.0001 0.0435 0.019];
options=sdpsettings('solver', 'cplex');
sdpvar x1 x2 x3
Obj = K(1,1)*x1^2 + K(1,2)*x1 + K(1,3) +...
      K(2,1)*x2^2 + K(2,2)*x2 + K(2,3) +...
      K(3,1)*x3^2 + K(3,2)*x3 + K(3,3) ;

br = sdpvar(13,96);
br(3,:) = x1 - LoadPower(4,:);
br(4,:) = SolarPower - LoadPower(5,:);
br(5,:) = - LoadPower(6,:);
br(2,:) = sum(br(3:5,:))-LoadPower(3,:);
br(1,:) = br(2,:) - LoadPower(2,:);
br(9,:) = WindPower - LoadPower(10,:);
br(8,:) = br(9,:) - LoadPower(9,:);
br(7,:) = x3 - LoadPower(8,:);
br(6,:) = br(8,:) + br(7,:) + x2 - LoadPower(7,:);
br(13,:) = - LoadPower(14,:);
br(12,:) = - LoadPower(13,:);
br(11,:) = br(12,:) + br(13,:) - LoadPower(12,:);
br(10,:) = br(11,:) - LoadPower(11,:);

% 计算机组经济调度出力（不考虑网损）
for i=1:96
    Cons =[ x1+x2+x3+WindPower(i)+SolarPower(i)==LoadSum(i);
            % 线路约束
            br(:,i) <= LineLimit * ones(13,1);
            -br(:,i) <= LineLimit * ones(13,1);
            x1>=mpc.gen(1,PMIN)*1e3;x1<=mpc.gen(1,PMAX)*1e3;
            x2>=mpc.gen(3,PMIN)*1e3;x2<=mpc.gen(3,PMAX)*1e3;
            x3>=mpc.gen(4,PMIN)*1e3;x3<=mpc.gen(4,PMAX)*1e3;];
    optimize(Cons,Obj,options);
    
    DGout(1,i) = double(x1);
    DGout(2,i) = double(x2);
    DGout(3,i) = double(x3);
end

% plot DG output
figure(2)
subplot(3,1,1)
plot(Tnode,DGout(1,:))
hold on
plot(Tnode,ones(1,96)*mpc.gen(1,PMIN)*1e3,'r-')
plot(Tnode,ones(1,96)*mpc.gen(1,PMAX)*1e3,'r-')
title('Output Power of Micro Turbine')
subplot(3,1,2)
plot(Tnode,DGout(2,:))
hold on
plot(Tnode,ones(1,96)*mpc.gen(3,PMIN)*1e3,'r-')
plot(Tnode,ones(1,96)*mpc.gen(3,PMAX)*1e3,'r-')
title('Output Power of Fuel Cell')
subplot(3,1,3)
plot(Tnode,DGout(3,:))
hold on
plot(Tnode,ones(1,96)*mpc.gen(4,PMIN)*1e3,'r-')
plot(Tnode,ones(1,96)*mpc.gen(4,PMAX)*1e3,'r-')
title('Output Power of Diesel Engine')

clear Obj Cons x1 x2 x3 ans options br
save('GridData.mat')