clc;clear;close all

load('GridData.mat')

P0 = zeros(Nl,96);
Q0 = zeros(Nl,96);
V0 = zeros(Nd,96);

genP = zeros(5,96);
genQ = zeros(5,96);

mpc = loadcase('node14case');
define_constants;

for i=1:96
    mpc.bus(:,[PD QD]) = [LoadPower(:,i) LoadQ(:,i)];
    mpc.bus(:, [PD, QD]) = mpc.bus(:, [PD, QD]) / 1e3;
    
    mpc.gen(:,PG) = [DGout(1,i);SolarPower(i);DGout(2,i);DGout(3,i);WindPower(i)];
    mpc.gen(:,PG) =  mpc.gen(:,PG) / 1e3 ;
    
    mpc.gen(:,QMAX) = [DGout(1,i);SolarPower(i);DGout(2,i);DGout(3,i);WindPower(i)];
    mpc.gen(:,QMAX) =  mpc.gen(:,QMAX) / 1e3 ;
    mpc.gen(:,QMIN) = -mpc.gen(:,QMAX);
    
    
    results = runpf(mpc);
    
    V0(:,i) = results.bus(:,VM);
    P0(:,i) = results.branch(:,PF) * 1e3;
    Q0(:,i) = results.branch(:,QF) * 1e3;
    
    genP(:,i) = results.gen(:,PG) * 1e3;
    genQ(:,i) = results.gen(:,QG) * 1e3;
    
    if i==50
        GenSave = results.gen(:,PG) * 1e3;
    end
    
end

save('BasePoint.mat','V0','P0','Q0','results','GenSave')

figure(1)
for i=1:5
    plot(genP(i,:))
    hold on
end
legend('MT','PV','DE','FC','WT')
title('Output of active power')

figure(2)
for i=1:5
    plot(genQ(i,:))
    hold on
end
legend('MT','PV','DE','FC','WT')
title('Output of reactive power')

figure(3)
for i= 1:13
    plot(P0(i,:))
    hold on
end
legend('1','2','3','4','5','6','7','8','9','10','11','12','13')
title('Branch Flow')

figure(4)
plot(sum(genP))
hold on
plot(LoadSum)
legend('GE','LD')

% plot wind solar power data
% subplot(3,1,1)
% plot(Tnode,WindPower)
% axis([0,24,0,inf])
% subplot(3,1,2)
% plot(Tnode,SolarPower)
% axis([0,24,0,inf])
% subplot(3,1,3)
% plot(Tnode,Load)
% hold on 
% plot(Tnode,SolarPower+WindPower)
% axis([0,24,0,inf])
% 
% figure(2)
% plot(Tnode,Load,'b')
% hold on
% plot(Tnode,WindPower+SolarPower+260,'r')
% plot(Tnode,WindPower+SolarPower+101,'g')






