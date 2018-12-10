clear;clc;close all

load('solution0.2.mat')
set(0,'defaultfigurecolor','w')
x1 = datenum('00:00');
x2 = datenum('24:00');
x = linspace(x1,x2,96);

figure
plot(x,sum(RealLoadPower),'r')
datetick('x',15)
xlabel('Time / hour:min','Fontname', 'Times New Roman','FontSize',13)
ylabel('Real-time Total Load Power / kW','Fontname', 'Times New Roman','FontSize',13)
set(gca,'Fontname', 'Times New Roman','FontSize',13);
set(gcf,'Units','centimeter','Position',[5,5,21,7]);

%%
close all

figure
plot(x,sum(LoadPower),'r');
datetick('x',15)
xlabel('Time / hour:min','Fontname', 'Times New Roman','FontSize',14)
ylabel('Forecast Total Load Power / kW','Fontname', 'Times New Roman','FontSize',14)
set(gcf,'Units','centimeter','Position',[5,5,21,7]);
set(gca,'Fontname', 'Times New Roman','FontSize',13);

figure
plot(x,WindPower,'r');
datetick('x',15)
xlabel('Time / hour:min','Fontname', 'Times New Roman','FontSize',14)
ylabel('Forecast Wind Power / kW','Fontname', 'Times New Roman','FontSize',14)
set(gcf,'Units','centimeter','Position',[5,5,21,7]);
set(gca,'Fontname', 'Times New Roman','FontSize',13);

figure
plot(x,SolarPower,'r');
datetick('x',15)
xlabel('Time / hour:min','Fontname', 'Times New Roman','FontSize',14)
ylabel('Forecast Solar Power / kW','Fontname', 'Times New Roman','FontSize',14)
set(gcf,'Units','centimeter','Position',[5,5,21,7]);
set(gca,'Fontname', 'Times New Roman','FontSize',13);

%%
close all
figure
plot(x,X_mpt_con_array(1,:))
hold on
plot(x,X_mpt_con_array(2,:))
plot(x,X_mpt_con_array(3,:))
datetick('x',15)
xlabel('Time / hour:min','Fontname', 'Times New Roman','FontSize',15)
ylabel('Forecast Solar Power / kW','Fontname', 'Times New Roman','FontSize',15)
set(gcf,'Units','centimeter','Position',[1,1,21,12]);
l1 = legend('Micro Turbine','Diesel Engine','Fuel Cell');
set(l1,'Fontname', 'Times New Roman','FontSize',12)
set(gca,'FontSize',12,'Fontname', 'Times New Roman')
