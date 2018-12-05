clc;clear;close all

load('BasePoint.mat');
load('GridData.mat');
mpc = loadcase('node14case');
define_constants;

% ����������½�
lb = [mpc.gen(1,PMIN);mpc.gen(3,PMIN);mpc.gen(4,PMIN);]*1e3;
ub = [mpc.gen(1,PMAX);mpc.gen(3,PMAX);mpc.gen(4,PMAX);]*1e3;

syms br1 br2 br3 br4 br5 br6 br7 br8 br9 br10 br11 br12 br13
syms th1 th2 th3 th4 th5
syms x1 x2 x3

% br = [br1;br2;br3;br4;br5;br6;br7;br8;br9;br10; br11 ;br12; br13];

% i ��·��� j ��ѹ�ڵ�� t ʱ��
brloss = @(bri,i,j,t) (P0(i,t)^2+2*P0(i,t)*(bri-...
P0(i,t))+Q0(i,t)^2)/V0(j,t)^2 * results.branch(i,BR_R) / 1e3;

t = 50;

% for t=1:96
	nd4=br3-brloss(br3,3,3,t)-LoadPower(4,t)*(1+th3)+x1;
	nd5=br4-brloss(br4,4,3,t)-LoadPower(5,t)*(1+th3)+SolarPower(t)*(1+th1);
	nd6=br5-brloss(br5,5,3,t)-LoadPower(6,t)*(1+th3);
	nd3=br2-brloss(br2,2,2,t)-LoadPower(3,t)*(1+th3)-br3-br4-br5;
	nd2=br1-brloss(br1,1,1,t)-LoadPower(2,t)*(1+th3)-br2;

	nd10=br9-brloss(br9,9,9,t)-LoadPower(10,t)*(1+th4)+WindPower(t)*(1+th2);
	nd9=br8 -brloss(br8,8,7,t)-LoadPower(9,t)*(1+th4)-br9;
 	nd8=br7 -brloss(br7,7,7,t)-LoadPower(8,t)*(1+th4)+x3; % V/th node not included
	nd7=br6 -brloss(br6,6,1,t)-LoadPower(7,t)*(1+th4)-br7-br8+x2;

	nd14=br13-brloss(br13,13,12,t)-LoadPower(14,t)*(1+th5);
	nd13=br12-brloss(br12,12,12,t)-LoadPower(13,t)*(1+th5);
	nd12=br11-brloss(br11,11,11,t)-LoadPower(12,t)*(1+th5)-br13-br12;
	nd11=br10-brloss(br10,10,1,t) -LoadPower(11,t)*(1+th5)-br11;
    nd1 = - LoadPower(1,t)*(1+th5) - br1 - br6 - br10;
    
    [br1,br2,br3,br4,br5,br6,br7,br8,br9,br10,br11,br12,br13]=...
		solve(nd1,nd2,nd3,nd4,nd5,nd6,nd7,nd9,nd10,nd11,nd12,nd13,nd14,...
		br1,br2,br3,br4,br5,br6,br7,br8,br9,br10,br11,br12,br13); 

balance = x1 + x2 + x3 + SolarPower(t)*(1+th1)+WindPower(t)*(1+th2)+...
			-brloss(br3,3,3,t)-brloss(br4,4,3,t)-brloss(br5,5,3,t)...
			-brloss(br2,2,2,t)-brloss(br1,1,1,t)-brloss(br9,9,9,t)...
			-brloss(br8,8,7,t)-brloss(br7,7,7,t)-brloss(br6,6,1,t)...
			-brloss(br13,13,12,t)-brloss(br12,12,12,t)-brloss(br11,11,11,t)...
			-brloss(br10,10,1,t)-LoadPower(4,t)*(1+th3)-LoadPower(5,t)*...
			(1+th3)-LoadPower(6,t)*(1+th3)-LoadPower(3,t)*(1+th3)...
			-LoadPower(2,t)*(1+th3)-LoadPower(1,t)*(1+th5)-LoadPower(9,t)...
			*(1+th4)-LoadPower(8,t)*(1+th4)-LoadPower(7,t)*(1+th4)...
			-LoadPower(14,t)*(1+th5)-LoadPower(13,t)*(1+th5)...
			-LoadPower(12,t)*(1+th5)-LoadPower(11,t)*(1+th5)...
			-LoadPower(10,t)*(1+th4);
balance = -balance;
% balance0 = x1 + x2 + x3 + SolarPower(t)*(1+th1)+WindPower(t)*(1+th2)...
%             -LoadPower(4,t)*(1+th3)-LoadPower(5,t)*...
% 			(1+th3)-LoadPower(6,t)*(1+th3)-LoadPower(3,t)*(1+th3)...
% 			-LoadPower(2,t)*(1+th3)-LoadPower(1,t)*(1+th5)-LoadPower(9,t)...
% 			*(1+th4)-LoadPower(8,t)*(1+th4)-LoadPower(7,t)*(1+th4)...
% 			-LoadPower(14,t)*(1+th5)-LoadPower(13,t)*(1+th5)...
% 			-LoadPower(12,t)*(1+th5)-LoadPower(11,t)*(1+th5)...
% 			-LoadPower(10,t)*(1+th4);
        
br = [br1, br2, br3, br4, br5, br6, br7, br8, br9, br10, br11, br12, br13]-LineLimit;

% ��ȡ����Ԫ��
ALineCons = zeros(Nl,3);
BLineCons = zeros(Nl,5);
bLineCons = zeros(Nl,1);

th = [th1,th2,th3,th4,th5];
x = [x1,x2,x3];

for i=1:Nl
    [coef, symb] = coeffs(br(i));
    for j=1:5
        tmp = find(symb==th(j), 1);
        if isempty(tmp) == 0
            BLineCons(i,j) = coef(tmp);
        end
    end

    for j=1:3
        tmp = find(symb==x(j), 1);
        if isempty(tmp) == 0
            ALineCons(i,j) = coef(tmp);
        end
    end
    
    tmp = find(symb==1, 1);
    if isempty(tmp) == 0
        bLineCons(i) = coef(tmp);
    end
end

BLineCons = -BLineCons;
bLineCons = -bLineCons;

for i = size(ALineCons,1):-1:1
    if sum(ALineCons(i,:) == [0,0,0])==3
        ALineCons(i,:) =[];
        bLineCons(i,:) =[];
        BLineCons(i,:) =[];
    end
end

% ����ƽ��
A_balance = zeros(1,3);
B_balance = zeros(1,5);

[coef, symb] = coeffs(balance);
for j=1:5
    tmp = find(symb==th(j), 1);
    if isempty(tmp) == 0
        B_balance(j) = coef(tmp);
    end
end
for j=1:3
    tmp = find(symb==x(j), 1);
    if isempty(tmp) == 0
        A_balance(j) = coef(tmp);
    end
end
tmp = find(symb==1, 1);
if isempty(tmp) == 0
    b_balance = double(coef(tmp));
end
B_balance = -B_balance;
b_balance = -b_balance;

AOutputCons = [eye(3);-eye(3)];
bOutputCons = [mpc.gen(1,PMAX);mpc.gen(3,PMAX);mpc.gen(4,PMAX);...
    -mpc.gen(1,PMIN);-mpc.gen(3,PMIN);-mpc.gen(4,PMIN);]*1e3;
BOutputCons = zeros(6,5);

A = [ALineCons;AOutputCons;A_balance];%-A_balance];
b = [bLineCons;bOutputCons;b_balance];%-b_balance];
pB =[BLineCons;BOutputCons;B_balance];%-B_balance];

% A = AOutputCons;
% b = bOutputCons;
% pB =BOutputCons;

% forecast error
% solar wind load1 load2 load3 
err = 1;
err = [err; err;err; err; err;];
Ath = [eye(5);-eye(5)];
bth = [err;err];
% parameters of the objective function
H = diag(K(:,1));
f = K(:,2);
c = sum(K(:,3));

problem = Opt('A',A,'pB',pB,'b',b,'H',H,'f',f,'c',c,'Ath',Ath,'bth',bth);

tic
R = mpt_call_mpqp(problem)
toc


% end
% x = [GenSave(1);GenSave(3);GenSave(4);];
% A*x
% b
% A*x<=b



