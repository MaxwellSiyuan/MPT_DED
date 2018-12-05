clc;clear;close all

load('BasePoint.mat');
load('GridData.mat');
mpc = loadcase('node14case');
define_constants;

% output 
lb = [mpc.gen(1,PMIN);mpc.gen(3,PMIN);mpc.gen(4,PMIN);]*1e3;
ub = [mpc.gen(1,PMAX);mpc.gen(3,PMAX);mpc.gen(4,PMAX);]*1e3;

syms br1 br2 br3 br4 br5 br6 br7 br8 br9 br10 br11 br12 br13
syms th1 th2 th3 th4 th5
syms x1 x2 x3

% br = [br1;br2;br3;br4;br5;br6;br7;br8;br9;br10; br11 ;br12; br13];

% Transmision Line loss
brloss = @(bri,i,j,t) (P0(i,t)^2+2*P0(i,t)*(bri-...
P0(i,t))+Q0(i,t)^2)/V0(j,t)^2 * results.branch(i,BR_R) / 1e3;

MPT_time = zeros(96);
Decision_time = zeros(96);
Option_time = zeros(96);

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

% without the Transmision loss
% balance0 = x1 + x2 + x3 + SolarPower(t)*(1+th1)+WindPower(t)*(1+th2)...
%             -LoadPower(4,t)*(1+th3)-LoadPower(5,t)*...
% 			(1+th3)-LoadPower(6,t)*(1+th3)-LoadPower(3,t)*(1+th3)...
% 			-LoadPower(2,t)*(1+th3)-LoadPower(1,t)*(1+th5)-LoadPower(9,t)...
% 			*(1+th4)-LoadPower(8,t)*(1+th4)-LoadPower(7,t)*(1+th4)...
% 			-LoadPower(14,t)*(1+th5)-LoadPower(13,t)*(1+th5)...
% 			-LoadPower(12,t)*(1+th5)-LoadPower(11,t)*(1+th5)...
% 			-LoadPower(10,t)*(1+th4);
        
br = [br1, br2, br3, br4, br5, br6, br7, br8, br9, br10, br11, br12, br13]-LineLimit;

% Extract the Coeffents 
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

% Power Balance
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
err = 0.1;
err = [err; err;err; err; err;];
Ath = [eye(5);-eye(5)];
bth = [err;err];
% parameters of the objective function
H = diag(K(:,1));
f = K(:,2);
c = sum(K(:,3));

problem = Opt('A',A,'pB',pB,'b',b,'H',H,'f',f,'c',c,'Ath',Ath,'bth',bth);
R = mpt_call_mpqp(problem);

MPT_time(t) = R.stats.solveTime;

Num_regions = R.xopt.Num;

% Generate the Region Table

MaxTable = zeros(Num_regions,5);
MinTable = zeros(Num_regions,5);

P = sdpvar(5,1);
% Optimaize in every dimision
for j=1:5
	Obj = P(j);
	for k=1:Num_regions
		Cons = [R.xopt.Set(k,1).A*P <= R.xopt.Set(k,1).b];
		optimize(Cons,Obj);
		MinTable(k,j)=double(Obj);
		optimize(Cons,-Obj);
		MaxTable(k,j)=double(Obj);
	end
end


% generate the real values of power

th1_ = unifrnd(-err(1),err(1));
th2_ = unifrnd(-err(2),err(2));
th3_ = unifrnd(-err(3),err(3));
th4_ = unifrnd(-err(4),err(4));
th5_ = unifrnd(-err(5),err(5));
th_ = [th1_;th2_;th3_;th4_;th5_;];

RealSolarPower = SolarPower*(1+th1_);
RealWindPower = WindPower*(1+th2_);

RealLoadPower=[LoadPower(1,t)*(1+th5_);
			   LoadPower(2,t)*(1+th3_);
			   LoadPower(3,t)*(1+th3_);
			   LoadPower(4,t)*(1+th3_);
			   LoadPower(5,t)*(1+th3_);
			   LoadPower(6,t)*(1+th3_);
			   LoadPower(7,t)*(1+th4_);
			   LoadPower(8,t)*(1+th4_);
			   LoadPower(9,t)*(1+th4_);
			   LoadPower(10,t)*(1+th4_);
			   LoadPower(11,t)*(1+th5_);
			   LoadPower(12,t)*(1+th5_);
			   LoadPower(13,t)*(1+th5_);
			   LoadPower(14,t)*(1+th5_);];

clc
tic          
          
Realth = ones(Num_regions,1) * th_';
Credit = ((Realth <= MaxTable)&(Realth >= MinTable))';
Credit = sum(Credit);
Credit = find(Credit==5);

if size(Credit,2)==1
	j = Credit
else
	for j = Credit
		if sum(R.xopt.Set(j,1).A*th_ <= R.xopt.Set(j,1).b) == size(th_,1)
			j
			break
		end
	end
end
x = R.mpqpsol.Fi{1,j}*th_+R.mpqpsol.Gi{1,j}
toc

% Conventional one-by-one judgment
tic
for j=1:Num_regions
	if sum(R.xopt.Set(j,1).A*th_ <= R.xopt.Set(j,1).b) == size(th_,1)
		j
		break
	end
end
x = R.mpqpsol.Fi{1,j}*th_+R.mpqpsol.Gi{1,j}
toc
