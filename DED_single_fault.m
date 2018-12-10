clc;clear;close all

load('BasePoint.mat');
load('GridData.mat');
mpc = loadcase('node14case');
define_constants;

% forecast error
err = 0.2;
err = ones(8,1)*err;

% output 
lb = [mpc.gen(1,PMIN);mpc.gen(3,PMIN);mpc.gen(4,PMIN);]*1e3;
ub = [mpc.gen(1,PMAX);mpc.gen(3,PMAX);mpc.gen(4,PMAX);]*1e3;

syms br1 br2 br3 br4 br5 br6 br7 br8 br9 br10 br11 br12 br13
syms th1 th2 th3 th4 th5 th6 th7 th8
syms x1 x2 x3

% br = [br1;br2;br3;br4;br5;br6;br7;br8;br9;br10; br11 ;br12; br13];

% Transmision Line loss
brloss = @(bri,i,j,t) (P0(i,t)^2+2*P0(i,t)*(bri-...
P0(i,t))+Q0(i,t)^2)/V0(j,t)^2 * results.branch(i,BR_R) / 1e3;

RealSolarPower = zeros(1,96);
RealWindPower =  zeros(1,96);
RealLoadPower = zeros(14,96);

% t = 50;
X_mpt_con_array = zeros(3,96);
X_mpt_con_array_f = zeros(2,96);
X_mpt_new_array = zeros(3,96);
X_yalmip_array = zeros(3,96);

MPT_time = zeros(1,96);
time_mpt_con_array = zeros(1,96);
time_mpt_new_array = zeros(1,96);
time_yalmip_array = zeros(1,96);

t = 53;
for t=1:96
    t
    
     % generate the real values of power
    rand('seed',t)
    th1_ = (rand(1,1)*2 - 1)*err(1);
    th2_ = (rand(1,1)*2 - 1)*err(2);
	th3_ = (rand(1,1)*2 - 1)*err(3);
	th4_ = (rand(1,1)*2 - 1)*err(4)/2;
	th5_ = (rand(1,1)*2 - 1)*err(5)/1;
	th6_ = (rand(1,1)*2 - 1)*err(6)/4;
	th7_ = (rand(1,1)*2 - 1)*err(7)/4;
	th8_ = (rand(1,1)*2 - 1)*err(8)/1;
    th_ = [th1_;th2_;th3_;th4_;th5_;th6_;th7_;th8_;];
    
    syms br1 br2 br3 br4 br5 br6 br7 br8 br9 br10 br11 br12 br13
	nd4=br3-brloss(br3,3,3,t)-LoadPower(4,t)*(1+th7)+x1;
	nd5=br4-brloss(br4,4,3,t)-LoadPower(5,t)*(1+th7)+SolarPower(t)*(1+th1);
	nd6=br5-brloss(br5,5,3,t)-LoadPower(6,t)*(1+th3);
	nd3=br2-brloss(br2,2,2,t)-LoadPower(3,t)*(1+th3)-br3-br4-br5;
	nd2=br1-brloss(br1,1,1,t)-LoadPower(2,t)*(1+th4)-br2;

	nd10=br9-brloss(br9,9,9,t)-LoadPower(10,t)*(1+th7)+WindPower(t)*(1+th2);
	nd9=br8 -brloss(br8,8,7,t)-LoadPower(9,t)*(1+th5)-br9; 
 	nd8=br7 -brloss(br7,7,7,t)-LoadPower(8,t)*(1+th7)+x3; 
	nd7=br6 -brloss(br6,6,1,t)-LoadPower(7,t)*(1+th7)-br7-br8+x2; % V/th node not included

	nd14=br13-brloss(br13,13,12,t)-LoadPower(14,t)*(1+th6);
	nd13=br12-brloss(br12,12,12,t)-LoadPower(13,t)*(1+th6);
	nd12=br11-brloss(br11,11,11,t)-LoadPower(12,t)*(1+th6)-br13-br12;
	nd11=br10-brloss(br10,10,1,t) -LoadPower(11,t)*(1+th6)-br11;
    nd1 = - LoadPower(1,t)*(1+th8) - br1 - br6 - br10;
    
    [br1,br2,br3,br4,br5,br6,br7,br8,br9,br10,br11,br12,br13]=...
		solve(nd1,nd2,nd3,nd4,nd5,nd6,nd8,nd9,nd10,nd11,nd12,nd13,nd14,...
		br1,br2,br3,br4,br5,br6,br7,br8,br9,br10,br11,br12,br13); 

    balance = x1 + x2 + x3 + SolarPower(t)*(1+th1)+WindPower(t)*(1+th2)+...
                -brloss(br3,3,3,t)-brloss(br4,4,3,t)-brloss(br5,5,3,t)...
                -brloss(br2,2,2,t)-brloss(br1,1,1,t)-brloss(br9,9,9,t)...
                -brloss(br8,8,7,t)-brloss(br7,7,7,t)-brloss(br6,6,1,t)...
                -brloss(br13,13,12,t)-brloss(br12,12,12,t)-brloss(br11,11,11,t)...
                -brloss(br10,10,1,t)-LoadPower(1,t)*(1+th8)-LoadPower(2,t)*(1+th4)-...
				LoadPower(3,t)*(1+th3)-LoadPower(4,t)*(1+th7)-...
				LoadPower(5,t)*(1+th7)-LoadPower(6,t)*(1+th3)-...
				LoadPower(7,t)*(1+th7)-LoadPower(8,t)*(1+th7)-...
				LoadPower(9,t)*(1+th5)-LoadPower(10,t)*(1+th7)-...
				LoadPower(11,t)*(1+th6)-LoadPower(12,t)*(1+th6)-...
				LoadPower(13,t)*(1+th6)-LoadPower(14,t)*(1+th6);

    balance = -balance;

    br = [br1, br2, br3, br4, br5, br6, br7, br8, br9, br10, br11, br12, br13]-LineLimit;

    % Extract the Coeffents 
    ALineCons = zeros(Nl,3);
    BLineCons = zeros(Nl,8);
    bLineCons = zeros(Nl,1);

    th = [th1,th2,th3,th4,th5,th6,th7,th8];
    x = [x1,x2,x3];

    for i=1:Nl
        [coef, symb] = coeffs(br(i));
        for j=1:8
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
    B_balance = zeros(1,8);

    [coef, symb] = coeffs(balance);
    for j=1:8
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
    BOutputCons = zeros(6,8);
    
    RampLimit = [80;140;110]/4;
    ARampCons = [eye(3);-eye(3)];
    if t~=1
        
        bRampCons = [RampLimit+X_mpt_con_array(:,t-1);RampLimit-X_mpt_con_array(:,t-1)];
    else
        bRampCons = [RampLimit+[55;65;53];RampLimit-[55;65;53]];
    end
    BRampCons = zeros(6,8);

%     A = [ALineCons;AOutputCons;A_balance;ARampCons];
%     b = [bLineCons;bOutputCons;b_balance;bRampCons];
%     pB =[BLineCons;BOutputCons;B_balance;BRampCons];

    A = [ALineCons;AOutputCons;A_balance;];%ARampCons];
    b = [bLineCons;bOutputCons;b_balance;];%bRampCons];
    pB =[BLineCons;BOutputCons;B_balance;];%BRampCons];

   
%     err = unifrnd(0,0.6,5,1);
    Ath = [eye(8);-eye(8)];
    bth = ones(size(Ath,1),1);
    % parameters of the objective function
    H = diag(K(:,1));
    f = K(:,2);
    c = sum(K(:,3));

%     problem = Opt('A',A,'pB',pB,'b',b,'H',H,'f',f,'c',c,'Ath',Ath,'bth',bth);
%     R = mpt_call_mpqp(problem);
%     
%     Num_regions = R.xopt.Num;
%     for j=1:Num_regions
%         if sum(R.xopt.Set(j,1).A*th_ <= R.xopt.Set(j,1).b) == size(R.xopt.Set(j,1).A,1)
%             j
%             break
%         end
%     end
%     X_mpt_con_array(:,t) = R.mpqpsol.Fi{1,j}*th_+R.mpqpsol.Gi{1,j};
%     R.mpqpsol.Fi{1,j}*th_+R.mpqpsol.Gi{1,j}

    Af = A(:,1:2);
    pBf = pB;
    bf = b;
    for i = size(Af,1):-1:1
        if sum(Af(i,:) == zeros(1,2))==2
            Af(i,:)=[];
            pBf(i,:)=[];
            bf(i,:)=[];
        end
    end
    Hf = diag(K(1:2,1));
    ff = K(1:2,2);
    cf = sum(K(1:2,3));
    
    problemf = Opt('A',Af,'pB',pBf,'b',bf,'H',Hf,'f',ff,'c',cf,'Ath',Ath,'bth',bth);
    Rf = mpt_call_mpqp(problemf);
    
    Num_regions_f = Rf.xopt.Num;
    for j=1:Num_regions_f
        if sum(Rf.xopt.Set(j,1).A*th_ <= Rf.xopt.Set(j,1).b) == size(Rf.xopt.Set(j,1).A,1)
    		j
            break
        end
    end
    X_mpt_con_array_f(:,t) = Rf.mpqpsol.Fi{1,j}*th_+Rf.mpqpsol.Gi{1,j};
    Rf.mpqpsol.Fi{1,j}*th_+Rf.mpqpsol.Gi{1,j}
end

save('solution0.2f.mat')

figure(1)
plot(time_mpt_new_array)
hold on
plot(time_mpt_con_array)
plot(time_yalmip_array)
title('The Comparasion of Computing time')
legend('Speed up Search Method with MPT',...
    'Conventional Search Method with MPT',...
   'Optimization with Yalmip')

figure
plot(X_yalmip_array(1,:))
hold on
plot(X_yalmip_array(2,:))
plot(X_yalmip_array(3,:))
title('Output Power of Dispatchble DGs with Yalmip')
legend('Micro Turbine','Diesel Engine','Fuel Cell')

figure
plot(X_mpt_con_array(1,:))
hold on
plot(X_mpt_con_array(2,:))
plot(X_mpt_con_array(3,:))
title('Output Power of Dispatchble DGs wiht MPT')
legend('Micro Turbine','Diesel Engine','Fuel Cell')

figure(4)
plot(X_mpt_new_array(1,:))
hold on
plot(X_mpt_new_array(2,:))
plot(X_mpt_new_array(3,:))
title('Output Power of Dispatchble DGs with New Searching method')
legend('Micro Turbine','Diesel Engine','Fuel Cell')

figure(5)
plot(sum(X_mpt_con_array)+RealSolarPower+RealWindPower)
hold on
plot(sum(RealLoadPower))


