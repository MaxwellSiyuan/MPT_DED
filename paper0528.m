clc
clear
% 24 hour power
windpower=[0.288;0.277;0.271;0.272;0.27;0.269;0.262;0.258;0.255;0.252;0.255;0.262;0.272;0.281;0.293;0.299;0.298;0.3;0.305;0.302;0.298;0.296;0.292;0.29];%1天的出力
Solar=[zeros(7,1);1;3;6;9;12;13;12.9;12.5;11;8;3;1;zeros(5,1)];
load=[140;150;155;160;165;170;175;190;240;220;240;250;240;220;200;180;170;185;200;240;225;190;160;145];
%  l= unifrnd(0,10,8,1);
% 发电机成本参数
K = [1165.06365581861,20.1650102578050,0.924920128230088;
    1503.46951881223,20.3902958223262,0.595329829887496;
    1127.37929730191,20.4561626228323,0.985326047262670;
    1190.31425094050,20.1901389926730,0.653150336520274;
    835.590089646348,20.4403593349635,0.571775051567062;
    1216.01943792928,20.1687841750247,0.754287260953161;
    1083.19888501040,20.0637903383511,0.577200288656170];
% K(:, 1) = unifrnd(600, 1530, 7, 1);
% K(:, 2) = unifrnd(20, 20.5, 7, 1);
% K(:, 3) = unifrnd(0.5, 1, 7, 1);
% 【？】K是什么 7*3
H=diag(K(:, 3)); % 成本二次项系数
f=K(:, 2); % 成本一次项系数
c=sum(K(:, 1)); % 成本常数

j=8; % 时间;上午 08:00

l= unifrnd(6,10,19,1);         %负荷：均匀分布[6,10] 19*1
windpower1=unifrnd(9,10,3,1);  %风能：均匀分布[9,10]  3*1
solar1=unifrnd(9,10,3,1);      %光能：均匀分布[9,10]  3*1

windpower2=windpower1*windpower(j)*17000/sum(windpower1); % 【？】随机构造一个风能
solar2=solar1*Solar(j)*220/sum(solar1); % 【？】随机构造一个光能

L=l*load(j)*50/(sum(l));  % 【？】随机构造一个负荷
windpower3=0.6*windpower2;% 【？】
windpower4=1.4*windpower2;% 【？】
solar3=0.6*solar2;  % 【？】
solar4=1.4*solar2;  % 【？】
L1=0.95*L;
L2=1.05*L;
    
Ath=[eye(6);-eye(6)];
bth=[windpower4(1);windpower4(2);windpower4(3);solar4(1)+0.01;solar4(2)+0.01;solar4(3)+0.01; % 【？】 +0.01？
    -windpower3(1);-windpower3(2);-windpower3(3);-solar3(1);-solar3(2);-solar3(3);];

% bth=[windpower4(1);windpower4(2);windpower4(3);solar4(1);solar4(2);solar4(3); % 【？】 +0.01？
%     -windpower3(1);-windpower3(2);-windpower3(3);-solar3(1);-solar3(2);-solar3(3);];
%【？】solar+0.01原因 ： 松弛以求解
%把负载的不确定性考虑进去
%做一个微网算例
%随机化负载
% M=140;
% l= unifrnd(0,10,8,1);
% L=l*M/(sum(l));
%可控机组的出力上下界【x:机组的出力】
lb=100*ones(7,1);
ub=2000*ones(7,1);
%线路约束
% 31*7 0-1矩阵
A1=[zeros(1,7);0 1,zeros(1,5); 0 1,zeros(1,5);0 1,zeros(1,5);
    1,1,zeros(1,5);1,1,zeros(1,5);1,1,zeros(1,5);1,1,zeros(1,5);
    1,1,zeros(1,5);1,1,zeros(1,5);1,1,zeros(1,5);1,1,zeros(1,5);
    zeros(1,7);zeros(1,3),1,zeros(1,3);zeros(1,3),1,zeros(1,3);zeros(1,3),1,zeros(1,3);
    zeros(1,3),1,zeros(1,3);zeros(1,3),1,zeros(1,3);zeros(1,3),1,zeros(1,3);zeros(1,2),1,1,zeros(1,3);
    ones(1,4),zeros(1,3);ones(1,4),zeros(1,3);ones(1,4),zeros(1,3);zeros(1,6),1;
    zeros(1,6),1;zeros(1,6),1;ones(1,4),zeros(1,2),1;zeros(1,5),1,0;
    zeros(1,4),1,1,0;zeros(1,4),1,1,0;zeros(1,4),1,1,0;];
% 31*25 矩阵
pB1=[zeros(1,16),-1,zeros(1,8);zeros(1,16),-1,zeros(1,8);zeros(1,15),-1,-1,zeros(1,8);zeros(1,14),-1,-1,-1,zeros(1,8);
    zeros(1,14),-1,-1,-1,zeros(1,8);zeros(1,3),1,zeros(1,10),-1,-1,-1,zeros(1,8);zeros(1,3),1,zeros(1,9),-1,-1,-1,-1,zeros(1,8);zeros(1,3),1,zeros(1,8),-1,-1,-1,-1,-1,zeros(1,8);
    zeros(1,2),1,1,zeros(1,8),-1,-1,-1,-1,-1,zeros(1,8);zeros(1,2),1,1,zeros(1,7),-1,-1,-1,-1,-1,-1,zeros(1,8); 0,1,1,1,zeros(1,7),-1,-1,-1,-1,-1,-1,zeros(1,8); 0,1,1,1,zeros(1,6),-1,-1,-1,-1,-1,-1,-1,zeros(1,8); 
    zeros(1,21),-1,zeros(1,3);zeros(1,21),-1,zeros(1,3);zeros(1,20),-1,-1,zeros(1,3);zeros(1,19),-1,-1,-1,zeros(1,3);
    zeros(1,4),1,zeros(1,14),-1,-1,-1,zeros(1,3);zeros(1,4),1,zeros(1,13),-1,-1,-1,-1,zeros(1,3);zeros(1,4),1,zeros(1,12),-1,-1,-1,-1,-1,zeros(1,3);zeros(1,4),1,zeros(1,12),-1,-1,-1,-1,-1,zeros(1,3);
    0,ones(1,4),zeros(1,4),-ones(1,8),-ones(1,5),zeros(1,3);0,ones(1,4),zeros(1,3),-ones(1,9),-ones(1,5),zeros(1,3);ones(1,5),zeros(1,3),-ones(1,9),-ones(1,5),zeros(1,3);zeros(1,25);
    zeros(1,24),-1;zeros(1,5),1,zeros(1,18),-1;ones(1,6),0,-1,-ones(1,9),-ones(1,5),zeros(1,2),-1;zeros(1,25);
    zeros(1,25);zeros(1,23),-1,0;zeros(1,22),-1,-1,0;];
A=[A1;-A1;ones(1,7);-ones(1,7)];   % 64*7 矩阵
pB=[-pB1;pB1;-ones(1,6),ones(1,19);ones(1,6),-ones(1,19);]; % 64*25 矩阵
b2=pB(1:64,7:25);
pB=pB(1:64,1:6);   % 64*6 矩阵
b=[4000*ones(31,1);4000*ones(31,1);0.01;0];  % 64*1 矩阵
b=b+b2*L;
problem = Opt('A',A,'pB',pB,'b',b,'lb',lb,'ub',ub,'H',H,'f',f,'c',c,'Ath',Ath,'bth',bth)%'Ae',Ae,'be',be,'pE',E)
tic
R = mpt_call_mpqp(problem)
toc
e=[];
for(i=1:3)
e(i,1)=unifrnd(windpower3(i),windpower4(i));
end
for(i=4:6)
e(i,1)=unifrnd(solar3(i-3),solar4(i-3));
end
% for(i=7:25)
%     e(i,1)=unifrnd(L1(i-6),L2(i-6));
% end
%加一个判断参量是否在区域内的部分
% r1=e<(20*ones(12,1));
% r2=e>(zeros(12,1));
% r3=r1.*r2;
% if(numel(r3,r3==0))
%     disp('The input is wrong!');
% else
% Ai=[];
% Bi=[];
[u1,u2]=size(R.mpqpsol.Fi);
tic
for i=1:u2
    r4=(R.xopt.Set(i,1).A)*e<=(R.xopt.Set(i,1).b);
    if(numel(r4, r4==0)==0)
        i
        break
    else
    end
end
% while(numel(r4, r4==1)&&i<=54)
%     i=i+1;
%     r4=(R.xopt.Set(i,1).A)*e>(R.xopt.Set(i,1).b);
% end
x=R.mpqpsol.Fi{1,i}*e+R.mpqpsol.Gi{1,i}
result=0.5*x'*H*x+f'*x+c
toc
%验证正确性
tic
P=sdpvar(7,1);
g1=0.5*P'*H*P+f'*P+c;
H1=[A*P<=b+pB*e;150*ones(7,1)<=P<=2000*ones(7,1)];
sol=solvesdp(H1,g1);
P=double(P)
g1
toc
