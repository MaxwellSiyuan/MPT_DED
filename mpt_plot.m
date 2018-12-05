clc
clear
close all

A = [2.364728 2.812932 1.255367
    0.237276 1.476349 1.116919
    0.336999 1.366638 1.570896
    2.788394 1.974039 2.897298
    2.869431 0.621541 1.354722
    1.257747 0.807349 0.070875
    0.761589 0.245812 0.483367
    0.274586 1.160139 1.593949
    1.900577 0.48503 0.376802
    1.264104 1.543067 2.520901
    0.204818 0.62668 2.765897
    0.527312 2.819617 0.163602];
b = [3.430588
    2.675002
    2.488142
    4.768894
    3.72822
    2.41143
    1.917835
    2.016847
    2.924795
    0.024295
    0.551397
    0.010254];
c = [3.811557 3.32477 4.593507]';
F = [-4.81047 -7.86886
    -1.18619 6.299455
    -8.74241 -0.05091
    -9.77271 -8.30604
    2.687734 7.96541
    6.896508 1.193105
    -6.62176 1.811549
    -3.903 -5.32411
    3.589349 -0.1948
    -8.58502 -9.7836
    4.681477 -7.74365
    5.33831 -1.55671];
Q = [1.96 0.63 0.4
    0.63 5.4 1.99
    0.328 0.212 1.26];
Ath = [eye(2); -eye(2)]; bth = 6*ones(4,1); 
% problem = Opt('A',A,'pB',F,'b',b,'lb',lb,'ub',ub,'H',Q,'f',c,'c',F,'Ath',Ath,'bth',bth)
problem = Opt('A',A,'pB',F,'b',b,'H',Q,'f',c,'Ath',Ath,'bth',bth)
figure('color','w')
solution = problem.solve
subplot(2,2,1)
solution.xopt.fplot('primal','position',1)
xlabel('\theta_{1} ')
ylabel('\theta_{2} ')
zlabel('x_1 ')
view(-25,65)

subplot(2,2,2)
solution.xopt.fplot('primal','position',2)
xlabel('\theta_{1} ')
ylabel('\theta_{2} ')
zlabel('x_2 ')
view(-25,65)

subplot(2,2,3)
solution.xopt.fplot('primal','position',3)
xlabel('\theta_{1} ')
ylabel('\theta_{2} ')
zlabel('x_3 ')
view(-25,65)

subplot(2,2,4)
solution.xopt.fplot('obj','showgrid',1==0)
xlabel('\theta_{1} ')
ylabel('\theta_{2} ')
zlabel('z\bf (\theta ) ')
view(-25,65)

% figure
% solution.xopt.fplot('primal','position',1)
% xlabel('\theta_{1} ')
% ylabel('\theta_{2} ')
% zlabel('x_1 ')


