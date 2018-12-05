t = 50;
brloss = @(i,j,t) (P0(i,t)^2+Q0(i,t)^2)/V0(j,t)^2 * results.branch(i,BR_R);

brloss(1,1,t)
P0(1,t)
Q0(1,t)