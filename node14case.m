function mpc = node14case


%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
% Pd Qd in kW
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	1	6.2	1.6	0	0	1	1	0	10	1	1.06	0.94
    2	1	8.2	2.5	0	0	1	1	0	10	1	1.06	0.94
    3	1	3.5	1.8	0	0	1	1	0	10	1	1.06	0.94
    4	2	9	5.8	0	0	1	1	0	10	1	1.06	0.94
    5	2	3.2	0.9	0	0	1	1	0	10	1	1.06	0.94
    6	1	9.5	3.4	0	0	1	1	0	10	1	1.06	0.94
    7	2	2.2	0.7	0	0	1	1	0	10	1	1.06	0.94
    8	3	2.4	1.2	0	0	1	1	0	10	1	1.06	0.94
    9	1	7.6	1.6	0	0	1	1	0	10	1	1.06	0.94
    10	2	5.8	2	0	0	1	1	0	10	1	1.06	0.94
    11	1	3.2	1.6	0	0	1	1	0	10	1	1.06	0.94
    12	1	8.7	6.7	0	0	1	1	0	10	1	1.06	0.94
    13	1	3.5	2.3	0	0	1	1	0	10	1	1.06	0.94
    14	1	2.4	1.2	0	0	1	1	0	10	1	1.06	0.94

];

%% generator data
% (KW)
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	4	0	0	30	-20	1	100	1	90	24	0	0	0	0	0	0	0	0	0	0	0
    5	0	0	50	-20	1	100	1	120	0	0	0	0	0	0	0	0	0	0	0	0
    7	0	0	50	-20	1	100	1	70	27	0	0	0	0	0	0	0	0	0	0	0
    8	0	0	70	-20	1	100	1	110	40	0	0	0	0	0	0	0	0	0	0	0
    10	0	0	70	-20	1	100	1	150	0	0	0	0	0	0	0	0	0	0	0	0
];

%% branch data
% r x in ohms
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
1	2	0.2511	0.493	0	0	0	0	0	0	1	-360	360 %1
2	3	0.1864	0.366	0	0	0	0	0	0	1	-360	360 %2
3	4	0.1941	0.3811	0	0	0	0	0	0	1	-360	360 %3
3	5	0.707	0.819	0	0	0	0	0	0	1	-360	360 %4
3	6	0.1872	0.6188	0	0	0	0	0	0	1	-360	360 %5
1	7	0.2351	0.7114	0	0	0	0	0	0	1	-360	360 %6
7	8	0.74	1.03	0	0	0	0	0	0	1	-360	360 %7
7	9	0.74	1.044	0	0	0	0	0	0	1	-360	360 %8
9	10	0.065	0.1966	0	0	0	0	0	0	1	-360	360 %9
1	11	0.1238	0.3744	0	0	0	0	0	0	1	-360	360 %10
11	12	1.155	1.468	0	0	0	0	0	0	1	-360	360 %11
12	13	0.5416	0.7129	0	0	0	0	0	0	1	-360	360 %12
12	14	0.526	0.591	0	0	0	0	0	0	1	-360	360 %13
];

%% convert branch impedances from Ohms to p.u.
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
Vbase = mpc.bus(1, BASE_KV) * 1e3;      %% in Volts
Sbase = mpc.baseMVA * 1e6;              %% in VA
mpc.branch(:, [BR_R BR_X]) = mpc.branch(:, [BR_R BR_X]) / (Vbase^2 / Sbase);

%% convert loads from kW to MW
% mpc.bus(:, [PD, QD]) = mpc.bus(:, [PD, QD]) / 1e3;

%% convert generators from kW to MW
mpc.gen(:,[4,5,9,10]) = mpc.gen(:,[4,5,9,10])/ 1e3;

