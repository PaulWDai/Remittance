% Life Cycle Remittance of Rural to Urban Migrants
% Author: Weifeng Dai
% Date: Sep 20
% Revised: no land endowment, and adding learning mechanism

clear all;
clc;
close all;

%% Parametes
% Calibrated Paremeters using fminsearch

baseline = [0.0628,0.3138,0.0933,0.0102,-0.9804,0.1569,1.9671,0.7886];
% baseline = [0.0893,1.0585,0.0626,0.0071,-0.3987,0.0014,4.4480,0.0543];
delta = baseline(1);
var_e = baseline(2);
alpha = baseline(3);
psi = baseline(4);
f = baseline(5);
theta = baseline(6);
var_k = baseline(7);
zeta = baseline(8);

% Pre-assigned Parameter

beta = .97;             % discount factor
gamma = 2;              % RRA parameter
wR = exp(8.85);         % rural wage
wRf = exp(8.85);        % 
wU = 0;                 % urban wage (all captured in the age-income profile)
r = .03;                % net return
eta0 = 8.30;            % intercept of migrant wage profile
eta1 = .098;            % first order parameter of migrant wage profile
eta2 = -0.00125;        % second order parameter of migrant wage profile


% Calibrated parameter
% zeta = 1;             % experience cost
% nu = 3;               % experience elasticity
%theta = .15;            % exogenous exit from urban
%delta = .06;            % death rate of parents per period
%alpha = .10;            % altruism weight
%psi = .01;              % altruism elasticity
%f = -1;                 % migration fixed cost
mu = 0;                 % land yielding  
%zeta = 0.8;             % land associated cost
lambda = 0.1;           % parameter of experience transitional matrix
c0 = 0;                 % non-homothetic threshold
c0f = 0;                % non-homothetic threshold

mean_k = 8.85;          % mean of initial capital value
%var_k = 2;              % variance of initial capital value
mean_e = 0;             % mean of transitory shock
%var_e = 0.3;            % variance of transitroy shock

std_y25 = .485;         % income variance at age 25
std_y26 = .543;         % income variance at age 26

% graph setting
set(0,'defaultTextInterpreter','latex');
% seed setting
rng('default');

%% Data and Emipirical Analysis (Import from csv Files)
%%% Life Cycle:
% The first column : age
% The second column: mean_log_remit (from data)
% The third column: mean_log_net_remit (from data)
% The fourth column: mean_no_parent (from data)
data_life_cycle = csvread('life_cycle_remit_parent.csv',1,0);

%%% Rural Age Distribution
% The first column : age
% The second column: no_rural in the sample (from data)
% The third column: proportion (from data)

data_age_rural = csvread('age_distribution_rural.csv',1,0);
no_rural = sum(data_age_rural(:,2));
data_age_rural(:,3) = data_age_rural(:,2)/no_rural;

%%% Migrant Age Distribution
% The first column : age
% The second column: no_migrant in the sample (from data)
% The third column: proportion (from data)

data_age_migrant = csvread('age_distribution_migrant.csv',1,0);
no_migrant = sum(data_age_migrant(:,2));
data_age_migrant(:,3) = data_age_migrant(:,2)/no_migrant;

%%% Remittance Distribution for Migrant
% The first column: remittance for group no_parent == 1
% The second column: remittance for group no_parent == 2
data_log_remit = csvread('log_remit_distribution.csv',1,0);
data_log_remit(data_log_remit==0)=nan;

%%% Consumption Profile for Migrant
% The first column: age
% The second column: wage profile
% The third column: saving profile
% The fourth column: consumption profile
data_log_consum = csvread('consum_profile.csv',1,0);

%% State Variables

% Number of parents
nmin = 10^-5;
nmax = 2;
% Age
amin = 24;
amax = 56; % Set a larger amax to avoid the inaccuracy from max(a+1,indexa)
% Land
Lmin = 0;
Lmax = 3.6;
% Transitory shock
emin = -2;
emax = 2;
% Capital
kmin = 0.0001;
kmax = 20;
% Grid Setting
indexn = 3;
indexa = 32;
indexL = 10;
indexe = 10;
indexk = 20;

ngrid = linspace(nmin,nmax,indexn);
agrid = linspace(amin,amax,indexa);
Lgrid = linspace(Lmin,Lmax,indexL);
egrid = linspace(emin,emax,indexe);
kgrid = linspace(kmin,kmax,indexk);

[n5D,a5D,L5D,e5D,k5D] =  ndgrid(ngrid,agrid,Lgrid,egrid,kgrid);

%% Transitional Matrix

indexstate = indexn*indexa*indexL*indexe;
tic

%%% Matrix for n %%%  
         % n = 0    n = 1       n = 2    
Trans_n = [1,       0,          0;...               % n = 0
           delta,   1-delta,    0;...               % n = 1
           delta^2, delta,      1-delta-delta^2];   % n = 2
       
Pn = repmat(Trans_n,1,1,indexa,indexL,indexe,...
                        indexa,indexL,indexe);
Pn = permute(Pn,[1,3,4,5,2,6,7,8]);
Pn = reshape(Pn,[indexstate,indexstate]);

%%% Matrix for e %%%

Pe = normcdf(egrid,mean_e,sqrt(var_e));
Pe = Pe/sum(Pe);
Pe = repmat(Pe,indexe,1);
Pe = repmat(Pe,1,1,indexn,indexa,indexL,...
                   indexn,indexa,indexL);
Pe = permute(Pe,[3,4,5,1,6,7,8,2]);
Pe = reshape(Pe,[indexstate,indexstate]);

%%% Matrix for a %%%

Pa = zeros(indexa,indexa);
for ia = 1:indexa
    Pa(ia,min(ia+1,indexa))=1;
end
Pa = repmat(Pa,1,1,indexn,indexL,indexe,...
                   indexn,indexL,indexe);
Pa = permute(Pa,[3,1,4,5,6,2,7,8]);
Pa = reshape(Pa,[indexstate,indexstate]);

%%% Matrix for L %%%

PL = diag(ones(1,indexL));
PL = repmat(PL,1,1,indexn,indexa,indexe,...
                   indexn,indexa,indexe);
PL = permute(PL,[3,4,1,5,6,7,2,8]);
PL = reshape(PL,[indexstate,indexstate]);

P = sparse(Pn.*Pa.*PL.*Pe);

toc

disp('------------------------------------------');
disp('--------Transitional Matrix: Done---------');
disp('------------------------------------------');

%% Instantaneous Utility Function
% All the k related matrix are taken exponoentially

n2D = reshape(n5D,[indexstate,indexk]);
a2D = reshape(a5D,[indexstate,indexk]);
L2D = reshape(L5D,[indexstate,indexk]);
e2D = reshape(e5D,[indexstate,indexk]);
k2D = exp(reshape(k5D,[indexstate,indexk]));

% Income profile
y2D = exp(eta0 + eta1*a2D + eta2*a2D.^2 + e2D); % take the exp afterwards 

% 3-D grid: state, k, kLead
n3D = repmat(n2D,1,1,indexk);
A3D = alpha*n3D.^(-psi);
k3D = repmat(k2D,1,1,indexk);
L3D = repmat(L2D,1,1,indexk);
e3D = repmat(e2D,1,1,indexk);
y3D = repmat(y2D,1,1,indexk);
kLead3D = permute(k3D,[1 3 2]);

% Agents staying in the Rural sector
cR3D = ... % consumption of rural agent (potential migrant)
    wR + (1+r)* k3D - kLead3D - c0;
cRf3D = wRf + mu * (exp(L3D)-1) - c0f;
%cRf3D = wRf - c0f;
cR3D(cR3D<=0) = 10^-5;
cRf3D(cRf3D<=0) = 10^-5;
uR3D = ... % utility of rural agent
    cR3D.^(1-gamma)/(1-gamma) ... 
    + A3D.*n3D.*cRf3D.^(1-gamma)/(1-gamma);

% Agents staying in the Urban sector
T3D = ... % transfer of urban agent to family member (total)
    (y2D + (1+r)*k3D -kLead3D - A3D.^(-1/gamma).*(wRf -c0f) -c0)./...
    (1+1./n3D.*A3D.^(-1/gamma)); % FOC, see appendix
T3D(T3D<=0) = 10^-5;
cU3D = ... % consumption of urban agent
     y2D + (1+r)*k3D - kLead3D - T3D - c0;
cU3D(cU3D<=0) = 10^-5;
cUf3D = T3D./n3D + wRf -c0f;
uU3D = ... % utility for urban agent
    cU3D.^(1-gamma)/(1-gamma) ...
    + A3D.*n3D.*cUf3D.^(1-gamma)/(1-gamma);

%% Value Function
% Set up
err = 1;
tol = 10^(-4);
maxism = 1000;
iter = 0;
VU2D = ones(indexstate,indexk);
VR2D = ones(indexstate,indexk);
V2D = ones(indexstate,indexk);

tic
while err>tol && iter<maxism
    iter = iter+1;
    V2Dlast = V2D;
    VR2Dlast = VR2D;
    VU2Dlast = VU2D;
    EV2D = P * V2D;
    %EV2D = V2D;
    %{
    VU2D = ... % urban agent value
        uU3D - (QU3D+f) + beta* permute(repmat(EV2D,1,1,indexk),[1 3 2]);
    %}
    VU2D = ... % urban agent value
        uU3D - f - zeta * (exp(L3D)-1) ...
        + beta* (1-theta)*permute(repmat(EV2D,1,1,indexk),[1 3 2])...
        + beta* theta * permute(repmat(VR2D,1,1,indexk),[1 3 2]);
    VR2D = ... % rural agent value
        uR3D + beta* permute(repmat(EV2D,1,1,indexk),[1 3 2]);
    [VU2D,IU] = max(VU2D,[],3); % IU for urban, document location of optK grid
    [VR2D,IR] = max(VR2D,[],3);
    V2D = max(VU2D,VR2D);
    err = max(abs(V2D-V2Dlast),[],'all');
end
toc
disp('------------------------------------------')
if err<tol && iter< maxism
    disp('Iteration is finished successfully');
else
    disp('Error on iteration!');
end
disp(['Error:               ' num2str(err)]);
disp(['Number of iteration: ' num2str(iter)]);
disp('------------------------------------------')

Urban = zeros(indexstate,indexk);
Vdiff = VU2D-VR2D;
Urban(Vdiff>0) = 1;
optk2D_U = kgrid(IU);
optk2D_R = kgrid(IR);

ksolve2D_U = optk2D_U;
ksolve2D_R = optk2D_R;

ksolve2D_U(Urban==0) = nan; % rule out 'not urban'
ksolve2D_R(Urban==1) = nan; % rule out 'not rural'

% Without the consideration of migration decision
ksolve6Duu = reshape(optk2D_U,[indexn,indexa,indexL,indexe,indexk]);
ksolve6Drr = reshape(optk2D_R,[indexn,indexa,indexL,indexe,indexk]);
% With the consideration of migration decision
ksolve6D_U = reshape(ksolve2D_U,[indexn,indexa,indexL,indexe,indexk]);
ksolve6D_R = reshape(ksolve2D_R,[indexn,indexa,indexL,indexe,indexk]);

%% Simulation: Parameter and Initial Distribution
no = 10000;             % no. of observation in the simulated sample
a0 = amin;              % age starts from amin
mean_L = 1.17;          % mean of log(1+LandSize)
std_L = 0.75;           % variance of log(1+LandSize)
percent_land0 = .14;

mean_Lpositive = 1.337;
std_Lpositive = 0.647;

% L0 = normrnd(mean_L,std_L,[1,no]);
L0 = zeros(1,no);
no_Lzero = round(no*percent_land0);
L0(1:no_Lzero) = 0;
L0(no_Lzero:no) = normrnd(mean_Lpositive,std_Lpositive,[1,length(no_Lzero:no)]);
L0(L0<0) = 0;
k0 = normrnd(mean_k,sqrt(var_k),[1,no]);

k0ind = round((k0-kmin)/(kmax-kmin)*(indexk-1))+1;

k_ind(1,:) = k0ind;
k_ind(k_ind<1)=1;
k_ind(k_ind>indexk)=indexk;
kU_ind(1,:) = k_ind(1,:);
kR_ind(1,:) = k_ind(1,:);

k_simu = zeros(indexa,no);
k_simu(1,:) = k0;
kU_simu = k_simu;
kR_simu = k_simu;

% Simulation of n
n_ind = simulate(dtmc(Trans_n),indexa-1,'X0',[0 0 no]);
n_simu = ngrid(n_ind);
% Simulation of a
a_ind = repmat((1:indexa)',1,no);
a_simu = repmat(agrid',1,no);

L_simu = repmat(L0,indexa,1);
L_ind = round((L_simu-Lmin)/(Lmax-Lmin)*(indexL-1))+1;
L_ind(L_ind<1) = 1; % adjust the index to include it in the original grid
L_ind(L_ind>indexL) = indexL;

% Simulation of e
e_simu = normrnd(mean_e,sqrt(var_e),[indexa,no]);
e_ind = round((e_simu-emin)/(emax-emin)*(indexe-1))+1;
e_ind(e_ind<1) = 1;
e_ind(e_ind>indexe) = indexe;

%{
% Simulation of m
m0_simu = normrnd(mean_x,sqrt(var_x),[1,no]);
m_simu = zeros(indexa,no);
m_simu(1,:) = m0_simu;
for i = 1:indexa-1
    m_simu(i+1,:) = (vgrid(i)* m_simu(i,:) + var_e*e_simu(i,:))...
        /(vgrid(i)+var_e);
end
m_ind = round((m_simu-mmin)/(mmax-mmin)*(indexm-1))+1;
m_ind(m_ind<1) = 1;
m_ind(m_ind>indexm) = indexm;
%}
%{
y_simu = repmat((eta0+eta1*agrid+eta2*agrid.^2)',1,no);
y_ind = round((y_simu-ymin)/(ymax-ymin)*(indexy-1))+1;
y_ind(y_ind<1) = 1; % adjust the index to include it in the original grid
y_ind(y_ind>indexy) = indexy;
%}
%% Simulation: Life Cycle of Asset Allocation of Heterogeneous Agents
% Life Cycle of Urban
indU = zeros(indexa,no);
for i = 1:indexa-1
    indU(i,:) = sub2ind(size(ksolve6Duu),...
        n_ind(i,:),a_ind(i,:),L_ind(i,:),e_ind(i,:),kU_ind(i,:));
    indTempU = indU(i,:);
    kU_simu(i+1,:) = ksolve6Duu(indTempU); % without consideration of decision
    kU_indTemp = round((kU_simu(i+1,:)-kmin)/(kmax-kmin)*(indexk-1))+1;
    kU_indTemp(kU_indTemp<1) = 1;
    kU_indTemp(kU_indTemp>indexk) = indexk;
    kU_ind(i+1,:) = kU_indTemp;
end

% Life Cycle of Rural
indR = zeros(indexa,no);
for i = 1:indexa-1
    indR(i,:) = sub2ind(size(ksolve6Drr),...
        n_ind(i,:),a_ind(i,:),L_ind(i,:),e_ind(i,:),kR_ind(i,:));
    indTempR = indR(i,:);
    kR_simu(i+1,:) = ksolve6Drr(indTempR); % without consideration of decision
    kR_indTemp = round((kR_simu(i+1,:)-kmin)/(kmax-kmin)*(indexk-1))+1;
    kR_indTemp(kR_indTemp<1) = 1;
    kR_indTemp(kR_indTemp>indexk) = indexk;
    kR_ind(i+1,:) = kR_indTemp;
end

% Life cycle of decision
Migr_simu = Urban(indU(1:indexa-1,:)); % from amin to amax-1
MigrRate = sum(Migr_simu,2)/no*(1-theta); % considering exogenous exit from urban

% Truncated age profile due to the expression of remittance T
kLeadU_simu = kU_simu(2:indexa,:);
kU_simu = kU_simu(1:indexa-1,:);
kR_simu = kR_simu(1:indexa-1,:);
n_simu = n_simu(1:indexa-1,:);
A_simu = alpha*n_simu.^(-psi);
y_simu = eta0+eta1*a_simu+eta2*a_simu.^2+e_simu;
y_simu = y_simu(1:indexa-1,:);
T_simu = ... % transfer of urban agent to family member (total) simulated
    (exp(y_simu) + (1+r)*exp(kU_simu) -exp(kLeadU_simu) - c0 -...
    A_simu.^(-1/gamma).*(wRf-c0f))./...
    (1+1./n_simu.*A_simu.^(-1/gamma)); % FOC, see appendix
T_M_simu = T_simu;
T_M_simu(Migr_simu==0) = nan;
T_M_simu(T_M_simu<0) = 0;
Mean_T_M = nanmean(T_M_simu,2);
Mean_log_T_M = nanmean(log(1+T_M_simu),2);



%% Calibration: Targeted Moments

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Migrants share %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% Migr_share = sum(Migr_simu,'all')/((indexa-1)*no);
Migr_share_age = mean(Migr_simu,2);
Migr_share = sum(Migr_share_age.* data_age_migrant(:,3))*(1-theta);

Migr_share_non_adjusted = nanmean(Migr_simu,'all')*(1-theta);
Data_Migr_2005 = 132.6*16.5/17.8 + 47*5.2/6.5; % unit: million
Data_Migr_2010 = 176.2*21.6/22.9+79.2*8.6/10.3;
Data_Migr_share = (Data_Migr_2005+Data_Migr_2010)/2/(714);
disp(['The share of migration, Data =' num2str(Data_Migr_share*100),'%, Model = ',num2str(Migr_share*100),'%'] );
disp(['The share of migration, Data =' num2str(Data_Migr_share*100),'%, Model = ',num2str(Migr_share_non_adjusted*100),'%'] );

disp('%%%%%%%%%%%%%%%%%% Remittance:over group of parents %%%%%%%%%%%%%%%%%%%%%%%')
T_M_simu_Parent1 = T_M_simu;
T_M_simu_Parent2 = T_M_simu;
T_M_simu_Parent0(n_simu>=1) = nan;
T_M_simu_Parent1(n_simu<0.99|n_simu>1.01) = nan;
T_M_simu_Parent2(n_simu<1.01) = nan;
log_T_M_simu_Parent1 = reshape(log(1+T_M_simu_Parent1),[1 no*(indexa-1)]);
log_T_M_simu_Parent2 = reshape(log(1+T_M_simu_Parent2),[1 no*(indexa-1)]);

mean_log_T_M_Parent1 = nanmean(log_T_M_simu_Parent1,'all');
mean_log_T_M_Parent2 = nanmean(log_T_M_simu_Parent2,'all');
disp(['Mean remittance (Parent = 1), Data = 5.02 ,Model = ',num2str(mean_log_T_M_Parent1)] );
disp(['Mean remittance (Parent = 2), Data = 5.10 ,Model = ',num2str(mean_log_T_M_Parent2)] );
Data_Mean_remittance_Parent_1 = 5.02;
Data_Mean_remittance_Parent_2 = 5.10;


%%%%%%%%%%%%%%%% Remittance: with or without land %%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%% Remittance over Land Properties %%%%%%%%%%%%%%%%%%%%%%')
%T_M_simu_Parent0 = T_M_simu;
T_M_simu_Land0 = T_M_simu;
T_M_simu_Land1 = T_M_simu;
T_M_simu_0(n_simu>=1) = nan;

disp('%%%%%%%%%%%%%%%%%%%%%%% Average Duration of Migration %%%%%%%%%%%%%%%%%%%%%%%%%');
YearMigr = sum(Migr_simu);
mean_YearMigr = mean(YearMigr);
Data_Migr_2005 = 132.6*16.5/17.8 + 47*5.2/6.5; % unit: million
Data_Migr_2010 = 176.2*21.6/22.9+79.2*8.6/10.3;
Data_Migr_share = (Data_Migr_2005+Data_Migr_2010)/2/(703.99);
disp(['Average Duration migration, Data = 7 ,Model = ',num2str(mean_YearMigr)]);

disp('%%%%%%%%%%%%%%%%%%% Distribution of Land Properties %%%%%%%%%%%%%%%%%%%%%%%');
L_Migr = L_simu;
L_Migr(Migr_simu ==0) =nan;
[len_L_Migr,wid_L_Migr] = size(L_Migr);
L_Migr = reshape(L_Migr,[1 len_L_Migr*wid_L_Migr]);
mean_L_Migr = nanmean(L_Migr,'all');
std_L_Migr = nanstd(L_Migr);
disp(['Mean land, Data =  1, Model = ',num2str(mean_L_Migr)]);
disp(['Standard deviation Land, Data = 0.45, Model =', num2str(std_L_Migr)]);

disp('%%%%%%%%%%%%%%%%% Number of Parents by Migrants Age %%%%%%%%%%%%%%%%%%%%%%');
n_simu_Migr = n_simu;
n_simu_Migr(Migr_simu==0) = nan;
mean_n_simu = nanmean(n_simu_Migr,2);
%mean_n_simu = mean_n_simu(1:indexa);
% To calibrate the death rate of delta
delta_loss = mean(abs(mean_n_simu-data_life_cycle(:,3)));
disp(['The deviation for delta',num2str(delta_loss)]);

disp('%%%%%%%%%%%%%%%%% Variance of simulated income %%%%%%%%%%%%%%%%%%%%%%');
y_simu_25 = y_simu(2,:);
y_simu_26 = y_simu(3,:);
std_y_simu_25 = nanstd((y_simu_25),0,'all');
std_y_simu_26 = nanstd((y_simu_26),0,'all');
disp(['Std income age 25, Data = ',num2str(std_y25), 'Model = ',num2str(std_y_simu_25)]);
disp(['Std income age 26, Data = ',num2str(std_y26), 'Model = ',num2str(std_y_simu_26)]);

disp('%%%%%%%%% Mean and variance of Remit, Consumption, Saving %%%%%%%%%%%%%%%')
cU_simu = exp(y_simu) + (1+r)*exp(kU_simu) - exp(kLeadU_simu) - T_simu;
cU_simu(cU_simu<0) = 10^-5;
log_cU_simu = log(cU_simu+1);
log_cU_simu(Migr_simu ==0) = nan;
mean_log_cU_simu= nanmean(log_cU_simu,'all');
disp(['Mean log(1+consumption), Data = 9.48 ,Model = ',num2str(mean_log_cU_simu)]);
Data_Mean_log_consumption = 9.48;

mean_log_T_M_simu = nanmean(log(T_M_simu+1),'all');
disp(['Mean log(1+Remit), Data = 4.73 ,Model = ',num2str(mean_log_T_M_simu)]);
Data_mean_log_remit = 4.73;

savingU_simu_Migr = exp(kLeadU_simu);
savingU_simu_Migr(Migr_simu==0) = nan;
mean_log_savingU_simu = nanmean(log(1+savingU_simu_Migr),'all');
disp(['Mean log(1+Saving), Data =  1.46 ,Model = ',num2str(mean_log_savingU_simu)]);
Data_mean_log_saving = 1.46;

std_log_T_M_simu = nanstd(log(T_M_simu+1),0,'all');
std_log_savingU_simu = nanstd(log(1+savingU_simu_Migr),0,'all');
std_log_cU_simu= nanstd(log_cU_simu,0,'all');

disp(['std of log(1+Remit), Data = 0.70, Model =', num2str(std_log_T_M_simu)]);
disp(['std of log(1+saving), Data = 3.15, Model =', num2str(std_log_savingU_simu)]);
disp(['std of log(1+Consumption), Data = 0.70, Model = ', num2str(std_log_cU_simu)]);
Data_std_log_consumption = 0.70;

[mean_n_simu,data_life_cycle(:,3)]
%% Figures
% Migration rate over ages

%{
figure(1);
subplot(3,2,1);
plot(agrid(2:indexa),MigrRate,...
    'LineWidth',2);
xlabel('Age');ylabel('Migration Rate');
title('Migration Rate');ylim([0,1]);
% Remittance pattern over ages
subplot(3,2,2);
scatter(agrid(2:indexa),data_life_cycle(:,2));hold on;
plot(agrid(2:indexa-1),Mean_log_T_M(2:indexa-1),...
    'LineWidth',2);
xlabel('Age');ylabel('Mean of $\log(1+Remittance)$');
% Remittance by parents group
subplot(3,2,3);
histogram(log_T_M_simu_Parent1,'Normalization','probability','BinWidth',1);
xlabel('$\log(1+T^M)$');ylabel('Density');
title('1 Parent')
subplot(3,2,4);
histogram(log_T_M_simu_Parent2,'Normalization','probability','BinWidth',1);
xlabel('$\log(1+T^M)$');ylabel('Density');
title('2 Parents')

log_T_M_simu_Parent1_nonzero = log_T_M_simu_Parent1;
log_T_M_simu_Parent1_nonzero(log_T_M_simu_Parent1_nonzero==0) = nan;
log_T_M_simu_Parent2_nonzero = log_T_M_simu_Parent2;
log_T_M_simu_Parent2_nonzero(log_T_M_simu_Parent2_nonzero==0) = nan;


subplot(3,2,5);
histogram(log_T_M_simu_Parent1_nonzero,...
    'Normalization','probability','BinWidth',1,'DisplayStyle','stairs');hold on;
histogram(data_log_remit(:,1),...
    'Normalization','probability','BinWidth',1,'DisplayStyle','stairs');

xlabel('$\log(1+T^M)$');ylabel('Density');
title('Non-Negative Remittance(1 Parent)')

subplot(3,2,6);
histogram(log_T_M_simu_Parent2_nonzero,...
    'Normalization','probability','BinWidth',1,'DisplayStyle','stairs');hold on;
histogram(data_log_remit(:,2),...
    'Normalization','probability','BinWidth',1,'DisplayStyle','stairs');
xlabel('$\log(1+T^M)$');ylabel('Density');
title('Non-negative Remittance(2 Parents)')
%}


figure(2);
scatter(agrid(2:indexa),data_life_cycle(:,2)/data_life_cycle(1,2),'filled');hold on;
plot(agrid(2:indexa-1),Mean_log_T_M(2:indexa-1)/Mean_log_T_M(2),...
    'LineWidth',2);
xlabel('Age','FontSize',16);
ylabel('Mean of $\log(1+Remittance)$','FontSize',16);
xlim([amin amax]);
legend('Data','Model','Interpreter','latex','FontSize',16);
title('Life Cycle Remittance: Data and Model','FontSize',16);

% Age distribution of Migrant
%{
figure(2)
no_migrant_age = data_age_rural(:,3).*MigrRate;
density_migrant_age = no_migrant_age/sum(no_migrant_age);
plot(agrid(2:indexa),density_migrant_age),hold on;
scatter(agrid(2:indexa),data_age_migrant(:,3));
%}
% Check for policy function

% median of grid
%{
figure(2);
medL = round(.5*(1+indexL));
mede = round(.5*(1+indexe));

subplot(1,3,1);
[AA_ak,KK_ak]= ndgrid(agrid,kgrid);
surf(AA_ak,KK_ak,permute(ksolve6D_U(1,:,medL,mede,:),[2,5,1,3,4]));view(0,90);
xlim([amin amax]);ylim([kmin kmax]);colorbar

subplot(1,3,2);
[AA_ak,KK_ak]= ndgrid(agrid,kgrid);
surf(AA_ak,KK_ak,permute(ksolve6D_U(2,:,medL,mede,:),[2,5,1,3,4]));view(0,90);
xlim([amin amax]);ylim([kmin kmax]);colorbar

subplot(1,3,3);
[AA_ak,KK_ak]= ndgrid(agrid,kgrid);
surf(AA_ak,KK_ak,permute(ksolve6D_U(3,:,medL,mede,:),[2,5,1,3,4]));view(0,90);
xlim([amin amax]);ylim([kmin kmax]);colorbar

Data_moment = [];
%}

age_mean_log_cU = nanmean(log_cU_simu,2);

figure(3);
scatter(amin+1:amax-1,data_log_consum(:,4)/data_log_consum(1,4),'filled'); hold on
plot(agrid(2:indexa-1),age_mean_log_cU(2:indexa-1)/age_mean_log_cU(2),...
    'LineWidth',2);
xlabel('Age','FontSize',16);
ylabel('Mean of $\log(1+Consumption)$','FontSize',16);
xlim([amin amax]);
legend('Data','Model','Interpreter','latex','FontSize',16);
title('Life Cycle Consumption: Data and Model','FontSize',16);
