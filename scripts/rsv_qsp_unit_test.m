close all
clear all
t_inoc_delay = 2.5;


%Normalization constants
T0 = 4e8; %cells
D00 = 1e6;

%Load parameters
sis_vc_vpop = load("vpop_sisunatovir_VC.mat");
VL_peak_PBO = sis_vc_vpop.vpop_peak_pbo;

load best_fit_Finhib.mat % loads best fit param vector as p_fit_new

% Initial parameter set
V0_PBO = p_est_new(end-2); %inoc amount
p_est_new = p_est_new([1:18 20 21]); %remove inoc from param vec.

%Initialize state variables
ICs = zeros(10,1);
%New model, Upper respiratory tract
%Include new state variables
%% assign variables
% Target Cells (T)
IC(1)= 1;
% Eclipse cells (E)
IC(2) = 0;
% Infectious cells (I)
IC(3) = 0;
% Virus Upper Respiratory (V)
IC(4) = V0_PBO;
% Refractory to infection (R)
IC(5) = 0;
% Resting innate immune cells (D0)
IC(6) = 1;
% Activated innate immune cells (D1)
IC(7) = 0;
% Damaged bystander cells
IC(8) = 0/T0;

%Initial conditions
IC_PBO_new = IC;

% Assign dose variables for ODE struct
% Other than being defined, does not matter as dose = 0 for unit tests
dose_type = 1;
%% No treatment - Dose = 0
drug.Dose = 0; %nmol;
drug.k_a = 0.5;  % /day
drug.t_half = 8.54; % day
drug.Vd = 1e4;     % L
drug.dose_times = (t_inoc_delay + [0.5:0.5:5])*24; %hours
mw_Finhib = 446.4;% g/mol
Finhib = 200e-3; %g
drug.Imax = 1;
drug.EC50 =  0.4; %nmol/L
drug.n = 1;
drug.dose_typ = dose_type;

%% UNIT TEST 1 - Nominal case of healthy response.
odefun = @(t,x)rsv_model(t,x,p_est_new,drug,drug);
OPTIONS = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t_U1,y_U1] = ode15s(odefun,linspace(0,(7+t_inoc_delay+40.5)*24,289),...
    IC,OPTIONS);

%% Unit test 2 -  No innate immune recognition of virus.
p_est_new_U2 = p_est_new;
%Unit test 2: No innate immune recognition of virus.
p_est_new_U2(7) = 0; %kappa_0
p_est_new_U2(9) = 0; %kappa_2
p_est_new_U2(14) = 0; %sigma_I
p_est_new_U2(18) = 0; %sigma_J
p_est_new_U2(19) = 0; %sigma_V

odefun = @(t,x)rsv_model(t,x,p_est_new_U2,drug,drug);
OPTIONS = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t_U2,y_U2] = ode15s(odefun,linspace(0,(10+t_inoc_delay+40.5)*24,289),...
    IC,OPTIONS);

%% Unit Test 3 - Sterile inflammatory response
IC_U3 = IC;
IC_U3(4) = 0;
IC_U3(8) = 0.1;
odefun = @(t,x)rsv_model(t,x,p_est_new,drug,drug);
OPTIONS = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t_U3,y_U3] = ode15s(odefun,linspace(0,(10+t_inoc_delay+40.5)*24,289),...
    IC_U3,OPTIONS);

%% Unit Test 4 - Sustained inflammatory response
IC_U4 = IC;
IC_U4(8) = 0;
p_est_new_U4 = p_est_new;
%Unit test 2: No innate immune recognition of virus.
p_est_new_U4(10) = 1.50*p_est_new_U4(10); %nu
p_est_new_U4(19) = 1.50*p_est_new_U4(19); %sigma_V

odefun = @(t,x)rsv_model(t,x,p_est_new_U4,drug,drug);
OPTIONS = odeset('RelTol',1e-3,'AbsTol',1e-6);
[t_U4,y_U4] = ode15s(odefun,linspace(0,(10+t_inoc_delay+40.5)*24,289),...
    IC_U4,OPTIONS);
%% Plots of all the unit tests


%% Combination figure
f = figure;
f.Position = [100 100 1140 600];
tlo = tiledlayout(2,4,'TileSpacing','Compact');
nexttile
plot(t_U1/24,y_U1(:,1),'LineWidth',3)
hold on
plot(t_U2/24,y_U2(:,1),'LineWidth',3)
hold on
plot(t_U4/24,y_U4(:,1),'LineWidth',3)
hold off
grid on
xlim([0 46.5+t_inoc_delay+0.5])
xlabel('Days since Infection')
ylabel({'Susceptible Cells (S)';'Normalized to initial baseline'})
xticks(2.5:5:49.5)
xticklabels({'0','5','10','15','20','25','30','35','40','45'})
grid on
set(gca,'Fontsize',12)

nexttile
plot(t_U1/24,y_U1(:,5),'LineWidth',3)
hold on
plot(t_U2/24,y_U2(:,5),'LineWidth',3)
hold on
plot(t_U4/24,y_U4(:,5),'LineWidth',3)
hold off
grid on
xlim([0 46.5+t_inoc_delay+0.5])
ylabel({'Refractory Cells (R)';'normalized to initial baseline'})
xlabel('Days since Infection')
xticks(2.5:5:49.5)
xticklabels({'0','5','10','15','20','25','30','35','40','45'})
grid on
set(gca,'Fontsize',12)

nexttile
semilogy(t_U1/24,y_U1(:,2)/max(y_U1(:,2)),'LineWidth',3)
hold on
semilogy(t_U2/24,y_U2(:,2)/max(y_U1(:,2)),'LineWidth',3)
hold on
semilogy(t_U4/24,y_U4(:,2)/max(y_U1(:,2)),'LineWidth',3)
hold off
grid on
xlim([0 46.5+t_inoc_delay+0.5])
ylim([1e-2 1e2])
ylabel({'Eclipsed Cells (E)';'Normalized to peak nominal response'})
xlabel('Days since Infection')
xticks(2.5:5:49.5)
xticklabels({'0','5','10','15','20','25','30','35','40','45'})
grid on
set(gca,'Fontsize',12)

nexttile
semilogy(t_U1/24,y_U1(:,3)/max(y_U1(:,3)),'LineWidth',3)
hold on
semilogy(t_U2/24,y_U2(:,3)/max(y_U1(:,3)),'LineWidth',3)
hold on
semilogy(t_U4/24,y_U4(:,3)/max(y_U1(:,3)),'LineWidth',3)
hold on
hold off
grid on
xlim([0 46.5+t_inoc_delay+0.5])
ylim([1e-2 1e3])
xlabel('Days since first dose')
ylabel({'Infectious Cells (I)';' Normalized to peak nominal response'})
xlabel('Days since Infection')
xticks(2.5:5:49.5)
xticklabels({'0','5','10','15','20','25','30','35','40','45'})
grid on
set(gca,'Fontsize',12)

nexttile
semilogy(t_U1/24,y_U1(:,4)*VL_peak_PBO ,'LineWidth',3)
hold on
semilogy(t_U2/24,y_U2(:,4)*VL_peak_PBO ,'LineWidth',3)
hold on
semilogy(t_U4/24,y_U4(:,4)*VL_peak_PBO ,'LineWidth',3)
hold off
grid on
xlim([0 46.5+t_inoc_delay+0.5])
ylim([1e0 1e9])
xlabel('Days since Infection')
ylabel('Virus Titer (V) (log10 PFUe/mL)')
xticks(2.5:5:49.5)
xticklabels({'0','5','10','15','20','25','30','35','40','45'})
grid on
set(gca,'Fontsize',12)

nexttile
p1 = plot(t_U1/24,y_U1(:,6)/max(y_U1(:,6)),'LineWidth',3);
hold on
p2 = plot(t_U2/24,y_U2(:,6)/max(y_U1(:,6)),'LineWidth',3);
hold on
p3 = plot(t_U4/24,y_U4(:,6)/max(y_U1(:,6)),'LineWidth',3);
hold off
grid on
xlim([0 46.5+t_inoc_delay+0.5])
ylabel({'Unactivated Immune Cells (D0)'; 'Normalized to peak nominal response'})
xlabel('Days since Infection')
xticks(2.5:5:49.5)
xticklabels({'0','5','10','15','20','25','30','35','40','45'})
grid on
set(gca,'Fontsize',12)

nexttile
plot(t_U1/24,y_U1(:,7)/max(y_U1(:,7)),'LineWidth',3);
hold on
plot(t_U2/24,y_U2(:,7)/max(y_U1(:,7)),'LineWidth',3);
hold on
plot(t_U4/24,y_U4(:,7)/max(y_U1(:,7)),'LineWidth',3);
hold off
grid on
xlim([0 46.5+t_inoc_delay+0.5])
ylabel({'Activated Immune Cells (D1)'; ' Normalized to peak normial response'})
xlabel('Days since Infection')
xticks(2.5:5:49.5)
xticklabels({'0','5','10','15','20','25','30','35','40','45'})
grid on
set(gca,'Fontsize',12)

nexttile
plot(t_U1/24,y_U1(:,8)/max(y_U1(:,8)),'LineWidth',3);
hold on
plot(t_U2/24,y_U2(:,8)/max(y_U1(:,8)),'LineWidth',3);
hold on
plot(t_U4/24,y_U4(:,8)/max(y_U1(:,8)),'LineWidth',3);
hold off
grid on
xlim([0 46.5+t_inoc_delay+0.5])
ylabel({'Damaged Bystander Cells (J)'; ' Normalized to peak nominal response'})
xlabel('Days since Infection')
xticks(2.5:5:49.5)
xticklabels({'0','5','10','15','20','25','30','35','40','45'})
grid on
set(gca,'Fontsize',12)
T = legend([p1(1),p2(1),p3(1)],'Nominal Reponse - U1','No Immune Response - U2',...
    'Sustained Inflammation - U3','Orientation','horizontal');
set(T,'location','southoutside')


print(gcf,'figs/supp_unit_test.png','-dpng','-r300');