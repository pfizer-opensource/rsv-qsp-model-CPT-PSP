function dy = rsv_model(t,x,params,drug_Finhib,drug_Ninhib)
% https://link.springer.com/content/pdf/10.1007/s00285-019-01364-1.pdf

dy=zeros(8,1);

%New model, Upper respiratory tract
%Include new state variables
%% assign variables
% Target Cells (T)
T= max(x(1),0);
% Eclipse cells (E)
E = max(x(2),0);
% Infectious cells (I)
I = max(x(3),0);
% Virus Upper Respiratory (V)
V = max(x(4),0);
% Refractory to infection (R)
R = max(x(5),0);
% Resting innate immune cells (D0)
D0 = max(x(6),0);
% Activated innate immune cells (D1)
D1 = max(x(7),0);
% Damaged bystander cells
J = max(x(8),0);

%% assign parameters

% infectivity rate   /day
beta = params(1);
% eclipse to infectious rate /day
k   = params(2);
% death rate from viral cytopathic effects
delta_I = params(3);
%virion clearance rate /day
c = params(4);
% virion production rate 1 /day
p = params(5);     
% refractory state reversion rate
delta_R = params(6);
% rate of transition to IFN induced refractory state
kappa0 = params(7);
% effect of innate immune response on the clearance of PAMP and DAMP
% expressing cells
kappa1 = params(8);
% effect of adaptive immune response on the clearance of infected cells
kappa2 = params(9);
% maximum bystander damage rate
nu = params(10);
%D1 concentration leading to half the bystander damage rate
K = params(11);
%death rate from cell damage
delta_J = params(12);
%resting immune cells replenishing rate or recruitment rate
lambda = params(13);
%innate immune cell activation rate
sigma_I = params(14);
%activated innate immune cell average death rate
delta_D = params(15);
%Time delay of adaptive immune response post-infection
tau = params(16);
%Hill coefficient for bystander cell damage
m = params(17);
sigma_J = params(18);
sigma_V = params(19);
k_alpha = params(20);

%First dose is 2.5 days post in
t_inoc_delay = 2.5*24;

%Define the delay step function for adaptive immune response
if t < tau + t_inoc_delay
    I_delay = 0;
else
    I_delay = 1;
end

%% sisunatovir PK & Drug Effect (F-protein inhibitor)
mw_Finhib = 446.4;
if drug_Finhib.Dose ~=0
    if drug_Finhib.dose_typ == 1 && drug_Finhib.Dose*mw_Finhib == 200e6
        %200mg sisunatovir mean PK
        %drug_conc_1 = max(PK_check(drug_Finhib,t),0);
        n_totalDoses = length(drug_Finhib.dose_times);
        drug_conc_1 = popPK_profile(t-drug_Finhib.dose_times(1), n_totalDoses,1);
    elseif drug_Finhib.dose_typ == 1 && drug_Finhib.Dose*mw_Finhib == 350e6
        %350mg sisunatovir mean PK
        %drug_conc_1 = max(PK_check(drug_Finhib,t),0);
        n_totalDoses = length(drug_Finhib.dose_times);
        drug_conc_1 = popPK_profile(t-drug_Finhib.dose_times(1), n_totalDoses,2);
    elseif drug_Finhib.dose_typ == 1 && drug_Finhib.Dose*mw_Finhib == 100e6
        %350mg sisunatovir mean PK
        %drug_conc_1 = max(PK_check(drug_Finhib,t),0);
        n_totalDoses = length(drug_Finhib.dose_times);
        drug_conc_1 = popPK_profile(t-drug_Finhib.dose_times(1), n_totalDoses,3);
    elseif drug_Finhib.dose_typ == 2
        %Loading Dose 1 LD 300mg, 200mg - 5 day treatmentsisunatovir mean PK from popPK
        n_LoadingDoses = 1;
        n_totalDoses = length(drug_Finhib.dose_times);
        drug_conc_1 = popPK_DrugConcentration(t-drug_Finhib.dose_times(1),n_LoadingDoses, n_totalDoses);
    elseif drug_Finhib.dose_typ == 3
        %Loading Dose 2 LD 300mg, 200mg - 5 day treatment sisunatovir mean PK from popPK
        n_LoadingDoses = 2;
        n_totalDoses = length(drug_Finhib.dose_times);
        drug_conc_1 = popPK_DrugConcentration(t-drug_Finhib.dose_times(1),n_LoadingDoses, n_totalDoses);
        
    elseif drug_Finhib.dose_typ == 4
        fu = 1;
        drug_conc_1 = fu*max(PK_check(drug_Finhib,t),0);
    
    end
    %Emax model for inhibition
    EC50  = drug_Finhib.EC50;
    Imax  = drug_Finhib.Imax;
    n_hill = drug_Finhib.n;
 
    prod_inhib_Fprotein = 1-(Imax*(drug_conc_1^n_hill))/((EC50)^n_hill+(drug_conc_1)^n_hill);
else
    prod_inhib_Fprotein = 1;
end

%zelicapavir PK & Drug Effect (N-protein inhibitor)
if drug_Ninhib.Dose ~=0
    %Phenomenological N-protein inhibitor
    if drug_Ninhib.dose_typ == 1
        drug_conc_1 = 0.06*max(PK_check(drug_Ninhib,t),0);
    elseif drug_Ninhib.dose_typ == 4
        fu = 1; %0.025;
        drug_conc_1 = fu*max(PK_check(drug_Ninhib,t),0);
    end  
    
    %Emax model for inhibition
    EC50  = drug_Ninhib.EC50;
    Imax  = drug_Ninhib.Imax;
    n_hill = drug_Ninhib.n;
    
    prod_inhib_Nprotein = 1-(Imax*(drug_conc_1^n_hill))/((EC50)^n_hill+(drug_conc_1)^n_hill);
else
    prod_inhib_Nprotein = 1;
end

%% RSV ODEs
%Target cell regeneration
tar_cell_prod = 0*T;
%Anti-viral term for target cell infection
alpha_TV = (k_alpha/(k_alpha+(J + I)));
%Target cell infection
tar_cell_infection         = alpha_TV*beta*T*V*prod_inhib_Fprotein;
%Eclipse cell transition to infectious cells
eclipse_to_infectious      = k*E;
%Clearance of infectious cells
infectious_clearance       = delta_I*I;
innate_immune_infection           = kappa1*D1*I;
adapt_immune_infection = kappa2*I_delay*I;
%Viral shedding
viral_production           = p*I*prod_inhib_Nprotein;
%viral shedding
viral_clearance            = c*V;
refract_on                 = kappa0*I*T;
refract_off                = delta_R*R;        
D0_base                    = lambda*(1 - D0);
%Add in log-sensing to the model
%immune_act                 = sigma_I*I*D0 + sigma_J*J*D0;%sigma*(I+J)*D0;
immune_death               = delta_D*D1;
release_PAMP               = (nu*(D1^m))/(D1^m + K^m);
byst_immune                = kappa1*D1*J;
byst_death                 = delta_J*J;

%log sensing to capture super-producers
immune_act                 = sigma_I*I*D0 + sigma_J*J*D0 + sigma_V*log(1+V)*D0;%sigma*(I+J)*D0;

% RHS of Equations
%% Target cells (T) %%
dy(1) = -tar_cell_infection - refract_on + refract_off + tar_cell_prod;
%% Eclipse cells (E) %%
dy(2) =  tar_cell_infection  - eclipse_to_infectious;
%% Infectious cells (I) %%
dy(3) =  -infectious_clearance  + eclipse_to_infectious - innate_immune_infection - adapt_immune_infection;
%% Virus (V) %%
dy(4) =  viral_production - viral_clearance;
%% Refractory to infection (R) %%
dy(5) = refract_on - refract_off;
%% Resting innate immune cells (D0) %%
dy(6) = D0_base - immune_act;
%% Activated innate immune cells (D1) %%
dy(7) =  immune_act - immune_death;
%% Damaged bystander cells %%
dy(8) = release_PAMP -  byst_death - byst_immune;