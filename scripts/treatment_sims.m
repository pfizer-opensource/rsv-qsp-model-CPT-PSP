load("vpop_EDP938_Peds.mat")
edp938_vc_vpop = load("vpop_EDP938_pt1_VC.mat");
edp938_vc_vpop2 = load("vpop_EDP938_pt2_VC.mat");
sis_vc_vpop = load("vpop_sisunatovir_VC.mat");

%%

all_params_ped = mitt3_all_params;
all_inoc_ped = mitt3_all_inoc;
vl_peak_pbos_ped = mitt3_all_peaks;

all_params_edp938 = edp938_vc_vpop.vpop_parameters;
vl_peak_pbos_edp938 = ones(size(all_params_edp938,1), 1)*edp938_vc_vpop.vpop_peak_pbo;
all_inoc_edp938 = edp938_vc_vpop.vpop_inoc;

all_params_edp938_2 = edp938_vc_vpop2.vpop_parameters;
vl_peak_pbos_edp938_2 = ones(size(all_params_edp938_2,1), 1)*edp938_vc_vpop2.vpop_peak_pbo;
all_inoc_edp938_2 = edp938_vc_vpop2.vpop_inoc;

all_params_sis = sis_vc_vpop.vpop_parameters;
all_inoc_sis = sis_vc_vpop.vpop_inoc;
vl_peak_pbos_sis = ones(length(all_inoc_sis), 1)*sis_vc_vpop.vpop_peak_pbo;

% PBO
%treatment times
length_treat_LD = [0];
length_treat_MD = [24 48 72 96]/24;

t_inoc_delay = 2.5;
edp938_fits = load('fit_EDP938_params.mat');
%EDP938 600 mg zelicapavir dose
drug_on = edp938_fits.drug1;
drug_on.Dose = 600*1e6; %ng
drug_on.LD = 600*1e6; %ng
drug_on.MD = 600*1e6; %ng

drug_on.Imax =  1;
drug_on.n = 3;
drug_on.EC50 =  (47.62/(9^(1/3))); %Median strain EC90 - used for patients

drug_on.dose_typ = 1;

drug = drug_on; %placebo with dose = 0
drug.Dose = 0;
drug.LD = 0; %ng
drug.MD = 0; %ng

xs_interp = 1:6:50*24;

t_treat_interv = [3,4,5,7]-1; % -1 to be from start of symptom onset, not peak

placebo_sims_peds = zeros(length(all_inoc_ped), length(xs_interp));
placebo_sims_edp938 = zeros(length(all_inoc_edp938), length(xs_interp));
placebo_sims_edp938_2 = zeros(length(all_inoc_edp938_2), length(xs_interp));
placebo_sims_sis = zeros(length(all_inoc_sis), length(xs_interp));
t_peaks_ped= zeros(length(all_inoc_ped), 1);
t_peaks_edp938= zeros(length(all_inoc_edp938), 1);
t_peaks_edp938_2= zeros(length(all_inoc_edp938_2), 1);
t_peaks_sis= zeros(length(all_inoc_sis), 1);

treat_sims3_peds = zeros(length(all_inoc_ped), length(xs_interp));
treat_sims3_edp938 = zeros(length(all_inoc_edp938), length(xs_interp));
treat_sims3_edp938_2 = zeros(length(all_inoc_edp938_2), length(xs_interp));
treat_sims3_sis = zeros(length(all_inoc_sis), length(xs_interp));

reduction_d3_ped = zeros(length(all_inoc_ped), length(t_treat_interv));
reduction_d3_edp938 = zeros(length(all_inoc_edp938), length(t_treat_interv));
reduction_d3_edp938_2 = zeros(length(all_inoc_edp938_2), length(t_treat_interv));
reduction_d3_sis = zeros(length(all_inoc_sis), length(t_treat_interv));

placebo_d3_ped = zeros(length(all_inoc_ped), length(t_treat_interv));
placebo_d3_edp938 = zeros(length(all_inoc_edp938), length(t_treat_interv));
placebo_d3_edp938_2 = zeros(length(all_inoc_edp938_2), length(t_treat_interv));
placebo_d3_sis = zeros(length(all_inoc_sis), length(t_treat_interv));

treatment_d3_ped = zeros(length(all_inoc_ped), length(t_treat_interv));
treatment_d3_edp938 = zeros(length(all_inoc_edp938), length(t_treat_interv));
treatment_d3_edp938_2 = zeros(length(all_inoc_edp938_2), length(t_treat_interv));
treatment_d3_sis = zeros(length(all_inoc_sis), length(t_treat_interv));

placebo_d5_ped = zeros(length(all_inoc_ped), length(t_treat_interv));
placebo_d5_edp938 = zeros(length(all_inoc_edp938), length(t_treat_interv));
placebo_d5_edp938_2 = zeros(length(all_inoc_edp938_2), length(t_treat_interv));
placebo_d5_sis = zeros(length(all_inoc_sis), length(t_treat_interv));

treatment_d5_ped = zeros(length(all_inoc_ped), length(t_treat_interv));
treatment_d5_edp938 = zeros(length(all_inoc_edp938), length(t_treat_interv));
treatment_d5_edp938_2 = zeros(length(all_inoc_edp938_2), length(t_treat_interv));
treatment_d5_sis = zeros(length(all_inoc_sis), length(t_treat_interv));

n_treat = length(t_treat_interv);

nominal_EC90 = 47.62; %nomimal, median strain EC90 - used for patient simulations

n_vpops = 5;
for j = 1:n_vpops
    
    %Load in respective VPops with appropriate scaling, & tolerances used for simulations
    if j == 1
        all_params = all_params_ped;
        all_inoc = all_inoc_ped;
        vl_peak_pbos = vl_peak_pbos_ped;
        scale_f = 1;
        tol = 1e-12;
    elseif j == 2
        all_params = all_params_edp938;
        all_inoc = all_inoc_edp938;
        vl_peak_pbos = vl_peak_pbos_edp938;
        scale_f = 1.4;
        tol = 1e-6;
    elseif j == 3
        all_params = all_params_sis;
        all_inoc =all_inoc_sis;
        vl_peak_pbos = vl_peak_pbos_sis;
        scale_f = 2.8;
        tol = 1e-6;
    elseif j == 4
        all_params = all_params_edp938_2;
        all_inoc = all_inoc_edp938_2;
        vl_peak_pbos = vl_peak_pbos_edp938_2;
        scale_f = 1.2;
        tol = 1e-6;
    end
    
    parfor i = 1:length(all_inoc)
        
        p_new = all_params(i,:);
        IC = zeros(8,1);
        IC(1) = 1;
        IC(4) = all_inoc(i);
        IC(6) = 1;
        
        VL_peak_PBO = vl_peak_pbos(i);
        
        if j == 1
            event_floor = 10;
            [curr_t_PBO,curr_y_PBO] = SolveBalances_rsv(0,50*24,0.001*24,p_new,...
                VL_peak_PBO, IC,drug, drug,tol,event_floor);
        else
            [curr_t_PBO,curr_y_PBO] = SolveBalances_rsv(0,50*24,0.001*24,p_new,...
                VL_peak_PBO, IC,drug, drug,tol);
        end
        
        model_VL_PBO = curr_y_PBO(1:end,4);
        
        model_VL_interp = interp1(curr_t_PBO,real(scale_f*log10(VL_peak_PBO*curr_y_PBO(:,4))), xs_interp);
        
        [curr_peak, curr_peak_ind] = max(real(scale_f*log10(VL_peak_PBO*curr_y_PBO(:,4))));
        
        t_peak =  curr_t_PBO(curr_peak_ind)/24;
        
        if j == 1
            placebo_sims_ped(i,:) = model_VL_interp;
            t_peaks_ped(i) = t_peak;
        elseif j == 2
            placebo_sims_edp938(i,:) = model_VL_interp;
            t_peaks_edp938(i) = t_peak;
        elseif j == 3
            placebo_sims_sis(i,:) = model_VL_interp;
            t_peaks_sis(i) = t_peak;
        elseif j == 4
            placebo_sims_edp938_2(i,:) = model_VL_interp;
            t_peaks_edp938_2(i) = t_peak;
        end
        
        for k=1:n_treat
            
            t_treat = t_treat_interv(k);
            loc_drug_on = drug_on;
            
            EC90 = nominal_EC90;
            loc_drug_on.EC50 =  (EC90/(9^(1/3))); %calcuate EC50, Hill n = 3
            
            loc_drug_on.dose_times_LD = (t_peak + t_treat + length_treat_LD)*24; %hrs
            loc_drug_on.dose_times_MD = (t_peak + t_treat + length_treat_MD)*24;
            loc_drug_on.dose_times = [loc_drug_on.dose_times_LD loc_drug_on.dose_times_MD];
            
            if j == 1
                event_floor = 10; %higher event handling floor cut-off for peds
                [curr_t_PBO2,curr_y_PBO2] = SolveBalances_rsv(0,50*24,0.001*24,p_new,...
                    VL_peak_PBO, IC,drug, loc_drug_on,tol, event_floor);
            else
                [curr_t_PBO2,curr_y_PBO2] = SolveBalances_rsv(0,50*24,0.001*24,p_new,...
                    VL_peak_PBO, IC,drug, loc_drug_on,tol);
            end

            model_VL_PBO = scale_f*curr_y_PBO2(1:end,4);
            
            model_VL_interp = scale_f*interp1(curr_t_PBO2,log10(VL_peak_PBO*curr_y_PBO2(:,4)), xs_interp);
            
            if k == 1
                if j == 1
                    treat_sims3_ped(i,:) = model_VL_interp;
                elseif j == 2
                    treat_sims3_edp938(i,:) = model_VL_interp;
                elseif j == 3
                    treat_sims3_sis(i,:) = model_VL_interp;
                elseif j == 4
                    treat_sims3_edp938_2(i,:) = model_VL_interp;
                end
            end
            
            [~, bl_ind] = min(abs(curr_t_PBO2/24 - t_peak - t_treat));
            [~, d3_post_dose_ind] = min(abs(curr_t_PBO2/24 - t_peak - t_treat - 2));
            [~, d5_post_dose_ind] = min(abs(curr_t_PBO2/24 - t_peak - t_treat - 4));
            
            
            day_x_treat_idxs = [bl_ind, d3_post_dose_ind,d5_post_dose_ind];
            vl_placebo = zeros(length(day_x_treat_idxs),1);
            vl_treat = zeros(length(day_x_treat_idxs),1);
            
            for kk = 1:length(vl_placebo)
                curr_idx = day_x_treat_idxs(kk);
                vl_placebo(kk) =  scale_f*log10(VL_peak_PBO*curr_y_PBO(curr_idx,4));
                vl_treat(kk) =  scale_f*log10(VL_peak_PBO*curr_y_PBO2(curr_idx,4));
            end
            
            
            blod = vl_treat < 1.;
            vl_treat(blod) = 1.;
            
            blod = vl_placebo < 1.;
            vl_placebo(blod) = 1.;

            if j == 1
                reduction_d3_ped(i,k) = vl_placebo(2)-vl_treat(2);
                reduction_d5_ped(i,k) = vl_placebo(3)-vl_treat(3);
                placebo_d3_ped(i,k) = vl_placebo(2);
                treatment_d3_ped(i,k) = vl_treat(2);
                placebo_d5_ped(i,k) = vl_placebo(3);
                treatment_d5_ped(i,k) = vl_treat(3);
            elseif j == 2
                reduction_d3_edp938(i,k) = vl_placebo(2)-vl_treat(2);
                reduction_d5_edp938(i,k) = vl_placebo(3)-vl_treat(3);
                placebo_d3_edp938(i,k) = vl_placebo(2);
                treatment_d3_edp938(i,k) = vl_treat(2);
                placebo_d5_edp938(i,k) = vl_placebo(3);
                treatment_d5_edp938(i,k) = vl_treat(3);
            elseif j == 3
                reduction_d3_sis(i,k) = vl_placebo(2)-vl_treat(2);
                reduction_d5_sis(i,k) = vl_placebo(3)-vl_treat(3);
                placebo_d3_sis(i,k) = vl_placebo(2);
                treatment_d3_sis(i,k) = vl_treat(2);
                placebo_d5_sis(i,k) = vl_placebo(3);
                treatment_d5_sis(i,k) = vl_treat(3);
            elseif j == 4
                reduction_d3_edp938_2(i,k) = vl_placebo(2)-vl_treat(2);
                reduction_d5_edp938_2(i,k) = vl_placebo(3)-vl_treat(3);
                placebo_d3_edp938_2(i,k) = vl_placebo(2);
                treatment_d3_edp938_2(i,k) = vl_treat(2);
                placebo_d5_edp938_2(i,k) = vl_placebo(3);
                treatment_d5_edp938_2(i,k) = vl_treat(3);
            end
            
            
        end
        
    end
end

%%
peak_shifted_xs_ped = xs_interp  - t_peaks_ped*24;
symptom_shifted_xs_ped = peak_shifted_xs_ped +24;

peak_shifted_xs_edp938 = xs_interp  - t_peaks_edp938*24;
symptom_shifted_xs_edp938 = peak_shifted_xs_edp938 +24;

peak_shifted_xs_edp938_2 = xs_interp  - t_peaks_edp938_2*24;
symptom_shifted_xs_edp938_2 = peak_shifted_xs_edp938_2 +24;

peak_shifted_xs_sis = xs_interp  - t_peaks_sis*24;
symptom_shifted_xs_sis = peak_shifted_xs_sis +24;

symptom_shifted_ts_interp = -10:0.1:21;
symptom_shifted_ys_pbo_ped = NaN*ones(length(t_peaks_ped), length(symptom_shifted_ts_interp));
symptom_shifted_ys_pbo_edp938 = NaN*ones(length(t_peaks_edp938), length(symptom_shifted_ts_interp));
symptom_shifted_ys_pbo_edp938_2 = NaN*ones(length(t_peaks_edp938_2), length(symptom_shifted_ts_interp));
symptom_shifted_ys_pbo_sis = NaN*ones(length(t_peaks_sis), length(symptom_shifted_ts_interp));

symptom_shifted_ys_treat3_ped = NaN*ones(length(t_peaks_ped), length(symptom_shifted_ts_interp));
symptom_shifted_ys_treat3_edp938 = NaN*ones(length(t_peaks_edp938), length(symptom_shifted_ts_interp));
symptom_shifted_ys_treat3_edp938_2 = NaN*ones(length(t_peaks_edp938_2), length(symptom_shifted_ts_interp));
symptom_shifted_ys_treat3_sis = NaN*ones(length(t_peaks_sis), length(symptom_shifted_ts_interp));

for i = 1:length(t_peaks_ped)
    if max(placebo_sims_ped(i,:)) < 1
        continue
    end
    
	symptom_shifted_ys_pbo_ped(i,:) = interp1(symptom_shifted_xs_ped(i,:)/24, placebo_sims_ped(i,:),symptom_shifted_ts_interp);
    
	symptom_shifted_ys_treat3_ped(i,:) = interp1(symptom_shifted_xs_ped(i,:)/24, treat_sims3_ped(i,:),symptom_shifted_ts_interp);
end

for i = 1:length(t_peaks_edp938)
    if max(placebo_sims_edp938(i,:)) < 1
        continue
    end
    
	symptom_shifted_ys_pbo_edp938(i,:) = interp1(symptom_shifted_xs_edp938(i,:)/24, placebo_sims_edp938(i,:),symptom_shifted_ts_interp);
    
 	symptom_shifted_ys_treat3_edp938(i,:) = interp1(symptom_shifted_xs_edp938(i,:)/24, treat_sims3_edp938(i,:),symptom_shifted_ts_interp);
end

for i = 1:length(t_peaks_sis)
    if max(placebo_sims_sis(i,:)) < 1
        continue
    end
    
	symptom_shifted_ys_pbo_sis(i,:) = interp1(symptom_shifted_xs_sis(i,:)/24, placebo_sims_sis(i,:),symptom_shifted_ts_interp);
	symptom_shifted_ys_treat3_sis(i,:) = interp1(symptom_shifted_xs_sis(i,:)/24, treat_sims3_sis(i,:),symptom_shifted_ts_interp);
end

for i = 1:length(t_peaks_edp938_2)
    if (max(placebo_sims_edp938_2(i,:)) < 1)
        continue
    end
    
 	symptom_shifted_ys_pbo_edp938_2(i,:) = interp1(symptom_shifted_xs_edp938_2(i,:)/24, placebo_sims_edp938_2(i,:),symptom_shifted_ts_interp);
	symptom_shifted_ys_treat3_edp938_2(i,:) = interp1(symptom_shifted_xs_edp938_2(i,:)/24, treat_sims3_edp938_2(i,:),symptom_shifted_ts_interp);
end


%% Ensure VL peaks increase across VPops
keep_ind_ped = ~(max(placebo_sims_ped') < 0.6);
keep_ind_edp938 = ~(max(real(placebo_sims_edp938)') < 0.6);
keep_ind_edp938_2 = ~(max(real(placebo_sims_edp938_2)') < 0.6);
keep_ind_sis = ~(max(real(placebo_sims_sis)') < 0.6);

red_d3_adult = [reduction_d3_edp938(keep_ind_edp938,:);reduction_d3_edp938_2(keep_ind_edp938_2,:); reduction_d3_sis(keep_ind_sis,:)];
red_d3_peds = reduction_d3_ped(keep_ind_ped,:);

red_d5_adult = [reduction_d5_edp938(keep_ind_edp938,:);reduction_d5_edp938_2(keep_ind_edp938_2,:); reduction_d5_sis(keep_ind_sis,:)];
red_d5_peds = reduction_d5_ped(keep_ind_ped,:);

plac_d3_adult = [placebo_d3_edp938(keep_ind_edp938,:);placebo_d3_edp938_2(keep_ind_edp938_2,:); placebo_d3_sis(keep_ind_sis,:)];
plac_d3_peds = placebo_d3_ped(keep_ind_ped,:);

treat_d3_adult = [treatment_d3_edp938(keep_ind_edp938,:);treatment_d3_edp938_2(keep_ind_edp938_2,:); treatment_d3_sis(keep_ind_sis,:)];
treat_d3_peds = treatment_d3_ped(keep_ind_ped,:);

plac_d5_adult = [placebo_d5_edp938(keep_ind_edp938,:);placebo_d5_edp938_2(keep_ind_edp938_2,:); placebo_d5_sis(keep_ind_sis,:)];
plac_d5_peds = placebo_d5_ped(keep_ind_ped,:);

treat_d5_adult = [treatment_d5_edp938(keep_ind_edp938,:);treatment_d5_edp938_2(keep_ind_edp938_2,:); treatment_d5_sis(keep_ind_sis,:)];
treat_d5_peds = treatment_d5_ped(keep_ind_ped,:);

%% Calculate peaks, slopes for final vpops
vpop_ped_slopes = ones(sum(keep_ind_ped),1);
vpop_ped_peaks = ones(sum(keep_ind_ped),1);
ped_idxs = find(keep_ind_ped);
for i = 1:sum(keep_ind_ped)
    try
        idx = ped_idxs(i);
        [slopels,b_ls,Vmid_index_ls,timevec_ls,Vmax_index]= terminal_slope_ls(symptom_shifted_ts_interp(91:10:end), symptom_shifted_ys_pbo_ped(idx,91:10:end));
        vpop_ped_slopes(i) = slopels;
        vpop_ped_peaks(i) = max(symptom_shifted_ys_pbo_ped(idx,:));
    catch
    end
    
end

edp938_slopes = ones(sum(keep_ind_edp938),1);
edp938_peaks = ones(sum(keep_ind_edp938),1);
edp938_idxs = find(keep_ind_edp938);
for i = 1:sum(keep_ind_edp938)
    try
        idx = edp938_idxs(i);
        [slopels,b_ls,Vmid_index_ls,timevec_ls,Vmax_index]= terminal_slope_ls(symptom_shifted_ts_interp(91:10:end), symptom_shifted_ys_pbo_edp938(idx,91:10:end));
        edp938_slopes(i) = slopels;
        edp938_peaks(i) = max(symptom_shifted_ys_pbo_edp938(idx,:));
    catch
    end
end

sis_slopes = ones(sum(keep_ind_sis),1);
sis_peaks = ones(sum(keep_ind_sis),1);
sis_idxs = find(keep_ind_sis);
for i = 1:sum(keep_ind_sis)
    try
        idx = sis_idxs(i);
        [slopels,b_ls,Vmid_index_ls,timevec_ls,Vmax_index]= terminal_slope_ls(symptom_shifted_ts_interp(91:10:end), symptom_shifted_ys_pbo_sis(idx,91:10:end));
        sis_slopes(i) = slopels;
        sis_peaks(i) = max(symptom_shifted_ys_pbo_sis(idx,:));
    catch
    end
    
end

edp938_2_slopes = ones(sum(keep_ind_edp938_2),1);
edp938_2_peaks = ones(sum(keep_ind_edp938_2),1);
edp938_2_idxs = find(keep_ind_edp938_2);
for i = 1:sum(keep_ind_edp938_2)
    try
        idx = edp938_2_idxs(i);
        [slopels,b_ls,Vmid_index_ls,timevec_ls,Vmax_index]= terminal_slope_ls(symptom_shifted_ts_interp(91:10:end), symptom_shifted_ys_pbo_edp938_2(idx,91:10:end));
        edp938_2_slopes(i) = slopels;
        edp938_2_peaks(i) = max(symptom_shifted_ys_pbo_edp938_2(idx,:));
    catch
    end
    
end

symptom_shifted_ys_pbo_healthy_adult = real([symptom_shifted_ys_pbo_edp938(keep_ind_edp938,:);symptom_shifted_ys_pbo_sis(keep_ind_sis,:) ]);
symptom_shifted_ys_treat3_healthy_adult = real([symptom_shifted_ys_treat3_edp938(keep_ind_edp938,:);symptom_shifted_ys_treat3_sis(keep_ind_sis,:) ]);

vpop_adult_slopes = [edp938_slopes;edp938_2_slopes; sis_slopes];
vpop_adult_peaks = [edp938_peaks;edp938_2_peaks; sis_peaks];

adult_vpop_slope_mean = -mean(vpop_adult_slopes);
adult_vpop_slope_std = std(vpop_adult_slopes);

ped_vpop_slope_mean = -mean(vpop_ped_slopes);
ped_vpop_slope_std =  std(vpop_ped_slopes);

adult_vpop_peak_mean = nanmean(max(symptom_shifted_ys_pbo_healthy_adult'));
adult_vpop_peak_std = nanstd(max(symptom_shifted_ys_pbo_healthy_adult'));

ped_vpop_peak_mean = nanmean(max(symptom_shifted_ys_pbo_ped(keep_ind_ped,:)'));
ped_vpop_peak_std = nanstd(max(symptom_shifted_ys_pbo_ped(keep_ind_ped,:)'));

%% Get dimensionalized/re-scaled parameters
%Constant parameters used for dimensionalizing
time = 24; %hours
T0 = 4e8; %cells
D00  = 1e6; %cells
 
n_vpops = 4;
scale_A = zeros(n_vpops,21);
for j = 1:n_vpops
    
    if j == 1
        all_params = all_params_ped;
        all_inoc = all_inoc_ped;
        vl_peak_pbos = vl_peak_pbos_ped;
        scale_f = 1;
        tol = 1e-12;
    elseif j == 2
        all_params = all_params_edp938;
        all_inoc = all_inoc_edp938;
        vl_peak_pbos = vl_peak_pbos_edp938;
        scale_f = 1.4;
        tol = 1e-6;
    elseif j == 3
        all_params = all_params_sis;
        all_inoc =all_inoc_sis;
        vl_peak_pbos = vl_peak_pbos_sis;
        scale_f = 2.8;
        tol = 1e-6;
    elseif j == 4
        all_params = all_params_edp938_2;
        all_inoc = all_inoc_edp938_2;
        vl_peak_pbos = vl_peak_pbos_edp938_2;
        scale_f = 1.2;
        tol = 1e-6;
    end
    
    A = [all_inoc all_params]';
 
    % %Mean VL - PBO
    V0 = (10^scale_f).*vl_peak_pbos(1);
    
    %Factors
    f1 = V0;
    f2 = time./V0;
    f3 = time;
    f4 = time;
    f5 = time;
    f6 = time*V0/T0;
    f7 = time;
    f8 = time/T0;
    f9 = time/D00;
    f10 = time;
    f11 = time*V0;
    f12 = time*D00;
    f13 = time;
    f14 = time;
    f15 = time/T0;
    f16 = time;
    f17 = 1/time;
    f18 = 1;
    f19 = time/T0;
    f20 = time./V0;
    f21 = time*T0;
 
    F = [f1 f2 f3 f4 f5 f6 f7 f8 f9 f10...
        f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21]';
 
    scale_A(j,:) = mean(F.*A,2);
%     scale_A_std(j,:) = std(F.*A,2);
 
    if j == 1
        all_params_scaled_ped = F.*A;
    elseif j == 2
        all_params_scaled_edp1 = F.*A;
    elseif j == 3
        all_params_scaled_sis = F.*A;
    elseif j == 4
        all_params_scaled_edp2 = F.*A;
    end
 
end

plac_d3_adult = [placebo_d3_edp938(keep_ind_edp938,:);placebo_d3_edp938_2(keep_ind_edp938_2,:); placebo_d3_sis(keep_ind_sis,:)];
plac_d3_peds = placebo_d3_ped(keep_ind_ped,:);

treat_d3_adult = [treatment_d3_edp938(keep_ind_edp938,:);treatment_d3_edp938_2(keep_ind_edp938_2,:); treatment_d3_sis(keep_ind_sis,:)];
treat_d3_peds = treatment_d3_ped(keep_ind_ped,:);

plac_d5_adult = [placebo_d5_edp938(keep_ind_edp938,:);placebo_d5_edp938_2(keep_ind_edp938_2,:); placebo_d5_sis(keep_ind_sis,:)];
plac_d5_peds = placebo_d5_ped(keep_ind_ped,:);

treat_d5_adult = [treatment_d5_edp938(keep_ind_edp938,:);treatment_d5_edp938_2(keep_ind_edp938_2,:); treatment_d5_sis(keep_ind_sis,:)];
treat_d5_peds = treatment_d5_ped(keep_ind_ped,:);


save("treatment_analysis", "adult_vpop_slope_mean","adult_vpop_slope_std",...
    "ped_vpop_slope_mean","ped_vpop_slope_std",...
    "adult_vpop_peak_mean","adult_vpop_peak_std",...
    "ped_vpop_peak_mean","ped_vpop_peak_std",...
    "symptom_shifted_ts_interp","symptom_shifted_ys_pbo_healthy_adult",...
    "symptom_shifted_ys_pbo_ped","red_d3_adult","red_d5_adult",...
    "red_d3_peds","red_d5_peds",...
    "plac_d3_adult","plac_d5_adult","plac_d3_peds","plac_d5_peds",...
    "treat_d3_adult","treat_d5_adult","treat_d3_peds","treat_d5_peds",...
     "all_params_scaled_ped", "all_params_scaled_edp1", ...
    "all_params_scaled_sis","all_params_scaled_edp2")

