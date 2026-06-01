function drug_conc =PK_check(drug,t)
sd_LD=0;
sd_MD = 0;


for i=1:numel(drug.dose_times_LD)

drug.Dose = drug.LD/drug.mwa;
temp_LD = single_dose_contribution(drug,drug.dose_times_LD(i), t);
sd_LD = sd_LD+temp_LD;
end

for j = 1:numel(drug.dose_times_MD)
drug.Dose = drug.MD/drug.mwa;
temp_MD = single_dose_contribution(drug,drug.dose_times_MD(j), t);
sd_MD = sd_MD + temp_MD;
end

drug_conc=sd_LD + sd_MD;

end


function sd = single_dose_contribution(drug,dose_time, t)
    %absoprtion rate constant
    k_a    = drug.k_a;%1/h
    %elimination rate
    k = drug.k_e; %hours
    % Volume of distribution
    Vd     = drug.Vd; %L
    %dose

    %duration of infusion
    T = t-dose_time;

    T(T<0) = 1e10;

    %Single dose formula - plasma concentration for oral administration
    %Initial concentration
    C0 = (drug.Dose/Vd); %nmol/L
    %Advantage for full PK model: simulate any dose
    sd = C0.*(k_a/(k_a-k)).*(exp(-k.*(T))-exp(-k_a.*(T)));

end