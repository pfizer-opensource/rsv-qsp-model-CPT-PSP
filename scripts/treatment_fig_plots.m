load("study_analysis.mat")
load("treatment_analysis.mat")

color_peds = [32,178,170]/256; %[4, 95, 95]/256;
color_adults = [0, 0, 0]/256; %[128, 0, 128]/256;

figure
subplot(1,3,2)

bar(18, mean(adult_slopes),"FaceColor",color_adults);
hold on 
bar(19, mean(ped_slopes), "FaceColor",color_peds);
hold on 

bar(20.5, adult_vpop_slope_mean,"FaceColor",color_adults,"FaceAlpha",0.3);
hold on 
bar(21.5, ped_vpop_slope_mean, "FaceColor",color_peds,"FaceAlpha",0.3);
hold on 


errorbar(18., mean(adult_slopes), std(adult_slopes), 'w','LineWidth',3);
hold on
errorbar(18, mean(adult_slopes), std(adult_slopes), 'k','LineWidth',2);
hold on
errorbar(19, mean(ped_slopes), std(ped_slopes), 'k','LineWidth',2);
hold on

errorbar(18.5,1.35,0, 0,0.5,0.5,'k')
hold on
text(18.45,1.4, '*','FontSize',20);
hold on

errorbar(20.5, adult_vpop_slope_mean, adult_vpop_slope_std, 'k','LineWidth',2);
hold on
errorbar(21.5, ped_vpop_slope_mean, ped_vpop_slope_std, 'k','LineWidth',2);
hold on



xticks([18, 19, 20.5, 21.5])
xticklabels({"Adult Challenge Studies (N = 6)", "Pediatric Studies (N = 9)",...
    "Simulated Healthy Adult VPop", "Simulated Pediatric Patient VPop"})

xtickangle(30)
ylim([0 1.6])

% title("VL Slopes Max to Half-Max (Manually Calculated)")
% title('Post-Peak Viral Load Slopes')
%ylabel("Negative terminal slope of viral shedding [log10 VL units/day]")
ylabel("Neg. Post-Peak Slope [log10 copies/(mL day)]")
grid on 
set(gca,'Fontsize',14)

subplot(1,3,1)
bar(18, mean(adult_peaks),"FaceColor",color_adults);
hold on 
bar(19, mean(ped_peaks), "FaceColor",color_peds);
hold on 

bar(20.5, adult_vpop_peak_mean,"FaceColor",color_adults,"FaceAlpha",0.3);
hold on 
bar(21.5, ped_vpop_peak_mean, "FaceColor",color_peds,"FaceAlpha",0.3);
hold on 

errorbar(18., mean(adult_peaks), std(adult_peaks), 'w','LineWidth',3);
hold on
errorbar(18, mean(adult_peaks), std(adult_peaks), 'k','LineWidth',2);
hold on
errorbar(19, mean(ped_peaks), std(ped_peaks), 'k','LineWidth',2);
hold on

errorbar(18.5,8,0, 0,0.5,0.5,'k')
hold on
text(18.45,8.2, '*','FontSize',20);
hold on

errorbar(20.5, adult_vpop_peak_mean, adult_vpop_peak_std, 'k','LineWidth',2);
hold on
errorbar(21.5, ped_vpop_peak_mean, ped_vpop_peak_std, 'k','LineWidth',2);
hold on


xticks([18, 19, 20.5, 21.5])
xticklabels({"Adult Challenge Studies (N = 6)", "Pediatric Studies (N = 9)",...
    "Simulated Healthy Adult VPop", "Simulated Pediatric Patient VPop"})

xtickangle(30)

% title("Peak or Baseline Viral Load")
ylabel("Peak Viral Load [log10 copies/mL]")
grid on 
set(gca,'Fontsize',14)

subplot(1,3,3)
p1 = plot(symptom_shifted_ts_interp, nanmean(symptom_shifted_ys_pbo_ped,1),'color', color_peds,'LineWidth',3); hold on

x_temp = [symptom_shifted_ts_interp(10:end)'; fliplr(symptom_shifted_ts_interp(10:end))'];
inBetween = [ prctile(symptom_shifted_ys_pbo_ped(:,10:end),2.5,1)'; ...
    fliplr(prctile(symptom_shifted_ys_pbo_ped(:,10:end),97.5,1))' ];
p6 = fill(x_temp,inBetween,'cyan');
p6.FaceColor = color_peds;
p6.FaceAlpha = 0.2;
p6.EdgeColor = 'none';
hold on

p2 = plot(symptom_shifted_ts_interp, nanmean(symptom_shifted_ys_pbo_healthy_adult,1),'color', color_adults,'LineWidth',3); hold on

x_temp = [symptom_shifted_ts_interp(50:end-50)'; fliplr(symptom_shifted_ts_interp(50:end-50))'];
inBetween = [ prctile(symptom_shifted_ys_pbo_healthy_adult(:,50:end-50),2.5,1)'; ...
    fliplr(prctile(symptom_shifted_ys_pbo_healthy_adult(:,50:end-50),97.5,1))' ];
p7 = fill(x_temp,inBetween,'cyan');
p7.FaceColor = color_adults;
p7.FaceAlpha = 0.2;
p7.EdgeColor = 'none';
hold on

xlabel("Time Since Symptom Onset (Days)")
ylabel("Viral Load (log10 copies/mL)")
ylim([1 10])
xlim([0 21])
legend([p1(1), p6(1), p2(1), p7(1)],"Pediatric Patient VPop\newlineSimulated Mean PBO","Pediatric 95% PI",...
    "Healthy Adult VPop\newlineSimulated Mean PBO","Healthy Adult 95% PI")

set(gca,'Fontsize',14)

set(gcf,'position',[20, 20, 1800, 500])
print(gcf,'figs/fig3.png','-dpng','-r300');


%%
figure
for k = 1:2
    subplot(1,2,k)
    if k == 1
        data = [mean(red_d3_adult,1); mean(red_d3_peds,1)]';
        errors = [std(red_d3_adult,1)/sqrt(length(red_d3_adult)); std(red_d3_peds,1)/sqrt(length(red_d3_peds))]';
    elseif k == 2
        data = [mean(red_d5_adult,1); mean(red_d5_peds,1)]';
        errors = [std(red_d5_adult,1)/sqrt(length(red_d5_adult)); std(red_d5_peds,1)/sqrt(length(red_d5_peds))]';
    end
    % Number of groups and number of bars in each group
    [numGroups, numBars] = size(data);
    % Plot the grouped bar chart
    hold on;
    hBar = bar(data);
    % Get the x positions of the bars
    xBar = nan(numGroups, numBars);
    for i = 1:numBars
        xBar(:, i) = hBar(i).XEndPoints;
    end
    % Plot the error bars
    for i = 1:numGroups
        for j = 1:numBars
            errorbar(xBar(i, j), data(i, j), errors(i, j), 'k', 'linestyle', 'none', 'CapSize', 10);
        end
    end
    % Labels and title
    xlabel('Zelicapavir Treatment Intervention Timing');
    legend('Healthy Adult', 'Pediatric Patients');
    hold off;
    
    hBar(1).FaceColor = color_adults;
    hBar(2).FaceColor = color_peds;
    %hBar(2).FaceAlpha = 0.3;
    
    grid on
    ylim([0 2])
    xticks([1 2 3 4])
    xticklabels({'3D PSO','4D PSO','5D PSO','7D PSO'})
    set(gca,'FontSize',14)
    
end
subplot(1,2,1)
ylabel('\Delta Viral Load Reduction Day 3 (log10 copies/mL)');
% title('Viral Load Reduction Day 3');
subplot(1,2,2)
ylabel('\Delta Viral Load Reduction Day 5 (log10 copies/mL)');
% title('Viral Load Reduction Day 5');

set(gcf,'position',[20, 20, 1200, 500])
print(gcf,'figs/fig6.png','-dpng','-r300');


%%
paramlistplain = ["Inoc","beta","k","deltaI","c","p","deltaR","kappa0","kappa1","kappa2","nu","K","deltaJ","lambda",...
     "sigmaI","deltaD","tau","m","sigmaJ","sigmaV","kalpha"];
 
figure
for i = 1:21
    subplot(7,3,i)
    
   if  any(i == [10, 13, 17, 18])
    mean_vals = [mean([all_params_scaled_sis(i,:), all_params_scaled_edp1(i,:), all_params_scaled_edp2(i,:)]), ...
        mean(all_params_scaled_ped(i,:))];
    std_vals = [std([all_params_scaled_sis(i,:), all_params_scaled_edp1(i,:), all_params_scaled_edp2(i,:)]), ...
        std(all_params_scaled_edp2(i,:))];
   else
    
    mean_vals = [mean(log10([all_params_scaled_sis(i,:), all_params_scaled_edp1(i,:), all_params_scaled_edp2(i,:)])), ...
        mean(log10(all_params_scaled_ped(i,:)))];
    std_vals = [std(log10([all_params_scaled_sis(i,:), all_params_scaled_edp1(i,:), all_params_scaled_edp2(i,:)])), ...
        std(log10(all_params_scaled_ped(i,:)))];
   end
    
    bh = bar(1:2, mean_vals,'facecolor', 'flat'); hold on
    
    clr = [0 0 0;
        3/255 177/255 169/255];
 
    bh.CData = clr;
    
    errorbar(1:2, mean_vals,std_vals, 'k', 'LineStyle','none')
    ylabel("log10 " + paramlistplain(i))
    if any(i == [10, 13, 17, 18])
        ylabel( paramlistplain(i))
    end
    
 
    xticks(1:2)
    xticklabels(["Healthy Adult","Pediatric Patients"])
    xtickangle(20)
    
    set(gca, 'FontSize', 12)
    grid on
    
end
 
set(gcf,'position',[20, 20, 1200, 900])
print(gcf,'figs/supp_params.png','-dpng','-r300');