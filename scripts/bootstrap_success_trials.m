load("treatment_analysis.mat")


Nvec      = [10 15 20 30 40 60 80 100 200];  % per-group sample sizes
nStudy    = 1000;    % number of simulated studies per N
alpha     = 0.05;
thresholds = [0.25 0.5 0.75 1];           % decision thresholds

rng(1);


figure('Color','w');
for i = 1:2
rng(1);

% Load in data for day 3 placebo and treated for adults and peds treated at 3D PSO
if i == 1
    subplot(2,1,1)
    x = plac_d3_peds(:,1) ; y = treat_d3_peds(:,1);
elseif i == 2
    subplot(2,1,2)
    x = plac_d3_adult(:,1) ; y = treat_d3_adult(:,1);
end


nN = numel(Nvec);
nT = numel(thresholds);

% Heatmap matrix: rows = thresholds, cols = N
Psuccess = zeros(nT, nN);

for iT = 1:nT
    thr = thresholds(iT);

    for iN = 1:nN
        N = Nvec(iN);
        success = false(nStudy,1);

        for s = 1:nStudy
            % --- Simulate one study ---
            xs = x(randi(numel(x), N, 1));
            ys = y(randi(numel(y), N, 1));

            % --- One-sided Welch CI vs threshold ---
            [~,~,ciObs] = ttest2(xs - thr, ys, ...
                'Vartype','unequal', ...
                'Tail','right', ...
                'Alpha',alpha);

            LCL = ciObs(1) + thr;
            success(s) = (LCL > thr);
        end

        % --- Operating characteristic ---
        Psuccess(iT,iN) = mean(success);
    end
end

%%
% Psuccess: nT x nN matrix
% thresholds: 1 x nT
% Nvec: 1 x nN



% Use evenly spaced x positions (1..nN) so columns are equal width
xPos = 1:numel(Nvec);

% Plot with explicit x/y coordinates
hImg = imagesc(xPos, thresholds, Psuccess);
ax = gca;
set(gca,'YDir','normal');  % so thresholds increase upward
ax.CLim = [0 1];
axis tight;

% --- Axes ticks/labels ---
% X: evenly spaced, but labeled with the actual N values
xticks(xPos);
xticklabels(string(Nvec));
xlabel('Sample Size (N)');

% Y: only show threshold values of interest as ticks
yticks(thresholds);
yticklabels(string(thresholds));
ylabel({'Success Threshold:'; 'log10 VL Placebo vs. Treatment'});

% --- Color/colormap ---
colormap(parula);
cb = colorbar;
cb.Label.String = 'P(Lower Confidence Limit > Threshold )';
% cb.Label.String = 'P( one-sided Welch CI LCL > threshold )';


if i == 1
    title({'Vitual Trial Success Heatmap:'; 'Day 3 Reduction in Pediatric Patients'});
elseif i == 2
    title({'Vitual Trial Success Heatmap:'; 'Day 3 Reduction in Healthy Adults'});
end

set(gca, ...
    'TickDir','out', ...
    'FontSize',12, ...
    'LineWidth',1, ...
    'Box','off');
    
% Annotate values inside cells (turn off if too busy)
doAnnotate = true;
if doAnnotate
    for iT = 1:numel(thresholds)
        for iN = 1:numel(Nvec)
            val = Psuccess(iT,iN);
            % choose contrasting text color based on cell value
            txtColor = 'w';
            if val > 0.75
                txtColor = 'k';
            end
            text(xPos(iN), thresholds(iT), sprintf('%.2f',val), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontSize',10, ...
                'Color',txtColor);
        end
    end
end

end
set(gcf,'position',[20, 20, 1000, 1000])
print(gcf,'figs/supp_bootstrap.png','-dpng','-r300');

