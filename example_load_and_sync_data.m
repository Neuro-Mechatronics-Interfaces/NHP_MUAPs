close all force;
clear; 
clc;

AUTO_SAVE_FIGURES = true; 
SUBJ = "Rawan";
YYYY = 2022;
MM = 10;
DD = 7;
BLOCK = 0;
ELIM = [20 80];
dt = datetime('now', 'TimeZone','local', 'Format', 'yyyy-MM-dd_hh-mm-ss');

x = io.load_tmsi(SUBJ, YYYY, MM, DD, ["A", "B"], BLOCK, '.poly5');

% Get the onset difference (seconds)
tdiff = seconds(datetime(x(2).date) - datetime(x(1).date));
fs = x(1).sample_rate;

% How many samples ahead of "A" (x(1)) is "B" (x(2))?
n_sync_samples = round(fs*tdiff);

% Remove the first `n_sync_samples` from "B".
iStart = n_sync_samples+1;
iEnd = min(size(x(2).samples,2)-n_sync_samples, size(x(1).samples, 2));
a_samples = 1:iEnd;
xA = x(1).samples(1:64, a_samples);
b_samples = (n_sync_samples+1):(n_sync_samples+iEnd);
xB = x(2).samples(1:64, b_samples);


fig = figure('Name', 'Sample RMS', 'Color', 'w', 'Position', [4014 -715  1317  420]);
L = tiledlayout(fig, 1, 2);
title(L, sprintf('%s\\_%04d\\_%02d\\_%02d\\_%d',SUBJ,YYYY,MM,DD,BLOCK), ...
    'FontName','Tahoma','Color','k'); 
ax = nexttile(L);
set(ax,'NextPlot','add','XColor','k','YColor','k','FontName','Tahoma');
title(ax,'A', 'FontName','Tahoma','Color',[0.3 0.3 0.3]);
rA = rms(xA,2)./1e3;
histogram(ax, rA);
xlabel(ax, 'RMS (mV)', 'FontName','Tahoma');
xline(ax, ELIM, 'r:', 'Label', 'Exclusions');
xlim(ax,[0 100]);
ax = nexttile(L);
set(ax,'NextPlot','add','XColor','k','YColor','k','FontName','Tahoma');
title(ax,'B', 'FontName','Tahoma','Color',[0.3 0.3 0.3]);
rB = rms(xB,2)./1e3;
histogram(ax, rB);
xlabel(ax, 'RMS (mV)', 'FontName','Tahoma');
xlim(ax,[0 100]);
xline(ax, ELIM, 'r:', 'Label', 'Exclusions');

if AUTO_SAVE_FIGURES
    default.savefig(fig, fullfile("R:/NMLShare/generated_data/human/TC/MUAPs", string(dt) ,"RMS"));
end

%%
all_ch = 1:64;
not_in_A = all_ch((rA < ELIM(1)) | (rA > ELIM(2)));
not_in_B = all_ch((rB < ELIM(1)) | (rB > ELIM(2)));
xA(not_in_A, :) = [];
xB(not_in_B, :) = [];
in_A = setdiff(all_ch, not_in_A);
in_B = setdiff(all_ch, not_in_B) + 64;

%%
yA = xA - mean(xA,2);
% yA = yA - mean(yA,1);

yB = xB - mean(xB,2);
% yB = yB - mean(yB,1);

[b,a] = butter(4, ([25 400])./4000, 'bandpass');
zB = filtfilt(b,a,yB')';
zA = filtfilt(b,a,yA')';

zB = zB - mean(zB,1);
zA = zA - mean(zA,1);

t = 0:(1/fs):((iEnd-1)/fs);
Y = [zA; zB];

%%
[coeff,score,latent,tsquared,explained,mu] = pca(Y');
fig = figure('Name', 'Pareto Plot', 'Color', 'w');
ax = axes(fig, 'NextPlot', 'add');
pareto(ax, explained);
if AUTO_SAVE_FIGURES
    default.savefig(fig, fullfile("R:/NMLShare/generated_data/human/TC/MUAPs", string(dt), "Pareto"));
end

%% The first PC is just noise.
fig = figure('Name', 'PC-1 Noise', 'Color', 'w');
spectrogram(score(:,1), 512, 256, 0:500, 4000, 'yaxis');
title("PC-1 is mains noise");
xlabel("Time");
ylabel("Amplitude (\muV)");
default.savefig(fig, fullfile("R:/NMLShare/generated_data/human/TC/MUAPs", string(dt), "PC1-Noise"));

%% Look next at PC-2 thru PC-5
fig = figure('Name', 'PC-2 thru PC-5', 'Color', 'w', 'Position', [4008        -996        1067         892]);
L = tiledlayout(fig, 'flow');
for ii = 2:5
    ax = nexttile(L);
    set(ax, 'NextPlot', 'add', 'FontName', 'Tahoma');
    axes(ax);
    spectrogram(score(:,ii), 512, 256, 0:500, 4000, 'yaxis');
    
    xlabel(ax, 'Time (mins)', 'FontName', 'Tahoma');
    ylabel(ax, 'Frequency (Hz)', 'FontName', 'Tahoma');
    ylim(ax, [0 300]);
    if ii == 5
        title(ax, sprintf('PC-%d (noise)', ii), 'FontName', 'Tahoma', 'Color', 'r');
    else
        title(ax, sprintf('PC-%d', ii), 'FontName', 'Tahoma', 'Color', 'k');
    end
end
if AUTO_SAVE_FIGURES
    default.savefig(fig, fullfile("R:/NMLShare/generated_data/human/TC/MUAPs", string(dt), "PC2-5"));
end

%% Look next at PC-6 thru PC-11
fig = figure('Name', 'PC-6 thru PC-11', 'Color', 'w', 'Position', [4008        -996        1067         892]);
L = tiledlayout(fig, 'flow');
for ii = 6:11
    ax = nexttile(L);
    set(ax, 'NextPlot', 'add', 'FontName', 'Tahoma');
    axes(ax);
    spectrogram(score(:,ii), 512, 256, 0:500, 4000, 'yaxis');
    
    xlabel(ax, 'Time (mins)', 'FontName', 'Tahoma');
    ylabel(ax, 'Frequency (Hz)', 'FontName', 'Tahoma');
    ylim(ax, [0 300]);
    if ii == 5
        title(ax, sprintf('PC-%d (noise)', ii), 'FontName', 'Tahoma', 'Color', 'r');
    else
        title(ax, sprintf('PC-%d', ii), 'FontName', 'Tahoma', 'Color', 'k');
    end
end
if AUTO_SAVE_FIGURES
    default.savefig(fig, fullfile("R:/NMLShare/generated_data/human/TC/MUAPs", string(dt), "PC6-11"));
end
%% Use thresholding on Y
THRESH_RMS = 1.5;
N_SAMPLES_PP = 1000;
SMOOTHING = 0.9;

fs = x(1).sample_rate;
pp = cell(size(Y,1),1);
tp = linspace(0, max(t), N_SAMPLES_PP);
Yp = nan(N_SAMPLES_PP, size(Y, 1));
% All I am doing is fitting a piecewise cubic polynomial to the
% instantaneous "lambda" (intensity) values (which are equivalent to
% 1/samples between each "spike"):
for ii = 1:size(Y, 1)
    s = abs(Y(ii,:)) > (rms(abs(Y(ii,:)), 2)*THRESH_RMS);
    ts = [1, find(s)./fs, size(S,2)]; % MUAP times in seconds
    dts = [0, 1/ts(2), 1./diff(ts(2:(end-1))), 0];
    i_remove = dts >= (0.35*fs); % This is to make sure that I didn't get spikes that were too close together
    ts(i_remove) = [];
    dts(i_remove) = [];
    pp{ii} = csaps(ts, dts, SMOOTHING); % Setting SMOOTHING lower makes the rate traces smoother.
    Yp(:, ii) = max(ppval(pp{ii}, tp), 0);
end
fig = figure('Name', 'MUAPs Rate', 'Color', 'w', ...
    'Position', [4014        -842         868         420]);
ax = axes(fig,'NextPlot', 'add', 'XColor', 'k', 'YColor', 'k');
plot(ax, tp, Yp);
title(ax, 'MUAPs Rates', 'FontName', 'Tahoma', 'Color', 'k');
subtitle(ax, sprintf('Piecewise Cubic Spline (Smoothing p = %3.2f)', SMOOTHING), 'FontName', 'Tahoma', 'Color', [0.5 0.5 0.5]); 
ylabel(ax, 'MUAPs/sec', 'FontName','Tahoma','Color','k');
xlabel(ax, 'Time (s)', 'FontName','Tahoma','Color','k');
if AUTO_SAVE_FIGURES
    default.savefig(fig, fullfile("R:/NMLShare/generated_data/human/TC/MUAPs", string(dt), "MUAP-Rates"));
end
%% Now, decompose those rates
[coeff,score,latent,tsquared,explained,mu] = pca(Yp);
fig = figure('Name', 'Rate PCs Pareto', 'Color', 'w', ...
    'Position', [4139        -838        1405         420]);

L = tiledlayout(fig, 2, 5);
ax = nexttile(L, 1, [1 1]);
pareto(ax, explained);
title("MUAPs Pareto", 'FontName', 'Tahoma', 'Color', 'k');
xlabel("Rate PC", 'FontName', 'Tahoma', 'Color', 'k');
h = findobj(ax.Children, 'Type', 'bar');
h.CData(1:6, :) = ax.ColorOrder(1:6,:);
h.FaceColor = 'flat';
h = findobj(ax.Children', 'Type', 'line');
h.LineWidth = 1.5;
h.Color = 'k';
h.MarkerIndices = 6;
h.Marker = 'o';
h.MarkerFaceColor = 'm';
cs = cumsum(explained);
text(ax, 6, cs(6)-10, sprintf('%5.1f%%', cs(6)), 'FontName', 'Tahoma', 'Color', 'm', 'FontWeight', 'bold');

ax = nexttile(L, 6, [1 1]);
plot(ax, tp, score(:, 1:6), 'LineWidth', 1.5);
xlabel(ax, "Time (s)", 'FontName', 'Tahoma', 'Color', 'k');
ylabel(ax, "Rate (MUAP/s)", 'FontName', 'Tahoma', 'Color', 'k');

ax = nexttile(L, 2, [2 4]);
set(ax, 'NextPlot', 'add', 'XColor', 'k', 'YColor', 'k', 'FontName', 'Tahoma', 'FontSize', 14);
S = score(:, 1:6) - min(score(:, 1:6), [], 1); 
S = S ./ max(S, [], 1);
for ii = 1:6
    scatter(ax, tp(iMax == ii), iMax(iMax == ii), 'Marker', 's', 'MarkerFaceColor', 'flat'); 
end
title(ax, 'Motor Decomp "State"', 'FontName', 'Tahoma', 'Color', 'k');
ylim(ax, [0 7]);
xlabel(ax, 'Time (s)', 'FontName', 'Tahoma', 'Color', 'k');
if AUTO_SAVE_FIGURES
    default.savefig(fig, fullfile("R:/NMLShare/generated_data/human/TC/MUAPs", string(dt), "MUAP-Rates-Pareto-State"));
end

%% Spatial distribution of PCs on electrodes
coeff_big = nan(128,6);
for ii = 1:6
    coeff_big(in_A, ii) = coeff(1:(numel(in_A)),ii);
    coeff_big(in_B, ii) = coeff((numel(in_A)+1):end, ii);
end
fig = figure('Name', 'Top-6 Rate PC Coefficients', 'Color', 'w', ...
    'Position', [4091        -821         920         420]);
L = tiledlayout(fig, 'flow');

for ii = 1:6
    ax = nexttile(L);
    set(ax, 'NextPlot', 'add', 'YDir', 'normal', ...
        'XTickLabel', 1:16:128, ...
        'XTick', 1:2:16, ...
        'XLim', [0.5 16.5], 'YLim', [0.5 8.5]);
    imagesc(ax, [1 16], [1 8], reshape(coeff_big(:, ii), 8, 16));
    colorbar(ax);
    xline(ax, 8.5, 'LineWidth', 2);
    title(ax, sprintf('PC-%d', ii), 'FontName', 'Tahoma', 'Color', 'k');
    ylabel(ax, 'Channel', 'FontName', 'Tahoma', 'Color', 'k');
end
title(L, 'Rawan MUAPs Rate Coefficients', 'FontName', 'Tahoma', 'Color', 'k');
subtitle(L, sprintf('(%5.1f%% of rate estimates explained)', cs(6)), 'FontName', 'Tahoma', 'Color', [0.5 0.5 0.5]);
if AUTO_SAVE_FIGURES
    default.savefig(fig, fullfile("R:/NMLShare/generated_data/human/TC/MUAPs", string(dt), "MUAP-Rates-Components"));
end
%%
