clear; 
clc;

SUBJ = "Rawan";
YYYY = 2022;
MM = 9;
DD = 19;
BLOCK = 2;
ELIM = [20 80];

x = io.load_tmsi(SUBJ, YYYY, MM, DD, ["A", "B"], BLOCK, '.poly5');

% Get the onset difference (seconds)
dt = seconds(datetime(x(2).date) - datetime(x(1).date));
fs = x(1).sample_rate;

% How many samples ahead of "A" (x(1)) is "B" (x(2))?
n_sync_samples = round(fs*dt);

% Remove the first `n_sync_samples` from "B".
iStart = n_sync_samples+1;
iEnd = min(size(x(2).samples,2)-n_sync_samples, size(x(1).samples, 2));
xA = x(1).samples(1:64,1:iEnd);
xB = x(2).samples(1:64, (n_sync_samples+1):(n_sync_samples+iEnd));


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
default.savefig(fig, "R:/NMLShare/generated_data/human/TC/MUAPs/RMS");

%%
xA((rA < ELIM(1)) | (rA > ELIM(2)), :) = [];
xB((rB < ELIM(1)) | (rB > ELIM(2)), :) = [];

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
default.savefig(fig, "R:/NMLShare/generated_data/human/TC/MUAPs/Pareto");

%% The first PC is just noise.
fig = figure('Name', 'PC-1 Noise', 'Color', 'w');
spectrogram(score(:,1), 512, 256, 0:500, 4000, 'yaxis');
title("PC-1 is mains noise");
xlabel("Time");
ylabel("Amplitude (\muV)");
default.savefig(fig, "R:/NMLShare/generated_data/human/TC/MUAPs/PC1-Noise");

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
default.savefig(fig, "R:/NMLShare/generated_data/human/TC/MUAPs/PC2-5");

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
default.savefig(fig, "R:/NMLShare/generated_data/human/TC/MUAPs/PC6-11");


