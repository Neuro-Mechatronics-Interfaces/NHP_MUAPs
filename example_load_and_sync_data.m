%EXAMPLE_LOAD_AND_SYNC_DATA  Load pilot Rawan recordings for MUAP initial tests.
close all force;
clear; 
clc;

AUTO_SAVE_FIGURES = true; 
SUBJ = "Rawan";
YYYY = 2022;
MM = 10;
DD = 10;
BLOCK = 0:4;
ELIM = [20 80];
dt = datetime('now', 'TimeZone','local', 'Format', 'yyyy-MM-dd_hh-mm-ss');

for iB = 1:numel(BLOCK)
    x = io.load_tmsi(SUBJ, YYYY, MM, DD, ["A", "B"], BLOCK(iB), '.poly5');
    block = sprintf("%s_%04d_%02d_%02d_%d", SUBJ, YYYY, MM, DD, BLOCK(iB));
    
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
    title(L, sprintf('%s\\_%04d\\_%02d\\_%02d\\_%d',SUBJ,YYYY,MM,DD,BLOCK(iB)), ...
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
        default.savefig(fig, fullfile("R:/NMLShare/generated_data/human/TC/MUAPs", block, string(dt), "RMS"));
    end
    
    %%
    all_ch = 1:64;
%     not_in_A = all_ch((rA < ELIM(1)) | (rA > ELIM(2)));
%     not_in_B = all_ch((rB < ELIM(1)) | (rB > ELIM(2)));
    not_in_A = [];
    not_in_B = [];
    xA(not_in_A, :) = [];
    xB(not_in_B, :) = [];
    in_A = setdiff(all_ch, not_in_A);
    in_B = setdiff(all_ch, not_in_B) + 64;
    n_samples = size(xA, 2);
    
    %%
    xA = xA - mean(xA,2);
    xB = xB - mean(xB,2);
    
    [b,a] = butter(4, ([25 400])./4000, 'bandpass');
    xB = filtfilt(b,a,xB')';
    xA = filtfilt(b,a,xA')';
    
    xB = xB - mean(xB,1);
    xB = reshape(xB, 8, 8, n_samples);
    xB = reshape([nan(1,8,n_samples); diff(xB, 2, 1); nan(1,8,n_samples)], 64, n_samples);
    xA = xA - mean(xA,1);
    xA = reshape(xA, 8, 8, n_samples);
    xA = reshape([nan(1,8,n_samples); diff(xA, 2, 1); nan(1,8,n_samples)], 64, n_samples);
    excvec = union(1:8:57, 8:8:64);
    in_A(excvec) = [];
    in_B(excvec) = [];

    t = 0:(1/fs):((iEnd-1)/fs);
    clear x;
    xA = xA(in_A,:);
    xB = xB(in_B-64,:);
    Y = [xA; xB];
    
    %% Use thresholding on Y
    THRESH_RMS = 2.0;
    N_SAMPLES_PP = 5000;
    SMOOTHING = 0.98;
    MAX_MUAPS_PER_SEC = 250; 
    
    pp = cell(size(Y,1),1);
    tp = linspace(0, max(t), N_SAMPLES_PP);
    Yp = nan(N_SAMPLES_PP, size(Y, 1));
    % All I am doing is fitting a piecewise cubic polynomial to the
    % instantaneous "lambda" (intensity) values (which are equivalent to
    % 1/samples between each "spike"):
    ts = cell(size(Y,1),1);
    for ii = 1:size(Y, 1)
        s = abs(Y(ii,:)) > (rms(abs(Y(ii,:)), 2)*THRESH_RMS);
        ts{ii} = [1./fs, find(s)./fs, size(Y,2)./fs]; % MUAP times in seconds
        dts = [0, 1/ts{ii}(2), 1./diff(ts{ii}(2:(end-1))), 0];
        i_remove = dts >= MAX_MUAPS_PER_SEC; % This is to make sure that I didn't get spikes that were too close together
        ts{ii}(i_remove) = [];
        dts(i_remove) = [];
        pp{ii} = csaps(ts{ii}, dts, SMOOTHING); % Setting SMOOTHING lower makes the rate traces smoother.
        Yp(:, ii) = max(ppval(pp{ii}, tp), 0);
    end
% 
%     nbuf = round(0.05 * N_SAMPLES_PP);
%     vec_keep = nbuf:(N_SAMPLES_PP - nbuf + 1);
%     
%     tp = tp(vec_keep);
%     Yp = Yp(vec_keep,:);

    iExc = isoutlier(rms(Yp,1));
    in_both = [in_A, in_B];
    in_both(iExc) = [];

    Yp = Yp(:, ~iExc);
    ts = ts(~iExc);  

    fig = figure('Name', 'MUAPs Rate', 'Color', 'w', ...
        'Position', [4014        -842         868         720]);
    L = tiledlayout(fig, 3, 1);
    ax = nexttile(L, 2, [1 1]);
    set(ax, 'NextPlot', 'add', 'XColor', 'k', 'YColor', 'k', 'XLim', [t(1), t(end)]);
    for ii = 1:numel(ts)
        tmpx = [ts{ii}; ts{ii}; nan(1, numel(ts{ii}))];
        tmpy = ones(3, numel(ts{ii})) .* [in_both(ii)-0.33; in_both(ii)+0.33; nan];
        line(ax, tmpx(:), tmpy(:), 'LineWidth', 1.25, 'Color', ax.ColorOrder(rem(ii-1, size(ax.ColorOrder,1))+1,:));
    end
    title(ax, 'MUAPs Raster', 'FontName', 'Tahoma', 'Color', 'k');
    ylabel(ax, 'Channel', 'FontName','Tahoma','Color','k');
    xlabel(ax, 'Time (s)', 'FontName','Tahoma','Color','k');

    ax = nexttile(L, 3, [1 1]);
    set(ax, 'NextPlot', 'add', 'XColor', 'k', 'YColor', 'k', 'XLim', [t(1), t(end)]);
    plot(ax, tp, Yp);
    title(ax, 'Cubic Spline Rates (MUAPs)', 'FontName', 'Tahoma', 'Color', 'k');
    subtitle(ax, sprintf('Piecewise Cubic Spline (Smoothing p = %3.2f)', SMOOTHING), 'FontName', 'Tahoma', 'Color', [0.5 0.5 0.5]); 
    ylabel(ax, 'Rate (MUAPs/sec)', 'FontName','Tahoma','Color','k');
    xlabel(ax, 'Time (s)', 'FontName','Tahoma','Color','k');

    ax = nexttile(L, 1, [1 1]);
    set(ax, 'NextPlot', 'add', 'XColor', 'k', 'YColor', 'k', 'XLim', [t(1), t(end)]);
    E = Y(~iExc, vec_keep);
    E = E./max(abs(E), [], 2) + in_both';
    plot(ax, tp, E);
    title(ax, strrep(block, '_', '\_'), 'FontName', 'Tahoma', 'Color', 'k');
    subtitle(ax, 'HD-EMG', 'FontName', 'Tahoma', 'Color', [0.6 0.6 0.6]);
    ylabel(ax, 'Channel', 'FontName','Tahoma','Color','k');
    xlabel(ax, 'Time (s)', 'FontName','Tahoma','Color','k');
    
    if AUTO_SAVE_FIGURES
        default.savefig(fig, fullfile("R:/NMLShare/generated_data/human/TC/MUAPs", block, string(dt), ...
            sprintf("MUAP-Rates_N-%d_THRESH-%d_SMOOTH-%d_MAX-MUAPS-%d", N_SAMPLES_PP, round(THRESH_RMS*10), round(SMOOTHING*10), round(MAX_MUAPS_PER_SEC))));
    end
    %% Now, decompose those rates
    [coeff,score,latent,tsquared,explained,mu] = pca(Yp);
    fig = figure('Name', 'Rate PCs Pareto', 'Color', 'w', ...
        'Position', [4139 -838 1405 420]);
    
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
    [~, iMax] = max(S, [], 2);
    S = S ./ max(S, [], 1);
    for ii = 1:6
        scatter(ax, tp(iMax == ii), iMax(iMax == ii), 'Marker', 's', 'MarkerFaceColor', 'flat'); 
    end
    title(ax, 'Motor Decomp "State"', 'FontName', 'Tahoma', 'Color', 'k');
    ylim(ax, [0 7]);
    xlabel(ax, 'Time (s)', 'FontName', 'Tahoma', 'Color', 'k');
    if AUTO_SAVE_FIGURES
        default.savefig(fig, fullfile("R:/NMLShare/generated_data/human/TC/MUAPs", block, string(dt),  "MUAP-Rates-Pareto-State"));
    end
    
    %% Spatial distribution of PCs on electrodes
    coeff_big = nan(128,6);
    for ii = 1:6
        coeff_big(in_both, ii) = coeff(1:(numel(in_both)),ii);
    %     coeff_big(in_A, ii) = coeff(1:(numel(in_A)),ii);
    %     coeff_big(in_B, ii) = coeff((numel(in_A)+1):end, ii);
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
        default.savefig(fig, fullfile("R:/NMLShare/generated_data/human/TC/MUAPs", block, string(dt), "MUAP-Rates-Components"));
    end
end