%% load dataset
% load('dataset_201.mat', 'df');      % load data of 201 brain regions
load('dataset_53.mat', 'df');       % or load data of 53 brain regions
% load('dataset_64.mat', 'df');       % or load data of 64 brain regions
load('onto.mat', 'onto');

ivol = cat(1, df.Volume{:});        % get volume of each region

% convert table to array
ds = df(:, {'HomeCage', 'ISO', 'Saline', 'Ket'});
ds = mat2cell(table2cell(ds), size(ds, 1), ones(size(ds, 2), 1));
ds = cellfun(@(x)(cat(1, x{:})), ds, 'UniformOutput', false);

nsub = cellfun(@(x)(size(x, 2)), ds);       % get N of each group
grp_sub = mat2cell(1:sum(nsub), 1, nsub);

da_all = cat(2, ds{:});     % concatenate 4 groups

% normalization
da_norm = log10(da_all+1);

iso_ctrl_ave = mean(da_norm(:, grp_sub{1}), 2);
iso_ctrl_var = std(da_norm(:, grp_sub{1}), 1, 2);
da_norm(:, [grp_sub{1}, grp_sub{2}]) = ...
    (da_norm(:, [grp_sub{1}, grp_sub{2}]) - iso_ctrl_ave) ./ iso_ctrl_var;
da_norm(iso_ctrl_var==0, [grp_sub{1}, grp_sub{2}]) = nan;

ket_ctrl_ave = mean(da_norm(:, grp_sub{3}), 2);
ket_ctrl_var = std(da_norm(:, grp_sub{3}), 1, 2);
da_norm(:, [grp_sub{3}, grp_sub{4}]) = ...
    (da_norm(:, [grp_sub{3}, grp_sub{4}]) - ket_ctrl_ave) ./ ket_ctrl_var;
da_norm(ket_ctrl_var==0, [grp_sub{3}, grp_sub{4}]) = nan;

%% p values
temp = da_norm(:, grp_sub{2});
nr = size(da_norm, 1);
[hsw, psr] = deal(nan(nr, 1));
for nn = 1:nr
    try
        hsw(nn) = swtest(temp(nn, :));
        psr(nn) = signrank(temp(nn, :));
    catch
        continue;
    end
end
[~, ptt] = ttest(temp');
ptt_adj = fdr_bh(ptt, 'pdep');
psr_adj = fdr_bh(psr, 'pdep');
[hsw, ptt(:), ptt_adj(:), psr, psr_adj];

%% plot parameters
clr = [225, 133, 97; 100, 176, 150]/255;
rname = cellfun(@(x)(sprintf('%s (%s)', onto.name{strcmp(onto{:, 2}, x)}, x)), ...
    df.Name, 'UniformOutput', false);

%% plot line chart for 53 areas
figure; set(gcf, 'Position', [100, 100, 1080, 760], ...
    'PaperOrientation', 'landscape'); 
subplot('Position', [0.05, 0.50, 0.9, 0.38]); hold on;

temp = da_norm(:, grp_sub{2});
[nr, nsub] = size(temp);
errorbar(1:nr, mean(temp, 2), std(temp, 0, 2)/sqrt(nsub), 'o-', ...
    'CapSize', 5, 'LineWidth', 0.7, 'MarkerSize', 6, ...
    'MarkerFaceColor', clr(1, :), 'Color', clr(1, :));

temp = da_norm(:, grp_sub{4});
[nr, nsub] = size(temp);
errorbar(1:nr, mean(temp, 2), std(temp, 0, 2)/sqrt(nsub), 'o-', ...
    'CapSize', 5, 'LineWidth', 0.7, 'MarkerSize', 6, ...
    'MarkerFaceColor', clr(2, :), 'Color', clr(2, :));

grid on;
set(gca, 'TickDir', 'out', 'TickLength', [0.005, 0.1], 'FontSize', 10, ...
    'XLim', [0, nr+1], 'YLim', [-2, 10], 'XTick', 1:nr, 'XTickLabel', ...
    rname, 'YTickLabelRotation', 90, 'XTickLabelRotation', 90);
legend({'ISO', 'KET'}, 'Box', 'off', 'Location', 'best');
ylabel('Z-score');

%% plot line chart for 201 areas
% part 1
figure; set(gcf, 'Position', [100, 100, 1080, 760], 'PaperOrientation', 'landscape'); 
subplot('Position', [0.05, 0.50, 0.9, 0.38]); hold on;

temp = da_norm(:, grp_sub{2});
[nr, nsub] = size(temp);
errorbar(1:nr, mean(temp, 2), std(temp, 0, 2)/sqrt(nsub), 'o-', ...
    'CapSize', 3, 'LineWidth', 0.7, 'MarkerSize', 4, ...
    'MarkerFaceColor', clr(1, :), 'Color', clr(1, :));

temp = da_norm(:, grp_sub{4});
[nr, nsub] = size(temp);
errorbar(1:nr, mean(temp, 2), std(temp, 0, 2)/sqrt(nsub), 'o-', ...
    'CapSize', 3, 'LineWidth', 0.7, 'MarkerSize', 4, ...
    'MarkerFaceColor', clr(2, :), 'Color', clr(2, :));

grid on;
set(gca, 'TickDir', 'out', 'TickLength', [0.003, 0.1], 'FontSize', 8, ...
    'XLim', [0.5, 101.5], 'YLim', [-4, 19], 'XTick', 1:nr, 'XTickLabel', ...
    rname, 'YTickLabelRotation', 90);
legend({'ISO', 'KET'}, 'Box', 'off', 'Location', 'best');
ylabel('Z-score');

% part 2
figure; set(gcf, 'Position', [100, 100, 1080, 760], 'PaperOrientation', 'landscape'); 
subplot('Position', [0.05, 0.50, 0.9, 0.38]); hold on;

temp = da_norm(:, grp_sub{2});
[nr, nsub] = size(temp);
errorbar(1:nr, mean(temp, 2), std(temp, 0, 2)/sqrt(nsub), 'o-', ...
    'CapSize', 3, 'LineWidth', 0.7, 'MarkerSize', 4, ...
    'MarkerFaceColor', clr(1, :), 'Color', clr(1, :));

temp = da_norm(:, grp_sub{4});
[nr, nsub] = size(temp);
errorbar(1:nr, mean(temp, 2), std(temp, 0, 2)/sqrt(nsub), 'o-', ...
    'CapSize', 3, 'LineWidth', 0.7, 'MarkerSize', 4, ...
    'MarkerFaceColor', clr(2, :), 'Color', clr(2, :));

grid on;
set(gca, 'TickDir', 'out', 'TickLength', [0.003, 0.1], 'FontSize', 8, ...
    'XLim', [101.5, 202], 'YLim', [-4, 19], 'XTick', 1:nr, 'XTickLabel', ...
    rname, 'YTickLabelRotation', 90);
legend({'ISO', 'KET'}, 'Box', 'off', 'Location', 'best');
ylabel('Z-score');
