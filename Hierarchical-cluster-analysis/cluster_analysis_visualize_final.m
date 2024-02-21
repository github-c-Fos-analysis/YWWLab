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
da_norm = (da_all+1)./ivol;  

iso_ctrl_ave = mean(da_norm(:, grp_sub{1}), 2);
da_norm(:, [grp_sub{1}, grp_sub{2}]) = ...
    da_norm(:, [grp_sub{1}, grp_sub{2}]) ./ iso_ctrl_ave;

ket_ctrl_ave = mean(da_norm(:, grp_sub{3}), 2);
da_norm(:, [grp_sub{3}, grp_sub{4}]) = ...
    da_norm(:, [grp_sub{3}, grp_sub{4}]) ./ ket_ctrl_ave;

da_norm = log10(da_norm);

%%
npc = 2;
hr = 0.6;
nr = size(df, 1);
da = da_norm(:, grp_sub{4});

[~, score] = pca(da, 'Centered', true);
Z = linkage(score(:, 1:npc), 'complete', 'euclidean');
C = cluster(Z, 'Cutoff', hr * max(Z(:, 3)), 'criterion', 'distance');
[~, nid, cout] = dendrogram(Z, nr);
[~, oo] = ismember(nid, cout);
[C, oo];
% [~, oo] = sort(oo);
% C = C(oo);

%% optimize clustering parameters
cri_name = {'Number of Clusters', 'Silhouette'};
grp_name = {'Home Cage', 'Isoflurane', 'Saline', 'Ketamine'};
grp_idx = [1, 2, 3, 4];
npc = 2;
hratio = (0.3:0.05:1);

% calculate cluster metrics
ngrp = length(grp_idx);
nh = length(hratio);
cri = nan(ngrp, nh, 2);
for gg = 1:ngrp
    da = da_norm(:, grp_sub{grp_idx(gg)});
    [~, score] = pca(da, 'Centered', true);
    X = score(:, 1:npc);
    Z = linkage(X, 'complete', 'euclidean');
    for kk = 1:nh
        C = cluster(Z, 'Cutoff', hratio(kk) * max(Z(:, 3)), 'criterion', 'distance');
        cri(gg, kk, 1) = max(C);
        eva = evalclusters(X, C, 'silhouette', 'Distance', 'Euclidean');
        cri(gg, kk, 2) = eva.CriterionValues;
    end
end

% plot
for cc = 1:length(cri_name)
    figure; set(gcf, 'Position', [100, 100, 480, 360]); hold on; 
    for gg = 1:ngrp
        plot(hratio, cri(gg, :, cc), '-o', 'LineWidth', 1, 'MarkerSize', 5);
    end
    xlabel('Ratio of Tree Height', 'FontSize', 10);
    ylabel(cri_name{cc}, 'FontSize', 10);
    set(gca, 'TickDir', 'out', 'LineWidth', 1.3, 'XLim', [0.25, 1.01]);
    legend(grp_name, 'Box', 'off', 'FontSize', 10, 'Location', 'bestoutside');
    title(cri_name{cc});
end

%% cluster distance matrix
grp_idx = [2, 4];
grp_name = {'Isoflurane', 'Ketamine'};
hratio = 0.5;
npc = 2;

nr = size(df, 1);
for gg = 1:length(grp_idx)
    da = da_norm(:, grp_sub{grp_idx(gg)});
    [~, score] = pca(da, 'Centered', true);
    Z = linkage(score(:, 1:npc), 'complete', 'euclidean');
    Cs = cluster(Z, 'Cutoff', hratio * max(Z(:, 3)), 'criterion', 'distance');
           
    for kk = 1:length(hratio)   
        figure; set(gcf, 'Position', [0, 0, 800, 620], 'PaperOrientation', 'landscape');

        subplot('Position', [0.1, 0.795, 0.6, 0.175]);
        [hs, nid, cout] = dendrogram(Z, nr); axis off; 
        set(hs, 'Color', 'k');
        title(sprintf('%s ratio: %.1f', grp_name{gg}, hratio(kk)), 'FontSize', 10);
        [~, oo] = ismember(nid, cout);
        [~, oo] = sort(oo);
        C = Cs(oo, kk);
        h = gca;
        
        subplot('Position', [0.1, 0.775, 0.6, 0.02]); hold on;
        for nn = 1:nr
            ic = df.Color{oo(nn)} / 255;
            patch([-0.5, 0.5, 0.5, -0.5, -0.5]+nn, [0, 0, 1, 1, 0], ...
                ic, 'EdgeColor', ic);
        end
        set(gca, 'XLim', get(h, 'XLim'), 'YLim', [0, 1.05]);
        axis off;

        subplot('Position', [0.1, 0.00, 0.6, 0.775]); hold on;
        dm = squareform(pdist(da(oo, :), 'euclidean'));
        imagesc(1:nr, 1:nr, dm, [0, 5]); colormap(Spectral);
        for cc = 1:max(C)
            c1 = find(C == cc, 1);
            c2 = find(C == cc, 1, 'last');
            plot([c1-0.5, c1-0.5, c2+0.5, c2+0.5, c1-0.5], ...
                [c1-0.5, c2+0.5, c2+0.5, c1-0.5, c1-0.5], ...
                'Color', 'k', 'LineWidth', 1.5);
        end
        for nn = 1:nr
            rname = onto.name{strcmp(onto{:, 2}, df.Name{oo(nn)})};
            text(nr+1, nn, sprintf('(%s) %s', df.Name{oo(nn)}, rname), ...
                'Rotation', 0, 'FontWeight', 'normal', ...
                'FontSize', 7, 'HorizontalAlignment', 'left', 'FontName', 'Arial');
        end
        set(gca, 'XLim', get(h, 'XLim'), 'YLim', get(h, 'XLim'), 'YDir', 'reverse');
        axis off; axis equal;
            
%         subplot('Position', [0, 0.00, 1, 0.775]); hold on;
%         dm = squareform(pdist(da(oo, :), 'euclidean'));
%         imagesc(1:nr, 1:nr, dm, [0, 5]); colormap(Spectral);
%         c = colorbar; 
%         c.Position = [0.1, 0.8, 0.8, 0.15];
%         c.Location = 'North';
%         c.Box = 'off';
%         c.TickDirection = 'out';
%         c.Ticks = (0:1:5);
%         c.LineWidth = 1.3;
%         c.FontSize = 10;
%         c.Label.String = 'Euclidean Distance';
    end
end

%% cluster Z-value hotmap
npc = 2;
hr = 0.5;

figure; set(gcf, 'Position', [100, 100, 650, 800]);

% Isoflurane
da = da_norm(:, grp_sub{2});
[~, score] = pca(da, 'Centered', true);
Z = linkage(score(:, 1:npc), 'complete', 'euclidean');
[nr, nsub] = size(da);
da = da ./ (std(da, 1, 1) / sqrt(nsub));

subplot('Position', [0.32, 0.02, 0.08, 0.85]);
[hs, nid, cout] = dendrogram(Z, nr, 'Orientation', 'left'); axis off;
set(hs, 'Color', 'k');
[~, oo] = ismember(nid, cout);
[~, oo] = sort(oo);
h = gca;

subplot('Position', [0.40, 0.02, 0.02, 0.85]); hold on;
for nn = 1:nr
    ic = df.Color{oo(nn)} / 255;
    rname = df.Name{oo(nn)};
    patch([0, 0, 1, 1, 0], [-0.5, 0.5, 0.5, -0.5, -0.5]+nn, ...
        ic, 'EdgeColor', ic);
    text(0.99, nn, rname, 'Rotation', 0, 'FontWeight', 'bold', ...
        'FontSize', 2, 'HorizontalAlignment', 'right', 'FontName', 'Arial');
end
set(gca, 'YLim', get(h, 'YLim'), 'XLim', [-0.05, 1.05]);
axis off;

subplot('Position', [0.42, 0.02, 0.05, 0.85]); hold on;
imagesc(1:nsub, 1:nr, da(oo, :), [-10, 10]); colormap(parula);
rr = find(mean(da(oo, :), 2) > 2);
r1 = rr(diff([-inf; rr]) > 1)';
r2 = rr(diff([rr; inf]) > 1)';
plot(reshape(repmat([0.5; nsub+0.5; nsub+0.5; 0.5; 0.5; nan], length(r1), 1), [], 1), ...
    reshape([-0.5+r1; -0.5+r1; 0.5+r2; 0.5+r2; -0.5+r1; nan+r2], [], 1), '-r');
rr = find(mean(da(oo, :), 2) < -2);
r1 = rr(diff([-inf; rr]) > 1)';
r2 = rr(diff([rr; inf]) > 1)';
plot(reshape(repmat([0.5; nsub+0.5; nsub+0.5; 0.5; 0.5; nan], length(r1), 1), [], 1), ...
    reshape([-0.5+r1; -0.5+r1; 0.5+r2; 0.5+r2; -0.5+r1; nan+r2], [], 1), '-b');
set(gca, 'XLim', [0.5, nsub+0.5], 'YLim', get(h, 'YLim'), 'YDir', 'normal');
axis off;
title('Iso');

% Ketamine
da = da_norm(:, grp_sub{4});
[~, score] = pca(da, 'Centered', true);
Z = linkage(score(:, 1:npc), 'complete', 'euclidean');
C = cluster(Z, 'Cutoff', hr * max(Z(:, 3)), 'criterion', 'distance');
[nr, nsub] = size(da);
da = da ./ (std(da, 1, 1) / sqrt(nsub));

subplot('Position', [0.60, 0.02, 0.08, 0.85]);
[hs, nid, cout] = dendrogram(Z, nr, 'Orientation', 'right'); axis off;
set(hs, 'Color', 'k');
[~, oo] = ismember(nid, cout);
[~, oo] = sort(oo);
h = gca;

subplot('Position', [0.58, 0.02, 0.02, 0.85]); hold on;
for nn = 1:nr
    ic = df.Color{oo(nn)} / 255;
    rname = df.Name{oo(nn)};
    patch([0, 0, 1, 1, 0], [-0.5, 0.5, 0.5, -0.5, -0.5]+nn, ...
        ic, 'EdgeColor', ic);
    text(0.01, nn, rname, 'Rotation', 0, 'FontWeight', 'bold', ...
        'FontSize', 2, 'HorizontalAlignment', 'left', 'FontName', 'Arial');
end
set(gca, 'YLim', get(h, 'YLim'), 'XLim', [-0.05, 1.05]);
axis off;

subplot('Position', [0.53, 0.02, 0.05, 0.85]); hold on;
imagesc(1:nsub, 1:nr, da(oo, :), [-10, 10]); colormap(parula);
rr = find(mean(da(oo, :), 2) > 2);
r1 = rr(diff([-inf; rr]) > 1)';
r2 = rr(diff([rr; inf]) > 1)';
plot(reshape(repmat([0.5; nsub+0.5; nsub+0.5; 0.5; 0.5; nan], length(r1), 1), [], 1), ...
    reshape([-0.5+r1; -0.5+r1; 0.5+r2; 0.5+r2; -0.5+r1; nan+r2], [], 1), '-r');
rr = find(mean(da(oo, :), 2) < -2);
r1 = rr(diff([-inf; rr]) > 1)';
r2 = rr(diff([rr; inf]) > 1)';
plot(reshape(repmat([0.5; nsub+0.5; nsub+0.5; 0.5; 0.5; nan], length(r1), 1), [], 1), ...
    reshape([-0.5+r1; -0.5+r1; 0.5+r2; 0.5+r2; -0.5+r1; nan+r2], [], 1), '-b');
set(gca, 'XLim', [0.5, nsub+0.5], 'YLim', get(h, 'YLim'), 'YDir', 'normal');
axis off;
title('Ket');

c = colorbar;
c.Position = [0.85, 0.40, 0.02, 0.20];
c.Box = 'off';
c.TickDirection = 'out';
c.LineWidth = 1;
c.FontSize = 10;
c.Label.String = 'Z value';
c.Label.Rotation = 270;
c.Label.VerticalAlignment = 'bottom';

