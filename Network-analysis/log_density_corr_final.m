%% load dataset
% load('dataset_201.mat', 'df');      % load data of 201 brain regions
% load('dataset_53.mat', 'df');       % or load data of 53 brain regions
load('dataset_64.mat', 'df');       % or load data of 64 brain regions
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
da_norm = log10(da_norm);

%% correlation matrix hotmap for 64 areas
grp_name = {'Home Cage', 'Isoflurane', 'Saline', 'Ketamine'};
grp_idx = [1, 2, 3, 4];

nr = size(df, 1);
for gg = 1:length(grp_idx)
    % calculate correlation
    da = 1-squareform(pdist(da_norm(:, grp_sub{gg}), 'correlation'));
    % sort each group by density distance
    oos = cell(6, 1);
    for cc = 1:1:6
        inc = find(cat(1, df.GroupID{:}) == cc);
        Z = linkage(da(inc, :), 'complete', 'euclidean');
        [~, nid, cout] = dendrogram(Z, length(inc)); close;
        [~, oo] = ismember(nid, cout);
        [~, oo] = sort(oo);
        oos{cc} = inc(oo);
    end
    oo = cat(1, oos{:});

    figure; set(gcf, 'Position', [100, 100, 800, 650], 'PaperOrientation', 'landscape');

    subplot('Position', [0.05, 0.10, 0.05, 0.80]); hold on;
    for nn = 1:nr
        ic = df.Color{oo(nn)} / 255;
        patch([0, 0, 1, 1, 0], [-0.5, 0.5, 0.5, -0.5, -0.5]+nn, ...
            ic, 'EdgeColor', ic);
        text(0.95, nn, df.Name{oo(nn)}, 'Rotation', 0, 'FontWeight', 'bold', ...
            'FontSize', 8, 'HorizontalAlignment', 'right', 'FontName', 'Arial');
    end
    set(gca, 'YLim', [0.5, nr+0.5], 'XLim', [-0.05, 1.05], 'YDir', 'reverse');
    axis off;

    subplot('Position', [0.10, 0.10, 0.65, 0.80]); hold on;
    imagesc(1:nr, 1:nr, da(oo, oo), [-1, 1]); colormap(flipud(RdYlBu));
    set(gca, 'XLim', [0.5, nr+0.5], 'YLim', [0.5, nr+0.5], 'YDir', 'reverse');
    title(sprintf('%s', grp_name{grp_idx(gg)}), 'FontSize', 10);
    axis off; axis equal;

    c = colorbar;
    c.Position = [0.80, 0.15, 0.03, 0.70];
    c.Box = 'off';
    c.TickDirection = 'out';
    c.LineWidth = 1.3;
    c.FontSize = 10;
    c.Label.String = 'Correlation';
    c.Label.Rotation = 270;
    c.Label.VerticalAlignment = 'bottom';

%     tbl = array2table(da(oo, oo), 'VariableNames', df.Name(oo), 'RowNames', df.Name(oo));
%     writetable(tbl, 'log_density_corr.xlsx', 'Sheet', grp_name{grp_idx  (gg)}, ...
%         'WriteRowNames', true);
end


