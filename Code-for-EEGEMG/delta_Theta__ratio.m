% load path file
% nex_idx = 8;

%
mat_idx = 7;
t1_idx = 2;
t2_idx = 4;
t3_idx = 5;
t4_idx = 6;
%
summary_co = [];
summary_ex = [];
toplot_co = [];
toplot_ex = [];
for fn = 1:size(fn_idx1,1)
    %load 
    load(fn_idx1{fn,mat_idx});
    dff = 25/size(recordingFile.eeg,1);
    idx_theta = round([6 10]/dff);
    idx_delta = round([0.5 4]/dff);
    eegRatio = sum(recordingFile.eeg(idx_delta(1):idx_delta(2),:),1)./...
        sum(recordingFile.eeg(idx_theta(1):idx_theta(2),:),1);
    figure
    win = [fn_idx1{fn,t1_idx}:fn_idx1{fn,t4_idx}];
    tmp = smooth(eegRatio,50);
    tmp = tmp(win);
    plot(win,tmp),hold on
    plot([fn_idx1{fn,t2_idx} fn_idx1{fn,t2_idx}],[0 max(tmp)])
    plot([fn_idx1{fn,t3_idx} fn_idx1{fn,t3_idx}],[0 max(tmp)])
    sum1 = mean(eegRatio(fn_idx1{fn,t1_idx}:fn_idx1{fn,t2_idx}));
    sum2 = mean(eegRatio(fn_idx1{fn,t2_idx}:fn_idx1{fn,t3_idx}));
    sum3 = mean(eegRatio(fn_idx1{fn,t3_idx}:fn_idx1{fn,t4_idx}));
    title([num2str(sum1),'  ',num2str(sum2),'  ',num2str(sum3)])
    hold off
    temp = [sum1 sum2 sum3];
    summary_co = cat(1,summary_co,temp);
    toplot_co = cat(1,toplot_co,tmp');
end
% colormap(lines(size(fn_idx,1)))
for fn = 1:size(fn_idx2,1)
    %load 
    load(fn_idx2{fn,mat_idx});
    dff = 25/size(recordingFile.eeg,1);
    idx_theta = round([6 10]/dff);
    idx_delta = round([0.5 4]/dff);
    eegRatio = sum(recordingFile.eeg(idx_delta(1):idx_delta(2),:),1)./...
        sum(recordingFile.eeg(idx_theta(1):idx_theta(2),:),1);
    figure
    win = [fn_idx2{fn,t1_idx}:fn_idx2{fn,t4_idx}];
    tmp = smooth(eegRatio,50);
    tmp = tmp(win);
    plot(win,tmp),hold on
    plot([fn_idx2{fn,t2_idx} fn_idx2{fn,t2_idx}],[0 max(tmp)])
    plot([fn_idx2{fn,t3_idx} fn_idx2{fn,t3_idx}],[0 max(tmp)])
    sum1 = mean(eegRatio(fn_idx2{fn,t1_idx}:fn_idx2{fn,t2_idx}));
    sum2 = mean(eegRatio(fn_idx2{fn,t2_idx}:fn_idx2{fn,t3_idx}));
    sum3 = mean(eegRatio(fn_idx2{fn,t3_idx}:fn_idx2{fn,t4_idx}));
    title([num2str(sum1),'  ',num2str(sum2),'  ',num2str(sum3)])
    hold off
    temp = [sum1 sum2 sum3];
    summary_ex = cat(1,summary_ex,temp);
    toplot_ex = cat(1,toplot_ex,tmp');
end
figure
plot(toplot_co','color',[0.5 1 0.83]),hold on
plot(toplot_ex','color',[1 0.75 0.8])
h1 = plot(mean(toplot_co,1),'g','lineWidth',2)
h2 = plot(mean(toplot_ex,1),'r','lineWidth',2)
legend([h1 h2],'CTRL','EXP')
axis tight
hold off