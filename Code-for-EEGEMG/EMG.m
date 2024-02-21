clear
clc
%%
mat_idx = 7;
t1_idx = 2;
t2_idx = 3;
t3_idx = 4;
t4_idx = 5;
%
toplot_co = [];
for fn = 1:size(fn_idx1,1)%读取行数逗号后的1表示行，列用2表示
    %load 
    load(fn_idx1{fn,mat_idx});
    tmp = sum(recordingFile.emg,1);
    figure
    win = [fn_idx1{fn,t1_idx}:fn_idx1{fn,t4_idx}];
%     tmp = smooth(eegRatio,50);
    tmp = tmp(win);
    plot(win,tmp),hold on
    plot([fn_idx1{fn,t2_idx} fn_idx1{fn,t2_idx}],[0 max(tmp)],'r--','lineWidth',1.5)
    plot([fn_idx1{fn,t3_idx} fn_idx1{fn,t3_idx}],[0 max(tmp)],'r--','lineWidth',1.5)
    hold off
    toplot_co = cat(1,toplot_co,tmp);
end
%%
toplot_ex = [];
for fn = 1:size(fn_idx2,1)
    %load 
    load(fn_idx2{fn,mat_idx});
    tmp = sum(recordingFile.emg,1);
    figure
    win = [fn_idx2{fn,t1_idx}:fn_idx2{fn,t4_idx}];
%     tmp = smooth(eegRatio,50);
    tmp = tmp(win);
    plot(win,tmp),hold on
    plot([fn_idx2{fn,t2_idx} fn_idx2{fn,t2_idx}],[0 max(tmp)],'r--','lineWidth',1.5)
    plot([fn_idx2{fn,t3_idx} fn_idx2{fn,t3_idx}],[0 max(tmp)],'r--','lineWidth',1.5)
    hold off
    toplot_ex = cat(1,toplot_ex,tmp);
end
%%
figure
plot(toplot_ex','color',[1 0.75 0.8]),hold on
h1 = plot(mean(toplot_ex,1),'r')
plot(-toplot_co','color',[0.5 1 0.83])
h2 = plot(-mean(toplot_co,1),'g')
hold off
axis tight
legend([h1 h2],'EXP','CTRL')