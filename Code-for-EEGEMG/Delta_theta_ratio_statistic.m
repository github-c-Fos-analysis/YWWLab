% load path file %Before and affter drug administration % delta/theta ratio
% differance between two groups.
% nex_idx = 7;
mat_idx = 7;
t1_idx = 2;
t2_idx = 4;
t3_idx = 5;
t4_idx = 6;
iso_idx = 4;%drug administration
%
summary1 = [];
summary2 = [];
summary3 = [];
summary4 = [];
toplot1 = [];
toplot2 = [];
toplot3 = [];
toplot4 = [];

for i = 1:2
    for fn = 1:size(fn_idx1,1) %change name of the data
        %load
        eval(['load(fn_idx',num2str(i),'{fn,mat_idx});']);%Merge characters
        dff = 25/size(recordingFile.eeg,1);
        idx_theta = round([6 10]/dff);
        idx_delta = round([0.5 4]/dff);
        eegRatio = sum(recordingFile.eeg(idx_delta(1):idx_delta(2),:),1)./...
            sum(recordingFile.eeg(idx_theta(1):idx_theta(2),:),1);
        eval(['win = [fn_idx',num2str(i),'{fn,t1_idx}:fn_idx',num2str(i),'{fn,t4_idx}];']);
        tmp = smooth(eegRatio,50);
        tmp = tmp(win);
        tmp = tmp';
        eval(['t = fn_idx',num2str(i),'{fn,iso_idx};']);
        %%temp = cat(2,mean(eegRatio(t-500:t)),mean(eegRatio(t:t+510)));
            temp = cat(2,mean(eegRatio(t-300:t)),mean(eegRatio(t:t+310)));
        eval(['summary',num2str(i),' = cat(1,summary',num2str(i),',temp);']);
        eval(['toplot',num2str(i),'= cat(1,toplot',num2str(i),',tmp);']);
    end
end

toplot = [];
for i = 1:2
    eval(['toplot = cat(1,toplot,mean(toplot',num2str(i),',1));']);
end
figure
plot(toplot')
legend('CTRL','EXP');
%summary1/2 means CTRL/EXP; (:,1)means 20min before drug
%administration;(:,2) means 20min after durg administration;
[H,P,~] = swtest(summary2(:,2)) %H=0 means abnormal distribution
[P,H] = ranksum(summary1(:,2),summary2(:,2))%summary1 CTRL, summary2 EXP
