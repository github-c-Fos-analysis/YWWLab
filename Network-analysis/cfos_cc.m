% cc heatmap and nodes for cytoscape(P<0.05,cc>0)
clc
clear
[~,~,raw] = xlsread('/Users/luna/Downloads/data.xlsx');
data = brainregionISOKeta(:,3:2:57);  
data= cell2mat(data); 
ISOmean=mean(data,2)
STD=std(data,0,2)

dataI = cell2mat(dataI); 
dataK= cell2mat(dataK);
cc = [];
P_value = [];
[cc,P_value] = corrcoef(data');
[cc,P_value] = corrcoef(dataK');

toclear = find(isnan(cc(:,1)));
cc(toclear,:)=[];
cc(:,toclear)=[];
P_value(toclear,:)=[];
P_value(:,toclear)=[]; 
%% 
CTXK=raw(1:19,13:18);
THK = raw(27:31,13:18);

CTXK = cell2mat(CTXK); 
THK= cell2mat(THK)

CTXI=raw(1:19,27:32);
THI = raw(27:31,27:32);

CTXI = cell2mat(CTXI); 
THI= cell2mat(THI);

Icc= [];
IP_value = [];

Kcc= [];
KP_value = [];

[Icc,IP_value] = corr(dataI');
[Kcc,KP_value] = corrcoef(dataK');

%% 

ZK= zscore(dataK,0,2);
ZI=zscore(dataI,0,2);


H=[]
H2=[]
P=[]

for i=1:19 
[h1,P1,~]=swtest(ZK(i,:)); 
[h2,P2,~]=swtest(ZI(i,:));
H(i,1)=h1
H(i,2)=h2
if h1==0 & h2==0 
[h,p,~]=ttest(ZK(i,:),ZI(i,:));

else
[p,h]=signrank(ZK(i,:),ZI(i,:));
end
P=cat(1,P,p);
mK=mean(ZK,2);
mI=mean(ZI,2)
scatter(mK,mI)
xlabel("Ketamine")
ylabel("ISO")
text(mk,mI,raw(:,1))
end


%%  Mean r analysis

[ISEM,Imean,Int]=Get_ms(Icc1)
[KSEM,Kmean,Knt]=Get_ms(Kcc1)
p = compare_correlation_coefficients(mean,r1,63,62)

figure
imagesc(Kcc)
load('/Users/luna/Desktop/Script/myColormap.mat')
colormap(myColormap)
colorbar


[row,col]=find(P_value<0.05);
len=length(row)
row1=[]
col1=[]
for i=1:len
    if cc(row(i),col(i))>0.82
      row1=[row1;row(i)] 
      col1=[col1;col(i)]
    end
end
%
len1=length(row1)
result=[]
for i=1:len1
    t=[raw(col1(i),1) raw(row1(i),1) cc(row1(i),col1(i)) raw(col1(i),2)]
    result=[result;t]
end
%% 

isnorm = 1;c
if isnorm
    data = data(:,1:2:end);
else
    data = data(:,2:2:end);
end


%%
[Ci,Q]=modularity_und(cc,1)
Z=module_degree_zscore(cc,Ci,0)
P=participation_coef(cc,Ci,0)
[deg] = degrees_und(cc)
BC=betweenness_bin(cc)

disp(char(result(:,1)))
disp(char(result(:,2)))
%%
figure
imagesc(cc)
load('/Users/Desktop/c-Fos20211221/HY cfos/myColormap.mat')
colormap(myColormap)
colorbar
tic
cgObj1 = clustergram(cc,'ColumnPDist','euclidean');
cgObj1.RowLabels = raw(str2num(cell2mat(cgObj1.RowLabels)),1);
cgObj1.ColumnLabels = raw(str2num(cell2mat(cgObj1.RowLabels)),1);

toc
set(gcf,'position',[10 10 1200 900])
%% 
n = 64
M = ones(n)
for i = 1:n
    for j = 1:n
        if P_value(i,j) < 0.05 & cc(i,j) > 0
            M(i,j) = cc(i,j)
        end
    end
end

[row,col]=find(M);
len=length(row)
result=[]
for i=1:len
    t=[raw(col(i),1) raw(row(i),1)  M(row(i),col(i)) raw(col(i),2)]
    result=[result;t]
end

