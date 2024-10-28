clear all; clc;
%% 請設定以下資料
Pup=[];          %%% 輸入 上游Protein 名字 N * 1 個 (column vector)
Gup=[];        %%% 輸入 上游基因 & microRNA 名字 M * 1 個 (column vector)
Gdn=[];          %%% 輸入 下游基因 & microRNA 名字 M * 1 個 (colum vector)
Appi=[];        %%% 輸入 PPI矩陣 應為對稱方陣 維度 N * N
Agrn=[];        %%% 輸入 GRN矩陣 非對稱方陣   維度 M (target Gene) * M (TF)

%% load data 
load('PPI_PNP.mat'); PPI_ID_normal = PPI_PNP;
load('GRN_PNP.mat'); GRN_ID_normal = GRN_PNP;
load('GRN_name.mat')
load('PPI_name.mat')
GRNlate_0_GRNlate = GRN_ID_normal;
PPIID_late_0 = PPI_ID_normal;
diff = size(GRN_ID_normal,1) - size(PPI_ID_normal,1);
% append zeros matrix
PPIID_late_0_p1=[PPIID_late_0(:,:) zeros(size(PPIID_late_0,1),diff) ];
PPIID_late_0_p2=[PPIID_late_0_p1(:,:); zeros(diff,size(PPIID_late_0_p1,2))];
PPIID_late_0=PPIID_late_0_p2; 
clear PPIID_late_0_p1 PPIID_late_0_p2
%% 1. 建A矩陣
% [Gene;mRNA;lncRNA]
GN=[]; PN=[]; Appi=[]; Agrn=[];
huname_GRN_unique = GRN_name;
GRN_ID_normal_genename = GRN_name;

UPindG=[]; UPindG=1:1:length(huname_GRN_unique); UPindM=find(strncmpi(huname_GRN_unique,'MIR',3));UPindG(UPindM)=[]; 
% concatenate microRNA name below gene name
Gup_f=cat(1,huname_GRN_unique(UPindG),huname_GRN_unique(UPindM));  PN=GRN_ID_normal_genename; GN=Gup_f; %clear huname_GRN_unique_f

Appi=PPIID_late_0;
Agrn_f=cat(2,GRNlate_0_GRNlate(:,UPindG),GRNlate_0_GRNlate(:,UPindM)); Agrn=cat(1,Agrn_f(UPindG,:),Agrn_f(UPindM,:)); %clear Agrn_f GRNlate_0_GRNlate
% sum Appi and Agrn in row direction and exclude the values with 0,respectively
ind=[];ind=find(sum(Appi~=0,2)==0); PN(ind)=[];Appi(ind,:)=[];%upstream:protein
ind=[];ind=find(sum(Agrn~=0,2)==0); GN(ind)=[];Agrn(ind,:)=[];%upstream:gene,linc,miR  

Gup=Gup_f;Gdn=[PN;GN]; A=[]; A=[Appi;Agrn]; %A=interaction matrix
ind=[]; ind=find(sum(A~=0,1)==0); Gup(ind)=[]; A(:,ind)=[]; %downstream

%% Principal Network Projection(PNP): sinhuname_uniquelar value decomposition (U,V)=sinhuname_uniquelar vector;(S)=sinhuname_uniquelar value
%  To see if setting an upper/lower bound
fprintf('Start PNP')
A=full(A); tmptmp=[]; for i=1:size(A,1);tmptmp=[tmptmp;A(i,find(A(i,:)~=0))'];end; LP33=sort(tmptmp,'descend');
lower=-5*10^(2); upper=5*10^(2);
A=A.*(A>lower)+(lower).*(A<=lower); A=A.*(A<upper)+(upper).*(A>=upper);

% 1. Acquire sinhuname_uniquelar vectors that represents 85% of the real network for projection 
[U1,S1,V1]=svd((A)); % U1=mxm S1=mxn V1=nxn 
S1d=diag(S1); % Sld=nx1
for i=1:length(S1d)
    S1acc(i)=sum((S1d(1:i).^2)./(sum(S1d.^2))); % cdf of S
end
ind=find(S1acc>=.85);
%%
  
% 2. Projection of the network on the selected sinhuname_uniquelar vector 
%(calculate the distance[2-norm] of the new constructed sinhuname_uniquelar vector space)
BB1=[];BB2=[];
V = V1(:,1:ind(1));gamma = A*V;
U = U1(:,1:ind(1));alpha =(U'*A)';      
gammaNaN=isnan(gamma);modify=[];modify=find(gammaNaN==1);gamma(modify)=0;
BB1=sqrt(sum(gamma.^2,2));%BB1則為下游蛋白與基因的projection vector length
alphaNaN=isnan(alpha);modify=[];modify=find(alphaNaN==1);alpha(modify)=0;
BB2=sqrt(sum(alpha.^2,2));%BB2為上游基因的projection vector length
save('Result/ AlphaGamma.mat','BB1','BB2', 'alpha','gamma','-v7.3')
save('Result/BB1.mat','BB1', '-v7.3')
save('Result/BB2.mat','BB2', '-v7.3')
save('Result/Gup.mat','Gup', '-v7.3')
save('Result/Gdn.mat','Gdn', '-v7.3')

