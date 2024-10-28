clc;
clear all;

[~,~,list]=xlsread('Expression.xlsx');
V=list;
V=cell2mat(V);
%linspace(start,end,linspace)
yy_EBV=spline([0 240 480 1440],V,linspace(0,1440,1440)); 

yyV_NAN=find(isnan(yy_EBV));%check
yyV_INF=find(isinf(yy_EBV));%check

V_max=max(V,[],2);
yy_EBV(find(yy_EBV<0))=0;

for j=1:size(V,1);
     yy_EBV(j,find(yy_EBV(j,:)>max(V(j,:))))=max(V(j,:));
 end
 yyV_max=max(yy_EBV,[],2);
 GRN_EXP=sparse(yy_EBV);
 save('GRN_EXP.mat','GRN_EXP')