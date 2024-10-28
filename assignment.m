clc;clear; 
load('GRN_name.mat');% load in the host gene names in your disease
%GRN_name = readmatrix('genesymbol.xlsx') 
load('TF_list_with_alias.mat');
load('receptor_DB.mat');load('lncRNAdatabase.mat');load('lncRNA_database.mat')
InRcp = [];InMir=[];InLnc = [];InTF = [];
%% If there are pathogen nodes in your nodes， please enter the viral protein entities below:
InVirus= {}; % Enter your node names 
fprintf('Start Arrangement');
for i = 1:length (GRN_name)
        if strncmpi(GRN_name(i),'LINC',4)||strncmpi(GRN_name(i),'LOC',3)|| any(strcmpi(GRN_name(i),lncRNA),'all')|| any(strcmpi(GRN_name(i),lncRNA_name),'all')
           InLnc = cat(1,InLnc,GRN_name(i));
        elseif strncmpi(GRN_name(i),'MIR',3)
           InMir= cat(1,InMir,GRN_name(i));
        elseif any(strcmpi(GRN_name(i),receptor_DB))
           InRcp = cat(1,InRcp,GRN_name(i));
        elseif any(strcmpi(GRN_name(i),TF_list_with_alias))
           InTF = cat(1,InTF,GRN_name(i));
        end
end
InGene = setdiff(GRN_name,union(union(union(InMir,InRcp),InLnc),InTF)); 
Protein = cat(1,InGene,InRcp,InTF);
Name = cat(1,Protein,InMir,InLnc);
%% Arrangement

fprintf('PPI Construction');
[PPI]=interMtx_20151228(Protein,Protein,1);
fprintf('GRN Construction');
[GRN_HH]=interMtx2_20160324(Name,Name,2);
%% Add Virus
if isempty(InVirus)
    PPI_name = Protein;
    GRN_name = Name;
    GRN = GRN_HH;
else
    PPI_name = cat(1,Protein,InVirus);
    GRN_name = cat(1,Name,InVirus);
    GRN = ones(length(GRN_name));
    GRN(1:length(GRN_HH),1:length(GRN_HH)) = GRN_HH;
end
num1 = sum(PPI,'all');num2 = sum(GRN,'all');
save('GRN_name.mat','GRN_name');save('PPI_name.mat','PPI_name');
save('PPI.mat','PPI');save('GRN.mat','GRN');
