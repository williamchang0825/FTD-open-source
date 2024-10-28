clear;clc;
tic % start time
%% loading data
load('GRN_EXP.mat');load('GRN.mat');load('GRN_name.mat');
interaction=GRN; gene_expre=GRN_EXP; Name = GRN_name;
nnn=size(GRN_EXP,2); % n times
%% ID
Final_reg_ability = double(zeros(size(interaction),'like',interaction));
basal = zeros(size(interaction,1),1);
[BB,index]=sort(sum(interaction(:,1:end)~=0,2));
id=find(BB~=0,1,'first');
fprintf('start ID(from %d to %d)\n',size(gene_expre,1)-id,0)
for j = id:size(gene_expre,1)
    i = index(j); % index of the instant protein
    mmm = size(gene_expre,1)-j; % the amount of protein left without processing
    All = 1:size(interaction,2); % full size matrix
    bind1=find(interaction(i,:)==1); % ppi information for 1 of the instant protein
    bind0=setdiff(All,bind1);
    bind = bind0;
    bindout = bind1;
    % initial AIC_value
    AIC_value_initial = 10000;
    [X,phi,theta,resnorm]=linear_regression(i,nnn,gene_expre,bind);
    AIC_value=AICValue(X,theta,resnorm);
    thetaout=theta;
    if length(bind1) > 1
        Bind1={};Theta={};AIC_value_pack={};
        binding=bind1;
        for m = 1:nnn % Foward
            for u = 1:length(binding)
                if isempty(Bind1)
                    bindtemp = binding(u);
                else
                    bindtemp = cat(2,Bind1{m-1},binding(u));
                end
                bind = setdiff(All,bindtemp);
                % linear regression
                [X,phi,theta_temp,resnorm_temp]=linear_regression(i,nnn,gene_expre,bind);
                pre_AIC_value=AICValue(X,theta_temp,resnorm_temp);
                if u==1
                    if ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                        AIC_value_temp = pre_AIC_value;
                        bindinfo = bindtemp;
                    else
                        AIC_value_temp = AIC_value_initial;
                    end
                else
                    if pre_AIC_value<=AIC_value_temp && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                        AIC_value_temp = pre_AIC_value;
                        bindinfo = bindtemp;
                    end
                end
            end
            Bind1{m}= bindinfo;
            Theta{m}= theta_temp;
            AIC_value_pack{m}=AIC_value_temp;
            binding = setdiff(bind1,Bind1{m});
        end
        [AIC_value,min_index] = min(cell2mat(AIC_value_pack));
        bindout= Bind1{min_index};
        if length(bindout) > 1
            clear bindinfo butheta_temp AIC_value_temp
            for y=1:length(bindout)-1 % Backward
                for z = 1:length(bindout)
                    bindtemp = bindout;
                    bindtemp(z)= [];
                    bind = setdiff(All,bindtemp);
                    % linear regression
                    [X,phi,theta_temp,resnorm_temp]=linear_regression(i,nnn,gene_expre,bind);
                    pre_AIC_value=AICValue(X,theta_temp,resnorm_temp);
                    if z == 1
                        if ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                        AIC_value_temp = pre_AIC_value;
                        bindinfo = bindtemp;
                        else 
                            AIC_value_temp = AIC_value_initial;
                        end
                    else
                        if pre_AIC_value<=AIC_value_temp && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                            AIC_value_temp = pre_AIC_value;
                            bindinfo = bindtemp;
                        end
                    end
                end
                bindout= bindinfo;
                Bind1{m+y}= bindinfo;
                Theta{m+y}= theta_temp;
                AIC_value_pack{m+y}=AIC_value_temp;
            end
        end
        [AIC_value,min_index] = min(cell2mat(AIC_value_pack));
        thetaout=Theta{min_index};
        bindout= Bind1{min_index};
    end
%     Information=cat(1,cat(2,Bind1,bindout),cat(2,AIC_value_pack,AIC_value));
    Difference{i} = setdiff(bindout,bind1);
    Difference2{i}= setdiff(bind1,bindout);
    basal(i) = thetaout(end);
    null = find(thetaout(1:end-3)==0);
    bindout(null)= [];thetaout(null)= [];
    for n = 1:length(bindout)
        Final_reg_ability(i,bindout(n))=thetaout(n);
    end
    A = full(Final_reg_ability);
    if mod(mmm,10)==0
        fprintf('%d\n',mmm)
    end
end
fprintf('Done\n')
fprintf('Interaction:[%6d ------> %-6d]\n',size(find(interaction~=0),1),size(find(A~=0),1))
fprintf('       Node:[%6d ------> %-6d]\n',size(interaction,1),length(find(sum(A,2)~=0)))
toc % elapsed time
%%
GRN_PNP = A;
GRN_Edge = A;
Basal_GRN = table(Name,basal);
save('Result/GRN_PNP.mat','GRN_PNP')
save('Result/GRN_Edge.mat','GRN_Edge')
writetable(Basal_GRN,'Result/Basal_GRN.txt')