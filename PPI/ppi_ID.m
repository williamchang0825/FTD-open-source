clear;clc;
%% loading data
load('PPI_EXP.mat');load('PPI.mat');load('PPI_name.mat');
interaction=PPI; gene_expre=PPI_EXP; Name = PPI_name;
nnn=size(PPI_EXP,2); % n times
Final_reg_ability = double(zeros(size(interaction),'like',interaction));%創造一個零矩陣跟PPI一樣大但是格式是double 用來存結果
basal = zeros(size(interaction,1),1);%創造一個跟PPI依樣大的column
[BB,index]=sort(sum(interaction(:,1:end)~=0,2));%把PPI每個row加起來再做排序 升冪排序成一個col:BB index存BB對應到原本排序前的位置
id=find(BB~=0,1,'first');%找到BB中第一個非0的數是第幾個放入id
%% ID
tic % 一個函數 主要用來計時 代表計時開始 下一個打toc代表計時結束 就能算總共花的時間
fprintf('start ID(from %d to %d)\n',size(gene_expre,1)-id,0)
start = ; % start order
for j = start:size(gene_expre,1)
    timer = tic;
    i = index(j); % index of the instant protein %j代表處理第幾個數 所以這邊i指的是從小到大第j個數在PPI的位置
    mmm = size(gene_expre,1)-j; % the amount of protein left without processing
    All = 1:size(interaction,2); % full size matrix
    bind1=find(interaction(i,:)==1); % ppi information for 1 of the instant protein%找到等於1的數的座標位置*是一個raw
    bind0=setdiff(All,bind1);%代表數字是0的位置
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
        for m = 1:nnn/2 % Foward
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
    basal(i,:) = thetaout(end);
    null = find(thetaout(1:end-3)==0);
    bindout(null)= [];thetaout(null)= [];
    for n = 1:length(bindout)
        Final_reg_ability(i,bindout(n))=thetaout(n);
    end
    A = full(Final_reg_ability);
    if mod(mmm,10)==0
        fprintf('%d (%d nodes --> %d nodes)  ',mmm,length(bind1),length(bindout))
        fprintf('Overall time = %f sec\n', toc(timer));
    end
    if mod(mmm,100)==0
        save PPI_ID_TEMP;
    end     
end
basal(find(basal>0))=1;basal(find(basal<0))=-1;
fprintf('Done\n')
fprintf('Interaction:[%6d ------> %-6d]\n',size(find(interaction~=0),1),size(find(A~=0),1))
fprintf('       Node:[%6d ------> %-6d]\n',size(interaction,1),length(find(sum(A,2)~=0)))
toc % elapsed time
 %%
PPI_PNP = A;
PPI_Edge = (A+A')./2;
Basal_PPI = table(Name,basal);
save('Result/PPI_PNP.mat','PPI_PNP', '-v7.3')
save('Result/PPI_Edge.mat','PPI_Edge', '-v7.3')
writetable(Basal_PPI,'Result/Basal_PPI.txt')