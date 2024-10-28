function [M]=candidate_network_construction(  p1 , p2 , p3 )

%% 參數名稱 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p1           : 
% p2           : 
% p3           : "PPI"     "GRN" 
% all_node : database 所擁有的gene 名  
% M1          :
% M            :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("======= Candidate_%s_construction started =======\n", p3)

%%  (human human) (human virus) (virus virus) 
if p3=="PPI"    
    MM         =  matfile('PPIDIP_BIND_BIOGRIDp_IntAct_MINT20160616.mat').PPI ; 
    all_node =   matfile('PPIDIP_BIND_BIOGRIDp_IntAct_MINT20160616.mat').Pu ;       
end
%%
if p3=="GRN"     
     MM           = matfile("GRNHTRIdb_ITFP_MRNTargetScanspstarBase2_bothCircuitDBcomplex20160324.mat").GRN;
     all_node   =  matfile("GRNHTRIdb_ITFP_MRNTargetScanspstarBase2_bothCircuitDBcomplex20160324.mat").Pu;
end
%% 


tic
index =  find( ismember( all_node , p1)  ) ;       
a= cell2table(  cat(2,  all_node(index) ,  num2cell(index) ) , 'VariableNames',  {'name' ,'index_in_database' } ) ;    
c= outerjoin( a , cell2table(p1 ,'VariableNames',{'name' })  ,'MergeKeys',true);   %   load( convertCharsToStrings( pwd  )  +"\database\GRN_cov2_human_database\"+"GRN_merge.mat"  ); MM=GRN_merge ; clear('GRN_merge'); 

c= table2cell(c) ;  c(:,1)=[] ;                                        % query gene 在 database index (nan: database 沒有)

M1= sparse( length(p1) , size(MM,2)  );
M1( find(~cellfun(@isnan, c)) , : )=MM(index, : ) ;  %  find(~cellfun(@isnan, c)) = 不為nan 在c 中的inde

M  = sparse( length(p1) , length(p2)   );
M( : , find(~cellfun(@isnan, c))    )=M1( :  , index  ) ;  
M=(M~=0);  % sparse to logical array

if p3=="PPI"
    M=(M | M');
    toc
else
    toc
end

fprintf("===== Candidate_%s_construction completed =======\n", p3)
end
