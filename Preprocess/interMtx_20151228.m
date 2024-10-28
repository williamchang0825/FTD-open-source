function [M]=interMtx_20151228(p1,p2,n)

switch(n)
    case 1
load PPIDIP_BIND_BIOGRIDp_IntAct_MINT20160616;MM=PPI;
    case 2
load GRNTFGnoGSEAsp;MM=GRN;
    case 3
load miTG_TargetScansp;MM=miTG;
end
M=sparse(length(p1),length(p2));
M1=sparse(length(p1),size(MM,2));
ind=[];

for i=1:length(p1);
    i
    ind{i,1}=[];
    if ~isempty(find(strcmpi(p1{i},Pu)));
        ind{i,1}=find(strcmpi(p1{i},Pu));
%     else
%         nn=[];
%         id1=[];
%         id1=find(strcmpi(p1{i},Gallnames));
%         if ~isempty(id1);
%             id1=mod(id1,size(Gallnames,1));
%             id1(find(id1==0))=size(Gallnames,1);
%             for j=1:length(id1);
%                 for k=1:size(Gallnames,2);
%                     if ~isempty(Gallnames{id1(j),k});
%                         if ~isempty(find(strcmpi(Gallnames{id1(j),k},Pu)));
%                             ind{i,1}=[ind{i,1} reshape(find(strcmpi(Gallnames{id1(j),k},Pu)),1,length(find(strcmpi(Gallnames{id1(j),k},Pu))))];
%                         end;
%                     end;
%                 end;
%             end;
%         else
%             ind{i,1}=[];
%         end;
    end;
    ind{i,1}=unique(ind{i,1});
    if ~isempty(ind{i});
        M1(i,:)=sum(MM(ind{i},:),1);
    end;
%     if n~=3;
%         for i=1:length(ind);
%             if ~isempty(ind{i});
%                 M(:,i)=sum(M1(:,ind{i}),2);
%             end;
%         end;
%     end;
end
if n~=3;
    for i=1:length(ind);
        if ~isempty(ind{i});
            M(:,i)=sum(M1(:,ind{i}),2);
        end;
    end;
end;
% if n==3;
%     indm=[];
%     for i=1:length(p2);
%         indm{i,1}=[];
%         if ~isempty(find(strcmpi(p2{i},Miu)));
%             indm{i,1}=find(strcmpi(p2{i},Miu));
%         else
%             id1=[];
%             id1=find(strcmpi(p2{i},Gallnames));
%             if ~isempty(id1);id1=mod(id1,size(Gallnames,1));
%                 id1(find(id1==0))=size(Gallnames,1);
%                 for j=1:length(id1);
%                     for k=1:size(Gallnames,2);
%                         if ~isempty(Gallnames{id1(j),k});
%                             if ~isempty(find(strcmpi(Gallnames{id1(j),k},Miu)));
%                                 indm{i,1}=[indm{i,1} reshape(find(strcmpi(Gallnames{id1(j),k},Miu)),1,length(find(strcmpi(Gallnames{id1(j),k},Miu))))];
%                             else
%                                 indm{i,1}=[];
%                             end;
%                         end;
%                     end;
%                 end;
%             end;
%         end;
%         indm{i,1}=unique(indm{i,1});
%         if ~isempty(indm{i});
%             M(:,i)=sum(M1(:,indm{i}),2);
%         end;
%     end;
% end
M=(M~=0);
if n==1;M=(M | M');end