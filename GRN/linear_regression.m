function [X,phi,theta,resnorm]=linear_regression(i,nnn,gene_expre,bind)
X=gene_expre(i,2:nnn)';
pi=gene_expre(i,1:nnn-1)'; % 3x1
phi=gene_expre(:,1:nnn-1)'; % 3x100
phi(:,bind) = []; % exclude protein without binding   
phi=cat(2,phi.*pi,pi,pi,ones(nnn-1,1)); % [Gi(t) Pi(t) 1]
phi = full(phi);
b=[0;1]; 
A= [zeros(1,size(phi,2)-3),-1,0,0;zeros(1,size(phi,2)-3),0,1,0];
options=optimset('LargeScale','off','Display','off');warning('off');
[theta,resnorm]=lsqlin(phi,X,A,b,[],[],[],[],[],options); % resnorm = (X-phi*theta)'*(X-phi*theta)
