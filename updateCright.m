function [Cright]=updateCright(Cright,B,X,A)
% this function execute a series of contraction to evaluate 
% sequentially the overlap depicted in the right of fig. 28 pag 133 
% of Schollwock, in the case in which an mpo is also present.

% it is the efficient way of evaluating expectation values for mpo's.

if isempty(X), X=reshape(eye(size(B,3)),[1,1,2,2]); end

Cright=contracttensors(A,3,2,Cright,3,3); %[aux1,aux2,phys1][aux3,aux4,aux5]->[aux1,phys1,aux3,aux4]
Cright=contracttensors(X,4,[2,4],Cright,4,[4,2]); %[aux5,aux6,phys2,phys3][aux1,phys1,aux3,aux4]->[aux5,phys2,aux1,aux3]
Cright=contracttensors(conj(B),3,[2,3],Cright,4,[4,2]);%[aux6,aux7,phys3][aux5,phys2,aux1,aux3]->[aux6,aux5,aux1]