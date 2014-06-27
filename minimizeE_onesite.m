% ************************ one-site optimization **************************

function [A,E,Heff] = minimizeE_onesite(hsetj,Hleft,Hright)

DAl = size(Hleft,1); 
DAr = size(Hright,1); 
d = 2;

% calculation of Heff

Heff = 0; 
Heff = contracttensors(Hleft,3,2,Hright,3,2); 
Heff = contracttensors(Heff,5,5,hsetj,3,3); 
Heff = permute(Heff,[1,3,5,2,4,6]); 
Heff = reshape(Heff,[DAl*DAr*d,DAl*DAr*d]); 

%
% optimization
%
options.disp = 0; 
[A,E] = eigs(Heff,1,'sr',options); 

A=reshape(A,[DAl,DAr,d]);
