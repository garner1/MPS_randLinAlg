function [Cstorage]=initCstorage(mpsB,mpoX,mpsA,N)
%this function seems to work with initial R-canonical state, 
%and sweeps procedure going first from L to R
%because of the updateCright function, and not updateCleft


%Cright is a tensor with 3-indeces because we are considering also mpo; 
% if it was only overlap of states it would have 2 indeces;
Cstorage=cell(1,N+1); 
Cstorage{1}=1; 
Cstorage{N+1}=1;
for i=N:-1:2
    if isempty(mpoX), X=[]; else X=mpoX{i}; end
    Cstorage{i}=updateCright(Cstorage{i+1},mpsB{i},X,mpsA{i}); 
end