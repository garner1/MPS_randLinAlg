function [mps,U]=prepare(mps,direction) 

% this function bring the mps into left-canonical form (direction='lr'; 
% 'rl' for right-canonical); the last 
% U is the norm of the state: sqrt(<mps|mps>); the output mps is normalized;
% Done according to pag 129 of Ulrich Schollwock, Annals of Physics 326, 96
N=length(mps);
switch direction    
    case 'lr'
        for i=1:N-1 
            [mps{i},U]=prepare_onesite(mps{i},'lr');
            mps{i+1}=contracttensors(U,2,2,mps{i+1},3,1);
        end
        [mps{N},U]=prepare_onesite(mps{N},'lr');
        
    case 'rl'
        for i=N:-1:2 
            [mps{i},U]=prepare_onesite(mps{i},'rl');
            mps{i-1}=contracttensors(mps{i-1},3,2,U,2,1);
            mps{i-1}=permute(mps{i-1},[1,3,2]); 
        end
        [mps{1},U]=prepare_onesite(mps{1},'rl');
end

