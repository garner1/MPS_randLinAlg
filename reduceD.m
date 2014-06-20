function [mpsB,norma,dist]=reduceD(mpsA,mpoX,DB,precision)

% PROVIDES THE VATIATIONALLY REDUCED MPSB, AFTER APPLYING MPOX TO MPSA;
% norma IS THE NORM OF THE OUTPUT STATE; dist IS THE SQUARED 2-NORM
% DISTANCE FROM THE EXACT STATE

N=length(mpsA); 
mpsB=truncate(mpomps(mpoX,mpsA),'rl',DB); 
% mpsB=createrandommps(N,DB,2);

% initialization of the storage 
Cstorage=initCstorage(mpsB,mpoX,mpsA,N);

% optimization sweeps 
while 1
    Kvalues=[];
    
    % ****************** cycle 1: j -> j+1 (from 1 to N-1) ****************
    for j=1:(N-1)
        Cleft=Cstorage{j};
        Cright=Cstorage{j+1};
        A=mpsA{j}; X=mpoX{j}; 
        [B,K]=reduceD2_onesite(A,X,Cleft,Cright);
        [B,U]=prepare_onesite(B,'lr');
        mpsB{j}=B; 
        Kvalues=[Kvalues,K];
        % storage-update
        Cstorage{j+1}=updateCleft(Cleft,B,X,A); 
    end
    
    % ****************** cycle 2: j -> j-1 (from N to 2) ******************
    for j=N:(-1):2
        Cleft=Cstorage{j};
        Cright=Cstorage{j+1};
        A=mpsA{j}; X=mpoX{j}; 
        [B,K]=reduceD2_onesite(A,X,Cleft,Cright);
        [B,U]=prepare_onesite(B,'rl');
        mpsB{j}=B; Kvalues=[Kvalues,K];
        % storage-update
        Cstorage{j}=updateCright(Cright,B,X,A); 
    end
    
    err=std(Kvalues)/abs(mean(Kvalues));
    
    if err<precision, 
        mpsB{1}=contracttensors(mpsB{1},3,2,U,2,1); 
        mpsB{1}=permute(mpsB{1},[1,3,2]);
        norma=sqrt(-K);
        dist=1+K;
        break;
    end
end


    
    