% this function provides the energy E, and ground state mps (with bond-dim D), of an input Hamiltonian (mpo)
% mpsB (optional if set to []) is a random initial state used to start the algorithm
function [E,mps]=minimizeE(mpo,D,precision,mpsB)

N = size(mpo); 
d = 2; 
mps = createrandommps(N,D,d); 
mps = prepare(mps);

% storage-initialization
Hstorage = initHstorage(mps,mpo,d);
if ~isempty(mpsB), Cstorage = initCstorage(mps,[],mpsB,N); end

% optimization sweeps 
while 1
    Evalues = [];

    % ****************** cycle 1: j -> j+1 (from 1 to N-1) **************** 
    for j = 1:(N-1)
        % projector-calculation 
        if ~isempty(mpsB)
            B = mpsB{j};
            Cleft = Cstorage{j};
            Cright = Cstorage{j+1}; 
        end

        % optimization
        Hleft = Hstorage{j};
        Hright = Hstorage{j+1};
        hsetj = mpo{j};
        [A,E] = minimizeE_onesite(hsetj,Hleft,Hright); 
        [A,U] = prepare_onesite(A,'lr');
        mps{j} = A; 
        Evalues = [Evalues,E];
 
        % storage-update 
        Hstorage{j+1} = updateCleft(Hleft{m},A,hsetj,A); 
        if ~isempty(mpsB) 
            Cstorage{j+1} = updateCleft(Cleft,A,[],B);
        end
    end
   
    % ****************** cycle 2: j -> j-1 (from N to 2) ****************** 
    for j = N:(-1):2
        % projector-calculation 
        if ~isempty(mpsB)
            B = mpsB{j};
            Cleft = Cstorage{j};
            Cright = Cstorage{j+1}; 
        end
        
        % minimization
        Hleft = Hstorage{j};
        Hright = Hstorage{j+1};
        hsetj = mpo{j}; 
        [A,E] = minimizeE_onesite(hsetj,Hleft,Hright); 
        [A,U] = prepare_onesite(A,'rl');
        mps{j} = A; 
        Evalues = [Evalues,E];

        % storage-update 
        Hstorage{j}=updateCright(Hright{m},A,hsetj,A); 
        if ~isempty(mpsB) 
            Cstorage{j}=updateCright(Cright,A,[],B);
        end
    end
    if (std(Evalues)/abs(mean(Evalues))<precision) 
        mps{1}=contracttensors(mps{1},3,2,U,2,1); 
        mps{1}=permute(mps{1},[1,3,2]);
        break;
    end
end





