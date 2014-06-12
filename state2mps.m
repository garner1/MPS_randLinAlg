function [mps,V]=state2mps(coeff,d,N)

% GIVEN AN INPUT STATE, DEFINED IN TERMS OF THE COEFFICIENTS WITH
% RESPECT TO A ORTHORNORMAL GLOBAL BASIS IN THE HILBERT SPACE, PROCUCES THE
% LEFT CANONICAL (OR RIGHT CANONICAL IN THE SECOND COMMENTED PART) MPS STATE ASSOCIATED TO IT; THE INPUT ARE THE VECTOR OF
% COEFFICIENTS AND THE LOCAL BASIS DIMENSION, AND THE NUMBER OF QUDITS; THE
% OUTPUT IS THE MPS REPRESENTATION WITH THE NORM OF THE INPUT STATE

% % d=2;N=10;
% % coeff=rand(d^N,1);
% % coeff=coeff/norm(coeff);

mps=cell(N,1);
V=coeff;r=1;
for ind=1:N
    V=reshape(V,[r*d,d^(N-ind)]);
    [U,S,V]=svd2(V);
    rold=r;r=size(S,1);
    V=S*V;
    U=reshape(U,[rold,d,r]);
    mps{ind}=permute(U,[1,3,2]);
end
% 
% U=coeff;r=1;
% for ind=1:N
%     U=reshape(U,[d^(N-ind),r*d]);
%     [U,S,V]=svd2(U);
%     rold=r;r=size(S,1);
%     U=U*S;
%     V=reshape(V,[r,d,rold]);
%     mps{ind}=permute(V,[1,3,2]);
% end
