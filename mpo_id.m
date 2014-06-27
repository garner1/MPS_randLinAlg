function [mpo]=mpo_id(N,d)

% CREATE THE IDENTITY MPO OPERATOR, ON N QUDITS WITH PHYSICAL DIMENSION d

mpo=cell(1,N);

for ind=1:N
    mpo{ind}=zeros(1,1,d,d);
    mpo{ind}(1,1,:,:)=eye(2); 
end