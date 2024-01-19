function [Consistency, core] = btd_corcond(T,Factors,L)

%This function is a prototype function for compute a core consistency diagnosis for BTD L,L,1
% 
% INPUT
% T: Tensor data array
% Factors: factor in CPD format as cell array
% L: rank Lr used in format of [L1, L2, .....Lr]
% OUTPUT
% The core consistency in from of percentage, with maximum number of 100%
%and the core tensor of the BTD(L,L,1)
% 
% This code was made by Noramon Dron under her PhD studentship and part of
% the work in preprint: https://arxiv.org/abs/2312.11151

    %% CPD format
    invers_A = pinv(Factors{1,1});
    invers_B = pinv(Factors{1,2});
    invers_c = pinv(Factors{1,3});

    core = nmodeproduct(T,invers_A,1);
    core = nmodeproduct(core,invers_B,2);
    core = nmodeproduct(core,invers_c,3);

    %create ideal identity matrix of size G and vectorize
    iden_size = size(core);
    %Identity = eye(iden_size);
    ident= zeros(iden_size);
    R=length(L);
    j=1;
    k=0;
    if R==1
        ident=eye(L);
    else
        for i=1:R
           k=k+L(i);
           ident(j:k,j:k,i) = eye(L(i));
           j=j+L(i);
        end
    end

    Identity = ident(:);%vectorize the ideal core tensor
    G=core(:); %vectorize the core tensor G

    %calculate the core consistency
    Gss=sum(G.^2); %sum of square of identity matrix
    Consistency=100*(1-sum((Identity-G).^2)/Gss);
    Consistency=real(Consistency);
%   end

end