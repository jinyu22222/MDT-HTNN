function [A,Tnew,Da]=tmult_MDT(T,M,n,Dt)
% The Tensor multiplication also named the n-mode multiplication
% turns for instance a 3-way array T(m x p x q) by multiplying M(r x p) 
% along second dimension into A(m x r x q)
%
% Input
% T :  tensor;
% M :  Jxsize(T,n)-matrix;
% n :  dimension to multiply;
% Dt : 
% Output:
% A :  tensor;
% Tnew : matrix given by M*matrizicing(T,n);-*
%Da : ?

Dm=size(M);
Da=Dt;
Da(n)=Dm(1);
Tn=unfold(T,n);
Tnew=M*Tn;
A=fold(Tnew,n,Da);

     

    
    
