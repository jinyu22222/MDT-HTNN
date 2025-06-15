function Z = tensor_allprod1(G,U,tr,Dg)
% Inputs:
% G : ones(size(T)), where T is  the input incomplete tensor;
% U : the cells of the duplication matrix;
% tr : tr = 0;?
% Dg : the size of each dimension;  is changing...

  N = length(Dg);
  Z = G;
  for n = 1:N
    if ~isempty(U{n})
      if tr == 0
        Z = tmult_MDT(Z,U{n},n,Dg);
        Dg(n) = size(U{n},1);
      else
        Z = tmult_MDT(Z,U{n}',n,Dg);
        Dg(n) = size(U{n},2);
      end
    end
  end
