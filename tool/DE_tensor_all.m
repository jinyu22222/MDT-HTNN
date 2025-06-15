function [ V S ] = DE_tensor_all( X, tau, S)
  % Delay Embedding for all directions (modes);
  % Inputs:
  % X : ones(size(T)), where T is  the input incomplete tensor;
  % tau : delay-embedding parameter:
  % S : the initialized duplication matrix, e.g. [];
  % Outputs:
  % V : the X*1 S(1)*2 S(2)......;
  % S : the cells of the duplication matrix;

  N  = ndims(X);
  II = size(X);
  if isempty(S)                     
    for n = 1:N
      S{n} = make_duplication_matrix(II(n),tau(n));
    end
  end  
  I2= II - tau + 1;
  JJ= [tau; I2];
  JJ = JJ(:);
  V = tensor_allprod1(X,S,0,size(X));
  V = reshape(full(V),[JJ']);

end
