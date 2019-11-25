function [res,g] = dualfune(beta,L,Y,E)

num = length(Y); 
K = size(beta,1);

CombinedL = E;
for k=1:K  
  CombinedL = CombinedL + beta(k)*L{k};
end

Z = CombinedL\Y;
res = Y'*Z;

g = zeros(K,1);
for k=1:K
  g(k) = - Z'*L{k}*Z;
end
