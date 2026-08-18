function vec = triu2vec(X,k,ind)
%TRIU2VEC
% VEC = TRIU2VEC(X)
% VEC = TRIU2VEC(X,K)
% VEC = TRIU2VEC(X,K,IND) returns 
%
% 

if nargin < 2
    k= 1 ;
end

if nargin > 2
    X = X(ind,ind);
end

indices = find(triu(ones(size(X)),k));

vec = X(indices);

return;

