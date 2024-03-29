function pcs = stackFastPCA(stackA,k);
%STACKFASTPCA
%PCS = STACKFASTPCA(STACK);
%PCS = STACKFASTPCA(STACK,K);
global stack;
if nargin < 2;k = 6;end;
    
if ndims(stack) > 2
    [ny,nx,nt]=size(stack);
    stack = reshape(stack,[ny*nx],[]);
end;

[nxy,nt]=size(stack);


pcs.spatial_av = mean(stack,2); 
pcs.temporal_av = (mean(stack,1)-mean(pcs.spatial_av))/mean(pcs.spatial_av);

fprintf('Removing spatial average\n');
% center, remove spatial average
spatial_av = mean(stack,2);

for iframe = 1:size(stack,2);
    stack(:,iframe) = stack(:,iframe)-spatial_av;
end

%stack = bsxfun(@minus,stack,spatial_av);
fprintf('Normalizing stack\n');
% normalize to ~ unit variance, assumes variance proportional to mean
pcs.spatial_norm = max(spatial_av,1);

for iframe = 1:size(stack,2);
    stack(:,iframe) = stack(:,iframe)./pcs.spatial_norm;
end
%stack = bsxfun(@rdivide,stack,pcs.spatial_norm);

fprintf('Removing temporal average\n');
% center, remove temporal average 
temporal_av = mean(stack,2);

for iframe = 1:size(stack,2);
    stack(:,iframe) = stack(:,iframe)-temporal_av;
end

%stack = bsxfun(@minus,stack,temporal_av);
fprintf('Computing principal components\n');
tic;[pcs.u,pcs.s,pcs.v]=pca(stack,k);toc

% de-normalize
pcs.u = bsxfun(@times,pcs.u,reshape(pcs.spatial_norm,nxy,1));

pcs.U = reshape(pcs.u,[ny,nx,k]);

return;
