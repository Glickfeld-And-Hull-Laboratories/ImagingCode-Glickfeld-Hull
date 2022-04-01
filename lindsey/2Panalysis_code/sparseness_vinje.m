function S = sparseness_vinje(resp_mat);
    %calculates sparseness of responses to multiple stimuli as in Vinje and
    %Gallant 2000. Values near 0 are dense, and near 1 are sparse.
    %resp_mat is ncells x nstim
    sz = size(resp_mat);
    x = (sum(resp_mat,2)./sz(2)).^2;
    y = sum(((resp_mat.^2)./sz(2)),2);
    z = 1-(1/sz(2));
    S = squeeze((1-(x./y))./z);
end