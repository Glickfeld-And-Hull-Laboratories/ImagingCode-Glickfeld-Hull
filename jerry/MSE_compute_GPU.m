function mse_values = MSE_compute_GPU(signal, m, r_frac, max_scale)
% MSE_compute_GPU: Multiscale Entropy computation using GPU acceleration
%
% Inputs:
%   signal    - 1D time series (vector)
%   m         - embedding dimension (e.g., 2)
%   r_frac    - similarity threshold (fraction of std, e.g., 0.15)
%   max_scale - maximum scale to compute MSE
%
% Output:
%   mse_values - vector of MSE values across scales

    if isrow(signal)
        signal = signal';
    end

    signal_gpu = gpuArray(signal);
    mse_values = zeros(1, max_scale, 'gpuArray');

    for scale = 1:max_scale
        cg_signal = coarse_grain_GPU(signal_gpu, scale);
        r = r_frac * std(cg_signal);  % Calculate r for each coarse-grained signal
        mse_values(scale) = sample_entropy_GPU_Block(cg_signal, m, r);
    end

    mse_values = gather(mse_values); % Bring back to CPU
end

function cg = coarse_grain_GPU(signal_gpu, tau)
% Coarse-grain the signal: average non-overlapping windows of size tau
    N = floor(length(signal_gpu) / tau);
    cg = mean(reshape(signal_gpu(1:N*tau), tau, N), 1)';
end

function SampEn = sample_entropy_GPU_Block(signal_gpu, m, r)
% GPU-based Sample Entropy calculation using blockwise comparisons
    N = length(signal_gpu);

    % Create embedded matrices
    X_m  = zeros(N - m + 1, m, 'gpuArray');
    X_m1 = zeros(N - m, m + 1, 'gpuArray');

    for i = 1:m
        X_m(:,i) = signal_gpu(i:N - m + i);
    end
    for i = 1:m+1
        X_m1(:,i) = signal_gpu(i:N - m + i - 1);
    end

    block_size = 512;  % Tune this based on your GPU capability
    count_m = 0;
    count_m1 = 0;

    % --- Count for m ---
    for i = 1:block_size:size(X_m,1)
        i_end = min(i + block_size - 1, size(X_m,1));
        Xi = X_m(i:i_end, :);

        dists = max(abs(Xi - permute(X_m, [3 2 1])), [], 2);
        dists = squeeze(dists);

        % Ignore self-match by setting diagonal distances to Inf
        if i <= size(X_m,1)
            mask = (i:i_end)'; 
            dists(sub2ind(size(dists), (1:length(mask))', mask - i + 1)) = inf;
        end

        count_m = count_m + sum(dists(:) < r);
    end

    % --- Count for m+1 ---
    for i = 1:block_size:size(X_m1,1)
        i_end = min(i + block_size - 1, size(X_m1,1));
        Xi = X_m1(i:i_end, :);

        dists = max(abs(Xi - permute(X_m1, [3 2 1])), [], 2);
        dists = squeeze(dists);

        if i <= size(X_m1,1)
            mask = (i:i_end)';
            dists(sub2ind(size(dists), (1:length(mask))', mask - i + 1)) = inf;
        end

        count_m1 = count_m1 + sum(dists(:) < r);
    end

    if count_m == 0 || count_m1 == 0
        SampEn = NaN;
    else
        SampEn = -log(double(count_m1) / double(count_m));
    end
end
