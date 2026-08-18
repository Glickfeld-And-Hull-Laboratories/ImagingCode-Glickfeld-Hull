function [opt_Az, opt_El] = findThreshStimLocation(fit_true_vec, goodfit_ind)

Az_centers = fit_true_vec(goodfit_ind, 5);
El_centers = fit_true_vec(goodfit_ind, 4);

threshold = 7.5; % deg
res = 0.5;       % grid resolution (deg)

az_range = min(Az_centers)-threshold : res : max(Az_centers)+threshold;
el_range = min(El_centers)-threshold : res : max(El_centers)+threshold;

[AZ, EL] = meshgrid(az_range, el_range);
n_covered = zeros(size(AZ));

for i = 1:numel(AZ)
    d = sqrt((Az_centers - AZ(i)).^2 + (El_centers - EL(i)).^2);
    n_covered(i) = sum(d < threshold);
end

[max_covered, idx] = max(n_covered(:));
opt_Az = AZ(idx);
opt_El = EL(idx);

fprintf('Optimal stimulus location: Az=%.2f deg, El=%.2f deg\n', opt_Az, opt_El)
fprintf('Cells within %.1f deg: %d / %d (%.0f%%)\n', ...
    threshold, max_covered, numel(Az_centers), 100*max_covered/numel(Az_centers))
end