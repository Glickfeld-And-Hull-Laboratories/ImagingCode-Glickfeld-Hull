function shifted = circshift_interp(data, shift_bins)
% circshift_interp  Fractional-bin circular shift via linear interpolation
%
%   shifted = circshift_interp(data, shift_bins)
%
% Input:
%   data       - 1 x N or M x N array (shifts along dim 2 if matrix)
%   shift_bins - scalar, number of bins to shift (positive = shift right)
%                Can be fractional (e.g., 2.7 bins)
%
% Output:
%   shifted    - same size as data, circularly shifted with interpolation
%
% Integer shifts match circshift exactly. Fractional shifts use linear
% interpolation between the two nearest integer-shifted versions.

    int_part  = floor(shift_bins);
    frac_part = shift_bins - int_part;

    if frac_part == 0
        % Pure integer shift — use built-in circshift
        shifted = circshift(data, int_part, 2);
    else
        % Interpolate between floor and ceil shifts
        shifted_lo = circshift(data, int_part, 2);
        shifted_hi = circshift(data, int_part + 1, 2);
        shifted = (1 - frac_part) * shifted_lo + frac_part * shifted_hi;
    end
end
