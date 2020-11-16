function hwhm = hwhmFromOriFit(f,orientations)
n = size(f,2);
[peak, peakInd] = max(f);
baseline = min(f);
amplitude = peak-baseline;
hwhm = nan(1,n);
for i = 1:n
    if all(f(:,i) == 0) || all(diff(f(:,i)) == 0)
        hwhm(i) = nan;
    else
        p = peakInd(i);
        f_temp = circshift(f(:,i),90-p);
        widthInd = f_temp >= (peak(i)-(0.5.*amplitude(i)));
        oriWidth = orientations(widthInd);
        hwhm(i) = oriWidth(end) - oriWidth(1);
    end
    if any(diff(oriWidth > 1)) & ~all(f(:,i) == 0)
        error('peak not centered enough')
    end
end

end