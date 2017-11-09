function [I_equal, LUT] = histogram_equal(I, T)

hist_in_org = imhist(I);

hist_in = hist_in_org((T+1):end); % take only relevant values according to threshold

norm_hist_in = (hist_in.')./sum(hist_in);
LUT = cumsum(norm_hist_in); % create suitable LUT

if (T>0)
    LUT_unchanged = [0:1:(T-1)] / 255; % LUT of values under threshold
    LUT = [LUT_unchanged'; (LUT_unchanged(end)+LUT*(255-T)/255)']'; % LUT of all values
end

I_equal = uint8(round(255*LUT(double(I)+1))); % process the image according to LUT

end