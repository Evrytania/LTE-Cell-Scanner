function r = filter_wo_tail(s, coef, convert_ratio)

if size(s, 1) == 1
    disp('filter_wo_tail: input signal must be column vector');
    return;
end

len_coef = length(coef);

if mod(len_coef, 2) == 0
    disp('filter_wo_tail: length of coef must be odd!');
    return;
end

if convert_ratio<=1 
    downsampling_ratio = round(1/convert_ratio);
    r = filter(coef, 1, [s; zeros(len_coef-1,1)]);
    r = r( (((len_coef-1)/2)+1) : (end-((len_coef-1)/2)) );
    r = r(1:downsampling_ratio:end);
else
    r = upsample(s, convert_ratio);
    r = filter(coef, 1, [r; zeros(len_coef-1,1)]);
    r = r( (((len_coef-1)/2)+1) : (end-((len_coef-1)/2)) );
end
