function [output] = FxFilt_Norm(input)
[a, b] = size(input);
temp = input - repmat(min(min(input)),a,b);
output = temp./max(max(temp));