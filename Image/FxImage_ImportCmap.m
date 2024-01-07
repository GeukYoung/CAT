function [cmap] = FxImage_ImportCmap(option)
if nargin == 0
    option = 1;
end

switch option
    case 1
        cmap = load('Cmap1.m');
    case 2
%         load();
end