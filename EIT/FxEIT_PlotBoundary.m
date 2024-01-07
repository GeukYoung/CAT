function [ROI_contour_outside] = FxEIT_PlotBoundary(bnd,color,width,op)
bnd(~isnan(bnd)) = 1;
bnd(isnan(bnd)) = 0;
if bnd(1,1) == 1
    bnd = ~bnd;
end

if nargin>3
se = strel('diamond',op);
bnd = imdilate(bnd,se);
end

[temp(1) temp(2)] = find(bnd,1);
ROI_contour_outside = bwtraceboundary(bnd,temp,'W',8,Inf,'counterclockwise');

if nargin < 2
    color = 'w';
    width = 2;
elseif nargin < 3
    width = 2;
end
hold on; plot(ROI_contour_outside([1:1:end 1],2),ROI_contour_outside([1:1:end 1],1),color,'LineWidth',width);
end

