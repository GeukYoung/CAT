function [RES50] = FxImage_CalRes50(Image,midValue)
% input
% Image : reconstruction image
% midValue : middle of recon Image value ex) imageslice case : 128
%
% output
% RES50 : 50% Resolution for fabric

if nargin < 2
    midValue = 0;
    disp('mid value set 0');
end

th_max50 = (max(max(Image))-midValue)*0.5+midValue;
Image(Image<th_max50) = 0;
Image(Image>=th_max50) = 1;
element_over_50 = sum(sum(Image));
total_element = size(Image,1)^2;
RES50 = (1-sqrt((element_over_50-1)/total_element))*100;
disp(['RES50 : ',num2str(RES50)]);