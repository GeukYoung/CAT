function [Image] = FxEIT_Tri2Grid(Element,Node,Sigma,nPixel,margin_FOV)
% input
%   Element : FEM Face
%   Node    : FEM Node
%   Sigma   : FEM Face value
% output
%   Image   : Grid Image

if nargin < 4
    nPixel = 256;
%     margin_FOV = 1.05;
    margin_FOV = 1;
elseif nargin < 5
%     margin_FOV = 1.05;
    margin_FOV = 1;
end

meshsize = max(max(abs(Node)))*margin_FOV;

for i = 1:size(Element,1)
    xy(i,:) = mean(Node(Element(i,1:3),:));
end

ti = -meshsize:(2*meshsize)/(nPixel-1):meshsize;
[qx,qy] = meshgrid(ti,ti);

if length(Element) ~= size(Sigma,1)
    Sigma = Sigma';
end

for i = 1:size(Sigma,2)
%     F = TriScatteredInterp(xy(:,1),xy(:,2),Sigma(:,i));
    F = scatteredInterpolant(xy(:,1),xy(:,2),Sigma(:,i));
%     interpolateSolution
%     F.Method = 'spline';
    Image(:,:,i) = F(qx,qy);
    disp([num2str(i) ' / ' num2str(size(Sigma,2))]);
end
