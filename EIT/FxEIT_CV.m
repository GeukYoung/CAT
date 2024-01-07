function [CV] = FxEIT_CV(Sigma,Mask,display_flag)
%% prepare grid image
% temp = size(Sigma);
% temp(temp==1) = [];
% Sigma_dim = length(temp);
% if Sigma_dim == 1
%     nPixel = 256; % default
%     Image = FxEIT_Tri2Grid(Element,Node,Sigma,nPixel);
%     if nargin > 3
%         display_flag = 1;
%     end
% elseif Sigma_dim == 2
%     Image = Sigma;
%     if nargin > 1
%         display_flag = 1;
%     end
% end
display_flag = 0;

Image = Sigma.*Mask;
%% find CV
temp = abs(reshape(Image,size(Image,1)*size(Image,2),1));
temp(isnan(temp)) = [];
temp(temp == 0) = [];
CV = std(temp)/mean(temp);
disp(['CV : ',num2str(CV)]);

%% Display
try
    if display_flag == 1
        imagesc(Image); axis image off; colormap hot;
        title(['CV : ',num2str(CV)]);
        set(gcf,'color','white');
    end
end

%%