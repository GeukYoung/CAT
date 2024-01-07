function [temp2] = FxEIT_YBar(Image,op1,op2)
% op2 = 1 : parellel sum
% op2 = 2 : mean

%     Image = Data1.TV_grid(:,:,1);
    if nargin == 1 % full pixel number
        temp = nansum(Image,2);  
        subplot(121);
        imagesc(Image);
        axis tight image;
        set(gca,'xtick',[], 'ytick',[]);
        set(gca,'Ydir','normal');
        subplot(122);
        barh(temp,1);
    else % 
        if nargin > 2
            if op2 == 1
            %% parellel sum
            half = round(size(Image,1)/2);
            Image = 1./Image;
            Image(isinf(Image)) = NaN;

            temp(:,1) = nansum(Image(:,1:half),2);    
            temp(:,2) = nansum(Image(:,half+1:end),2);
            
            temp2(:,1) = nansum(reshape(temp(:,1),op1,size(temp,1)/op1));
            temp2(:,2) = nansum(reshape(temp(:,2),op1,size(temp,1)/op1));
            
%             temp2(:,2) = temp(1:op1:end,2) + temp(2:op1:end,2);

            temp2 = 1./temp2;
            temp2(isinf(temp2)) = 0;
            %% mean
            elseif op2 == 2 
                half = round(size(Image,1)/2);
                Image(isinf(Image)) = NaN;
                
                temp(:,1) = nanmean(Image(:,1:half),2);
                temp(:,2) = nanmean(Image(:,half+1:end),2);
                
                temp2(:,1) = nanmean(reshape(temp(:,1),op1,size(temp,1)/op1));
                temp2(:,2) = nanmean(reshape(temp(:,2),op1,size(temp,1)/op1));
                
                temp2(isinf(temp2)) = 0;
            end
        else
            %% normal sum
            half = round(size(Image,1)/2);
            Image(isinf(Image)) = NaN;

            temp(:,1) = nansum(Image(:,1:half),2);    
            temp(:,2) = nansum(Image(:,half+1:end),2);
            
            temp2(:,1) = nansum(reshape(temp(:,1),op1,size(temp,1)/op1));
            temp2(:,2) = nansum(reshape(temp(:,2),op1,size(temp,1)/op1));
            
            temp2(isinf(temp2)) = 0;
        end
            
%         figure;
%         subplot(121);
%         imagesc(1./Image);
%         axis tight image;
%         set(gca,'xtick',[], 'ytick',[]);
%         set(gca,'Ydir','normal');
%         subplot(122);
%         barh(-temp(:,1),1); hold on;
%         barh(temp(:,2),1);
    end

end


