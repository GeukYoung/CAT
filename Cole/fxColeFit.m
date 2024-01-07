function [Rinf, R0 , X0, Y0, R , Error] = fxColeFit(Real, Quad)
nFreq = length(Real);
ang=0:0.01:2*pi;
% x0(1) = max(Real)+abs(mean(Real));
% y0(1) = min(Quad)-abs(mean(Quad));
x0(1) = (mean(Real));
y0(1) = (mean(Quad));
iter = 1000;

for i = 1:length(Real)
    dist(i) = pdist([x0(1) y0(1);Real(i) Quad(i)],'euclidean');
end
r(1) = min(dist);

error(1) = 0;
for j = 1:length(Real)
    error(1) = error(1) + (pdist([x0(1) y0(1);Real(j) Quad(j)])-r(1))^2;
end

for i = 1:iter-1
    % x,y,r update
    temp1 = 0; temp2 = 0; temp3 = 0;
    for j = 1:length(Real)
        temp1 = temp1 + (Real(j) - x0(i))/pdist([x0(i) y0(i);Real(j) Quad(j)]);
        temp2 = temp2 + (Quad(j) - y0(i))/pdist([x0(i) y0(i);Real(j) Quad(j)]);
        temp3 = temp3 + pdist([x0(i) y0(i);Real(j) Quad(j)]);
    end
    x0(i+1) = (1/length(Real))*(sum(Real) - r(i)*temp1);
    y0(i+1) = (1/length(Real))*(sum(Quad) - r(i)*temp2);
    r(i+1) = (1/length(Real))*temp3;
    
    % calcurate error
    error(i+1) = 0;
    for j = 1:length(Real)
        error(i+1) = error(i+1) + (pdist([x0(i+1) y0(i+1);Real(j) Quad(2)])-r(i+1))^2;
    end
    
    if (abs((error(end)-error(end-1))/error(end-1)) < 0.01) && (error(end) < 0.01)
        break;
    end
end

% % display
% subplot(121);
% for j = 1:length(Real)
%     plot(Real(j),Quad(j),'r*');
%     hold on;
% end
% xp=r(i+1)*cos(ang);
% yp=r(i+1)*sin(ang);
% plot(x0(i+1),y0(i+1),'b*');
% plot(x0(i+1)+xp,y0(i+1)+yp);
% %     axis([-1 1 -1 1]);
% axis square;
% text(-0.9,0.9,['Iter : ',num2str(i+1)]);
% text(-0.9,0.75,['Error : ',num2str(error(i+1))]);
% hold off;
% 
% subplot(122);
% plot(error,'LineWidth',2);
% axis([1 i+1 0 inf]);
% axis square;
% drawnow;
 
R0 = x0(end) + sqrt(r(end)^2 - y0(end)^2);
Rinf = x0(end) - sqrt(r(end)^2 - y0(end)^2);
alpha = 1 - (2/pi)*asin(abs(y0(end))/r(end));
Error = error(end);

X0 = x0(end);
Y0 = y0(end);
R = r(end);
end


% function [Rinf, R0 , x0, y0, r] = fxCoreFit(Real, Quad)
% nFreq = length(Real)
% ang=0:0.01:2*pi;
% x0(1) = max(Real)+abs(mean(Real));
% y0(1) = min(Quad)-abs(mean(Quad));
% iter = 1000;
% 
% for i = 1:length(Real)
%     dist(i) = pdist([x0(1) y0(1);Real(i) Quad(i)],'euclidean')
% end
% r(1) = min(dist);
% 
% error(1) = 0;
% for j = 1:length(Real)
%     error(1) = error(1) + (pdist([x0(1) y0(1);Real(j) Quad(j)])-r(1))^2;
% end
% 
% for i = 1:iter-1
%     % x,y,r update
%     temp1 = 0; temp2 = 0; temp3 = 0;
%     for j = 1:length(Real)
%         temp1 = temp1 + (Real(j) - x0(i))/pdist([x0(i) y0(i);Real(j) Quad(j)]);
%         temp2 = temp2 + (Quad(j) - y0(i))/pdist([x0(i) y0(i);Real(j) Quad(j)]);
%         temp3 = temp3 + pdist([x0(i) y0(i);Real(j) Quad(j)]);
%     end
%     x0(i+1) = (1/length(Real))*(sum(Real) - r(i)*temp1);
%     y0(i+1) = (1/length(Real))*(sum(Quad) - r(i)*temp2);
%     r(i+1) = (1/length(Real))*temp3;
%     
%     % calcurate error
%     error(i+1) = 0;
%     for j = 1:length(Real)
%         error(i+1) = error(i+1) + (pdist([x0(i+1) y0(i+1);Real(j) Quad(2)])-r(i+1))^2;
%     end
%     
% %     % display
% %     subplot(121);
% %     for j = 1:length(Real)
% %         plot(Real(j),Quad(j),'r*');
% %         hold on;
% %     end
% %     xp=r(i+1)*cos(ang);
% %     yp=r(i+1)*sin(ang);
% %     plot(x0(i+1),y0(i+1),'b*');
% %     plot(x0(i+1)+xp,y0(i+1)+yp);
% %     %     axis([-1 1 -1 1]);
% %     axis square;
% %     text(-0.9,0.9,['Iter : ',num2str(i+1)]);
% %     text(-0.9,0.75,['Error : ',num2str(error(i+1))]);
% %     hold off;
% %     
% %     subplot(122);
% %     plot(error,'LineWidth',2);
% %     axis([1 i+1 0 inf]);
% %     axis square;
% %     
% %     drawnow;
% %     %     pause(0.001);
% end
% R0 = x0(end) + sqrt(r(end)^2 - y0(end)^2);
% Rinf = x0(end) - sqrt(r(end)^2 - y0(end)^2);
% alpha = 1 - (2/pi)*asin(abs(y0(end))/r(end));
% Error = error(end);
% 
% end
