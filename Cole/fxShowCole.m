function fxShowCole(X0, Y0, R, option_color, width)
% display
ang=0:0.01:2*pi;

% for j = 1:length(Real)
%     plot(Real(j),Quad(j),'r*');
%     hold on;
% end
xp=R*cos(ang);
yp=R*sin(ang);
% plot(X0,Y0,'b*');
if nargin == 3
    plot(X0+xp,Y0+yp);
elseif nargin == 4
    plot(X0+xp,Y0+yp,'Color',option_color);
elseif nargin == 5
    plot(X0+xp,Y0+yp,'Color',option_color,'LineWidth',width);
end
%     axis([-1 1 -1 1]);
axis square;
