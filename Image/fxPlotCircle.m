function fxPlotCircle(X0, Y0, R, option_color)
% display
ang=0:0.01:2*pi;

% for j = 1:length(Real)
%     plot(Real(j),Quad(j),'r*');
%     hold on;
% end
xp=R*cos(ang);
yp=R*sin(ang);
% plot(X0,Y0,'b*');
hold on; 
if nargin == 3
    plot(X0+xp,Y0+yp);
elseif nargin == 4
   plot(X0+xp,Y0+yp,'Color',option_color,'linewidth',2);
end
%     axis([-1 1 -1 1]);
axis square;
