function [Theta] = FxIMU_Filt(Acc,Gyro,fs,ani)
if size(Acc,2) == 3
    Acc = Acc';
end
if size(Gyro,2) == 3
    Gyro = Gyro';
end

time = (1:length(Acc))'/fs;
ts = 1/fs;
% Acc = Acc / norm(Acc);
Gyro = Gyro./5;
%% Cal Angle
% Angular velocity to Angle
Gyroscope2Angle = (cumsum(Gyro')/250)';

% Acc to Euler Angle
AccAng(1,:) = rad2deg(atan(Acc(2,:)./sqrt(Acc(1,:).^2 + Acc(3,:).^2)));
AccAng(2,:) = rad2deg(atan(-Acc(1,:)./sqrt(Acc(2,:).^2 + Acc(3,:).^2)));
AccAng(3,:) = rad2deg(atan(deg2rad(sqrt(Acc(1,:).^2 + Acc(2,:).^2)./Acc(3,:))));

% Complimentaly filter
for j = 1:3
    temp_AccAng = AccAng(j,:);
    temp_Gyro = Gyro(j,:);
    StartPoint = 1;
    a = 0.1;
    ThetaOut = zeros(StartPoint,1);
    int_temp = 0;
    for i = (StartPoint + 1):length(temp_AccAng)
        temp = -1/a*(-temp_AccAng(i) + ThetaOut(i-1)) + temp_Gyro(i);
        ThetaOut(i) = ThetaOut(i-1) + temp*ts;
    end
    Theta(j,:) = ThetaOut;
end

%% Plot result
figure('Name', 'Sensor Data'); set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
subplot(3,2,2);
% axis(1) = subplot(3,1,1);
hold on;
plot(time, Gyro(1,:), 'r');
plot(time, Gyro(2,:), 'g');
plot(time, Gyro(3,:), 'b');
legend('X', 'Y', 'Z','Location','northwest'); xlim([0 120]);
xlabel('Time (s)');
ylabel('Angular rate (deg/s)');
title('Gyroscope'); 
xlim([-inf inf]);
hold off;

subplot(3,2,1);
% axis(2) = subplot(3,1,2);
hold on;
plot(time, Acc(1,:), 'r');
plot(time, Acc(2,:), 'g');
plot(time, Acc(3,:), 'b');
legend('X', 'Y', 'Z','Location','northwest'); xlim([0 120]);
xlabel('Time (s)');
ylabel('Acceleration (g)');
title('Accelerometer');
xlim([-inf inf]);
hold off;

subplot(3,2,4);
hold on;
plot(time, Gyroscope2Angle(1,:), 'r');
plot(time, Gyroscope2Angle(2,:), 'g');
plot(time, Gyroscope2Angle(3,:), 'b');
legend('X', 'Y', 'Z','Location','northwest'); xlim([0 120]);
xlabel('Time (s)');
ylabel('Gyro Angle (degree)');
title('Gyro Angle');
xlim([-inf inf]);
hold off;

subplot(3,2,3);
hold on;
plot(time, AccAng(1,:), 'r');
plot(time, AccAng(2,:), 'g');
plot(time, AccAng(3,:), 'b');
legend('X', 'Y', 'Z','Location','northwest'); xlim([0 120]);
xlabel('Time (s)');
ylabel('Acc Angle (degree)');
title('Acc Angle');
xlim([-inf inf]);
hold off;

subplot(3,2,5);
hold on;
plot(time, Theta(1,:), 'r');
plot(time, Theta(2,:), 'g');
plot(time, Theta(3,:), 'b');
legend('X', 'Y', 'Z','Location','northwest'); xlim([0 120]);
xlabel('Time (s)');
ylabel('Filtout Angle (degree)');
title('Filtout Angle');
xlim([-inf inf]);
hold off;

%% Animate
if nargin > 3
    figure; set(gcf,'units','normalized','outerposition',[0.1 0.05 0.45 0.9]);
%     subplot(1,3,[2 3]);
%     subplot(3,2,[1 2 3 4]);
%     subplot(3,2,[5 6]);
    subplot(3,3,7:9);
%     plot(1:length(Theta),Theta(1,:),1:length(Theta),Theta(2,:),1:length(Theta),Theta(3,:)); hold on;
    plot(1:length(Theta),Theta(1,:),'r'); hold on;
    plot(1:length(Theta),Theta(2,:),'g');
    plot(1:length(Theta),Theta(3,:),'b');
    temp_bar = plot([1 1],[-300 300],'r');
    legend('roll','pitch','yaw','Location','northeast');
%     xlim([2000 9000])
    ylim([-90 90]); xlim([-inf inf])
    set(gca,'ytick',[-90 0 90]);
    set(gca,'xtick',[]);
    set(gca,'YTickLabel',{'-\pi/2', '0' ,'\pi/2'});
    title('Angle data','fontsize',15);
    % set(gca,'FontSize',12,'Font','symbol');
    
%     subplot(3,2,[1 2 3 4]);
    subplot(3,3,1:6);
    hold on;
    view(-6.7,35);
    axis equal;
    hold on; grid on;
    set(gca,'Xlim',[-2 2]);
    set(gca,'Ylim',[-2 2]);
    set(gca,'Zlim',[-2 2]);
    % set(gca,'XDir','reverse');
    set(gcf,'renderer','opengl');
    xlabel('x'); ylabel('y'); zlabel('z');
    at1 = [0; 0; 2.5]*4;
    at2 = [0; 2.5; 0]*4;
    at3 = [2.5; 0; 0]*4;
    st =  [0; 0; 0]*3;
    pat1 = plot3([st(1) at1(1)],[st(2) at1(2)], [st(3) at1(3)],'b--','linewidth',0.1);
    pat2 = plot3([st(1) at2(1)],[st(2) at2(2)], [st(3) at2(3)],'g--','linewidth',0.1);
    pat3 = plot3([st(1) at3(1)],[st(2) at3(2)], [st(3) at3(3)],'r--','linewidth',0.1);
    
    an1 = [0; 0; 10;];
    an2 = [0; 10; 0;];
    an3 = [10; 0; 0;];
    sz =  [0; 0; 0;];
    plot3([sz(1) an1(1)],[sz(2) an1(2)], [sz(3) an1(3)],'b-','linewidth',0.1);
    plot3([sz(1) an2(1)],[sz(2) an2(2)], [sz(3) an2(3)],'g-','linewidth',0.1);
    plot3([sz(1) an3(1)],[sz(2) an3(2)], [sz(3) an3(3)],'r-','linewidth',0.1);
    an1 = [0; 0; -10;];
    an2 = [0; -10; 0;];
    an3 = [-10; 0; 0;];
    sz =  [0; 0; 0;];
    plot3([sz(1) an1(1)],[sz(2) an1(2)], [sz(3) an1(3)],'b-','linewidth',0.1);
    plot3([sz(1) an2(1)],[sz(2) an2(2)], [sz(3) an2(3)],'g-','linewidth',0.1);
    plot3([sz(1) an3(1)],[sz(2) an3(2)], [sz(3) an3(3)],'r-','linewidth',0.1);
    
    Vert = [0 0 0.25; 1 0 0.25; 1 1 0.25; 0 1 0.25; 0 0 0.75; 1 0 0.75; 1 1 0.75; 0 1 0.75];
    % Vert = -Vert;
    Vert(:,1) = (Vert(:,1)-0.5);
    Vert(:,2) = (Vert(:,2)-0.5);
    Vert(:,3) = (Vert(:,3)-0.5);
    % Vert = -Vert;
    Faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
    ptch.Vertices = Vert;
    ptch.Faces = Faces;
    ptch.FaceVertexCData = copper(6);
    ptch.FaceColor = 'flat';
    patch_handle = patch(ptch);
    title('Animation','fontsize',15);
    % 320600:10:321000
    
    set(gcf,'color','white');
    cnt = 1;

    for i = 1:round(length(Theta)/1000):length(Theta)
        set(temp_bar,'Xdata',[i i]);
        
%         temp = [222300 268700 299600 516500]
%         i = temp(4)
        sz = sin(deg2rad(0));
        cz = cos(deg2rad(0));
        sy = sin(deg2rad(Theta(2,i)));
        cy = cos(deg2rad(Theta(2,i)));
        sx = sin(deg2rad(Theta(1,i)));
        cx = cos(deg2rad(Theta(1,i)));
        
        Cx = [1 0 0; 0 cx sx; 0 -sx cx];
        Cy = [cy 0 -sy; 0 1 0; sy 0 cy];
        Cz = [cz sz 0; -sz cz 0; 0 0 1];
        
        Cbn = Cx*Cy*Cz;
        
        s_ = st;
        a1_ = Cbn'*at1;
        a2_ = Cbn'*at2;
        a3_ = Cbn'*at3;
        
        sz = sin(deg2rad(0));
        cz = cos(deg2rad(0));
        sy = sin(deg2rad(-Theta(2,i)));
        cy = cos(deg2rad(-Theta(2,i)));
        sx = sin(deg2rad(-Theta(1,i)));
        cx = cos(deg2rad(-Theta(1,i)));
        
        Cx = [1 0 0; 0 cx sx; 0 -sx cx];
        Cy = [cy 0 -sy; 0 1 0; sy 0 cy];
        Cz = [cz sz 0; -sz cz 0; 0 0 1];
        
        Cbn = Cx*Cy*Cz;
        
        Vert_ = Vert;
        for j=1:size(Vert,1)
            Vert_(j,:) = (Vert(j,:)*Cbn');
        end
        set(patch_handle,'Vertices',Vert_);
        
        set(pat1,'Xdata',[s_(1) a1_(1)]);
        set(pat1,'Ydata',[s_(2) a1_(2)]);
        set(pat1,'Zdata',[s_(3) a1_(3)]);
        set(pat2,'Xdata',[s_(1) a2_(1)]);
        set(pat2,'Ydata',[s_(2) a2_(2)]);
        set(pat2,'Zdata',[s_(3) a2_(3)]);
        set(pat3,'Xdata',[s_(1) a3_(1)]);
        set(pat3,'Ydata',[s_(2) a3_(2)]);
        set(pat3,'Zdata',[s_(3) a3_(3)]);
        drawnow;
        disp(num2str(i));
        
        cnt = cnt + 1;
    end
end
end