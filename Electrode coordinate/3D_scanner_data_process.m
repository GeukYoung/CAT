[path_1, path_2] = uigetfile('.ply');
Temp = pcread(strcat(path_2,path_1));
% Temp = ply_read(strcat(path_2,path_1));

Raw_Data.P = double(Temp.Location);
Raw_Data.Color = Temp.Color;

clear BT_elec idx_Elec h t hFig; close all;
% figure;
hFig = figure;
showPointCloud(Raw_Data.P,Raw_Data.Color,'MarkerSize',100); xlabel('x'); ylabel('y'); zlabel('z'); hold on;
BT_tag = {'E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12','E13','E14','E15','E16',};
for i = 1:16
    BT_elec(16-i+1) = uicontrol(hFig, 'Style', 'togglebutton','units','normalized', ...
        'String', BT_tag(16-i+1),'position',[0 1/16*(i-1) 1.5/16 1/16]);
end
% BT_stop = uicontrol(hFig, 'Style', 'pushbutton', 'String', 'Stop','position',[10 10+50*(17-1) 100 50]);

%% Elec pos click
% clickA3DPoint(Raw_Data.P',Raw_Data.Color');
for cnt = 1:16
    En_idxBT(cnt) = BT_elec(cnt).Value;
end

if sum(En_idxBT) > 1
    errordlg('Please select only 1 button','Multi select');
elseif sum(En_idxBT) == 0
    errordlg('Please select elec ch','None select');
else
    En_idxBT = En_idxBT.*[1:16];
    En_idxBT = sum(En_idxBT);
    
    temp = callbackClickA3DPoint([],[],Raw_Data.P');
%     if exist(h{En_idxBT}) == []
% %         exist(strcat('h{',num2str(En_idxBT),'}.Visible')) == 1
%         delete(h{En_idxBT});
%         delete(t{En_idxBT});
%     end
    h{En_idxBT} = plot3(Raw_Data.P(temp,1),Raw_Data.P(temp,2),Raw_Data.P(temp,3),'ro','MarkerFaceColor','r','markersize',10);
    t{En_idxBT} = text(Raw_Data.P(temp,1),Raw_Data.P(temp,2),Raw_Data.P(temp,3),strcat('E',num2str(En_idxBT)));
    idx_Elec(En_idxBT,1) = temp;
    idx_Elec(En_idxBT,2) = 1;
end

for cnt = 1:16
    BT_elec(cnt).Value = 0;
end

%%
sum_elec_pos = zeros(1,3);
for cnt = 1:16
    if idx_Elec(cnt,2) == 1
        elec_pos(cnt,:) = Raw_Data.P(idx_Elec(cnt,1),:);
        sum_elec_pos = sum_elec_pos + elec_pos(cnt,:);
    end
end
mean_elec_pos = sum_elec_pos/sum(idx_Elec(:,2));
clear sum_elec_pos;

for cnt = 1:16
    if idx_Elec(cnt,2) == 0
        elec_pos(cnt,:) = mean_elec_pos;
    end
end

plot3(elec_pos(:,1),elec_pos(:,2),elec_pos(:,3),'o')

figure;
for i = 1:16
    hold on; plot3(elec_pos(i,1),elec_pos(i,2),elec_pos(i,3),'ro','Markersize',15);
    text(elec_pos(i,1),elec_pos(i,2),elec_pos(i,3),num2str(i));
end
axis square equal;

%% 3D to 2D interpolation
EltdCollection = elec_pos;
EltdNo = size(EltdCollection,1); 
Center = mean(EltdCollection);
%%%%%%%%%%%   PlaneInfor(1)*x +PlaneInfor(2)*y + PlaneInfor(3) = z
PlaneInfor = [EltdCollection(:,1:2) ones(EltdNo,1)]\EltdCollection(:,3);
NormalPlane = [PlaneInfor(1:2); -1]/norm([PlaneInfor(1:2); -1]);
if NormalPlane(3)<0
    NormalPlane = - NormalPlane;   %% NormalPlane(1)*x+NormalPlane(2)*y+NormalPlane(3)*z+d =0;    
end
d = - sum(NormalPlane'.*Center);
NormalCos = NormalPlane(2)/norm(NormalPlane(1:2));
NormalSin = NormalPlane(1)/norm(NormalPlane(1:2));
Rz = [NormalCos -NormalSin 0; NormalSin NormalCos  0; 0 0 1];
Rx = [1 0 0; 0 NormalPlane(3) -sqrt(1-NormalPlane(3)^2); 0  sqrt(1-NormalPlane(3)^2) NormalPlane(3)];
EltdZLevel = Rx*Rz*(EltdCollection'-repmat(Center',[1,EltdNo]));
EltdZLevel = EltdZLevel';
ScaleSize = max(EltdZLevel(:,1));
EltdZLevel = 1/ScaleSize*EltdZLevel;
figure,plot3(EltdZLevel(:,1),EltdZLevel(:,2),EltdZLevel(:,3),'r*');view(2); axis equal
Pos = EltdZLevel;
Pos(:,3) = [];

for i = 1:size(Pos,1)
     axis square; axis([-2 2 -2 2]);
     text(EltdZLevel(i,1),EltdZLevel(i,2),num2str(i), 'horizontal','left', 'vertical','bottom')
end

hFig = figure;
h = scatter(Pos(:,1),Pos(:,2),'filled');
% Pos(11:16,:) = Pos(7:12,:);
for i = 1:size(Pos,1)
    axis square; axis([-2 2 -2 2]);
     text(Pos(i,1),Pos(i,2),num2str(i), 'horizontal','left', 'vertical','bottom')
end

%% Point confirm
clear h t hFig; close all;
hFig = figure;
hold on;
for cnt = 1:16
    if (cnt == 1) || (cnt == 5) || (cnt == 9) || (cnt == 13)
        h{cnt} = plot(Pos(cnt,1),Pos(cnt,2),'ro','MarkerFaceColor','r');
    else
        h{cnt} = plot(Pos(cnt,1),Pos(cnt,2),'bo','MarkerFaceColor','b');
    end
    t{cnt} = text(Pos(cnt,1),Pos(cnt,2),0.05,strcat(num2str(cnt)));
end
axis square; axis([-1.5 1.5 -1.5 1.5]);

En_Elec = ~idx_Elec(:,2)' .* [1:16];
En_Elec(En_Elec == 0) = [];

BT_tag = {'E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12','E13','E14','E15','E16',};
for i = En_Elec
    BT_elec(i) = uicontrol(hFig, 'Style', 'togglebutton','units','normalized', ...
        'String', BT_tag(i),'position',[0 1/16*(16-(i-1)) 1.5/16 1/16]);
end
% BT_stop = uicontrol(hFig, 'Style', 'pushbutton', 'String', 'Stop','position',[10 10+50*(17-1) 100 50]);

%% Elec pos click
% clickA3DPoint(Raw_Data.P',Raw_Data.Color');
clear En_idxBT;
for cnt = En_Elec
    En_idxBT(cnt) = BT_elec(cnt).Value;
end

if sum(En_idxBT) > 1
    errordlg('Please select only 1 button','Multi select');
elseif sum(En_idxBT) == 0
    errordlg('Please select elec ch','None select');
else
    En_idxBT = En_idxBT.*[1:length(En_idxBT)];
    En_idxBT = sum(En_idxBT);
    
    temp = ginput(1);
    
    h{En_idxBT}.XData = temp(1);
    h{En_idxBT}.YData = temp(2);
    t{En_idxBT}.Position = [temp 0];
    Pos(En_idxBT,:) = temp;
end

for cnt = En_Elec
    BT_elec(cnt).Value = 0;
end

%%
figure;
scatter(Pos(:,1),Pos(:,2),'filled');
% Pos(11:16,:) = Pos(7:12,:);
for i = 1:size(Pos,1)
    axis square; axis([-2 2 -2 2]);
     text(Pos(i,1),Pos(i,2),num2str(i), 'horizontal','left', 'vertical','bottom')
end

% tile elec
% 1번 9번의 기울어진 각도를 계산하여 전체 전극위치를 보정함
temp = Pos; % 16 electrode point
axis square; scatter(temp(:,1),temp(:,2)); axis square; axis([-1.5 1.5 -1.5 1.5]);
for i = 1:size(temp,1)
    text(temp(i,1),temp(i,2), ['  ',num2str(i)]);
end
title('before');

tp1 = 1;
tp2 = 9;
rotAngle = -atan((temp(tp1,2)-temp(tp2,2))/(temp(tp1,1)-temp(tp2,1))); rad2deg(rotAngle)
% rotAngle = rotAngle-0.5*pi;
% tp1 = 1;
% tp2 = 10;
% rotAngle = -atan(temp(tp1,1)-temp(tp2,1))/(temp(tp2,2)-temp(tp1,2));
% rotAngle = rotAngle + 0.54*pi;

temp2(:,1) = temp(:,1).*cos(rotAngle) - temp(:,2).*sin(rotAngle);
temp2(:,2) = temp(:,1).*sin(rotAngle) + temp(:,2).*cos(rotAngle);

temp2(:,1) = -temp2(:,1);
temp2(:,1) = temp2(:,1) - (max(temp2(:,1))+min(temp2(:,1)))/2;
temp2(:,2) = temp2(:,2) - (max(temp2(:,2))+min(temp2(:,2)))/2;
temp2 = temp2.*1/max(max(abs(temp2)))
figure;
axis square; scatter(temp2(:,1),temp2(:,2)); axis square; grid off; axis([-1.5 1.5 -1.5 1.5]);
for i = 1:size(temp2,1)
    text(temp2(i,1),temp2(i,2), ['  ',num2str(i)]);
end
title('after');

elec_pos2 = temp2;