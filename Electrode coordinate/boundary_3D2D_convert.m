clear all
% Data.E_posi_raw = elec_pos;
% EltdCollection = Data.E_posi_raw;
EltdCollection = elec_pos;
%% 3D to 2D interpolation
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

% OrderChecking = angle(complex(BoundPos(:,1),BoundPos(:,2)));
% [RightOrder,IndexOrder] = sort(OrderChecking-OrderChecking(1),'descend');
% FirstEltdIndex = find(RightOrder==0);
% ReorderBoundPos = BoundPos(IndexOrder([FirstEltdIndex:EltdNo 1:(FirstEltdIndex-1)]),:);
% Rx2 = [cos(OrderChecking(1)) sin(OrderChecking(1));-sin(OrderChecking(1)) cos(OrderChecking(1))];
% ReorderBoundPos = Rx2*ReorderBoundPos(:,1:2)';
% Pos = ReorderBoundPos';
%Pos2 = Pos(end:1);

% Pos = Pos([1 16:-1:2],:);
figure;
scatter(Pos(:,1),Pos(:,2),'filled');
% Pos(11:16,:) = Pos(7:12,:);
for i = 1:size(Pos,1)
    axis square; axis([-2 2 -2 2]);
     text(Pos(i,1),Pos(i,2),num2str(i), 'horizontal','left', 'vertical','bottom')
end
% for i = 7:10
%     i = 15
%     Pos(i,:) = ginput(1)'; hold on;
%     scatter(Pos(i,1),Pos(i,2),'filled');
% end

%% tile elec
% 5번 13번의 기울어진 각도를 계산하여 전체 전극위치를 보정함

temp = Pos; % 16 electrode point
axis square; scatter(temp(:,1),temp(:,2)); axis square; axis([-1.5 1.5 -1.5 1.5]);
for i = 1:size(temp,1)
    text(temp(i,1),temp(i,2), ['  ',num2str(i)]);
end
title('before');

tp1 = 5;
tp2 = 13;
rotAngle = -atan(temp(tp1,1)-temp(tp2,1))/(temp(tp2,2)-temp(tp1,2));

% tp1 = 1;
% tp2 = 10;
% rotAngle = -atan(temp(tp1,1)-temp(tp2,1))/(temp(tp2,2)-temp(tp1,2));
% rotAngle = rotAngle + 0.54*pi;

temp2(:,1) = temp(:,1).*cos(rotAngle) - temp(:,2).*sin(rotAngle);
temp2(:,2) = temp(:,1).*sin(rotAngle) + temp(:,2).*cos(rotAngle);

temp2(:,1) = -temp2(:,1);
figure;
axis square; scatter(temp2(:,1),temp2(:,2)); axis square; grid off; axis([-1.5 1.5 -1.5 1.5]);
for i = 1:size(temp2,1)
    text(temp2(i,1),temp2(i,2), ['  ',num2str(i)]);
end
title('after');

elec_pos = temp2;