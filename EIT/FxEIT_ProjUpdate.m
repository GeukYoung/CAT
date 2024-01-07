function [Data2] = FxEIT_ProjUpdate(Data)
nData = size(Data,2);
Data2 = [];
for i = 1:nData-1
    temp = repmat(Data(:,i),1,16);
    temp2 = Data(:,i+1);
    for j = 1+1:16
        temp(1:13*(j-1),j) = temp2(1:13*(j-1));
    end
    Data2 = [Data2 temp];
end
