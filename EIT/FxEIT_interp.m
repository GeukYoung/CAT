function [EIT_data2] = FxEIT_interp(EIT_data,intp_num)
for i = 1:size(EIT_data,1)
    EIT_data2(i,:) = interp1(1:size(EIT_data,2),EIT_data(i,:),1:1/intp_num:size(EIT_data,2));
end
