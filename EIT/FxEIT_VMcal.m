function [Data_VMcal] = FxEIT_VMcal(Data,Data_R)
ch = 16;
Data_R = Data_R(:,9);
if length(Data_R) < 256
    idx_full = 1:ch^2;
    idx_full(FxEIT_mask(ch)) = [];
    temp_R = zeros(ch^2,size(Data_R,2));
    temp_R(idx_full,:) = Data_R;
    Data_R = temp_R
end

temp = Data_R
temp = reshape(temp,16,16);
VM_factor = zeros(ch,ch);
for j = 1:16
    for k = 1:16
%         VM_factor(j,k) = abs(((temp(j,k) - temp(k,j))/((temp(j,k) + temp(k,j))/2)) * 100);
        VM_factor(j,k) = (temp(j,k)+temp(k,j))/(2*temp(j,k));
    end
end
VM_factor = reshape(VM_factor,256,1);
VM_factor(FxEIT_mask(16)) = [];

if size(Data,1) == 256
    Data(FxEIT_mask(16),:) = [];
end

if size(Data_R,1) == 256
    Data_R(FxEIT_mask(16),:) = [];
end

% imdl= mk_common_model('d2d1c',16);
% img_1 = mk_image(imdl);
% stim = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
% img_1.fwd_model.stimulation = stim;
% vh = fwd_solve(img_1);
Data_VMcal = Data.*repmat(VM_factor,1,size(Data,2));

% Data_VMcal = Data.*repmat((vh.meas./Data_R(:,10)),1,size(Data,2));
