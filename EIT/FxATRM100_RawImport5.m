function [EIT_data] = FxATRM100_RawImport5(EIT_data_path)
scan_num = 512; % 1 binary data have 512 scan data
meas_num = 256; % 256 measurement
Header_num = 72; % 72 byte header of 1 data
data_num = 1536;

cd(EIT_data_path);
cnt = 1;

file_path = strcat(EIT_data_path,'\','EIT_0.bin');
fid = fopen(file_path);
raw_data = fread(fid,'uint8');
scan_num = size(raw_data,1)/(Header_num + data_num);
raw_data = reshape(raw_data,meas_num*6+Header_num,scan_num);


fclose(fid);
raw_data(1:72,:) = [];

for j = 1:scan_num
    temp = raw_data(:,j);
    temp = reshape(temp,6,256);
    temp(1,:)=temp(1,:)-128;
    temp(1:2,:) = temp([2 1],:);
    temp(3,:) = temp(4,:).*256 + temp(3,:); % VM용 L H 변경
    temp(4,:) = temp(6,:).*256 + temp(5,:); % VM용 L H 변경
    temp(5:6,:) = [];
    temp(temp>32767) = temp(temp>32767) - 65536;
    temp = temp';
    %         temp = temp([1:16:256 2:16:256 3:16:256 4:16:256 5:16:256 6:16:256 7:16:256 8:16:256 9:16:256 10:16:256 11:16:256 12:16:256 13:16:256 14:16:256 15:16:256 16:16:256],:);
    scan_data(:,:,cnt) = temp;
    EIT_data(:,cnt) = sqrt(scan_data(:,3,cnt).^2 + scan_data(:,4,cnt).^2);
    %         OV(:,cnt) = temp(:,2);
    clear temp;
    cnt = cnt + 1;
end
clear raw_data;
