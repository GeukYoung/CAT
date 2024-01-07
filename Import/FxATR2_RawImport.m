function [EIT_data] = FxATR2_RawImport(EIT_data_path)
meas_num = 256; % 256 measurement
Header_num = 168; % 168 byte header of 1 data
data_num = 1536;

cd(EIT_data_path);
cnt = 1;
header_cnt=1;

file_path = strcat(EIT_data_path,'\','EIT_0.bin');
fid = fopen(file_path);
raw_data = fread(fid,'uint8');
scan_num = size(raw_data,1)/(Header_num + data_num);
raw_data = reshape(raw_data,meas_num*6+Header_num,scan_num);


fclose(fid);
header_data=raw_data(1:Header_num,:);
raw_data(1:Header_num,:) = [];

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
    EIT_data.EIT_raw(:,cnt) = sqrt(scan_data(:,3,cnt).^2 + scan_data(:,4,cnt).^2);
    
    temp2 = header_data(:,j);
    header_cnt=5;
    EIT_data.Amp(:,cnt)=typecast(uint8(temp2(header_cnt:header_cnt+7)),'double');  
    header_cnt=13;
    for k=1:16
        EIT_data.Gain(k,cnt) = typecast(uint8(temp2(header_cnt+8*(k-1):header_cnt+8*(k-1)+7)),'double');  
    end
    
    EIT_data.Timestamp(1,cnt)=typecast(uint8(temp2(151:152)),'uint16');    
    EIT_data.Timestamp(2,cnt)=temp2(153);
    EIT_data.Timestamp(3,cnt)=temp2(154);
    EIT_data.Timestamp(4,cnt)=temp2(155);
    EIT_data.Timestamp(5,cnt)=temp2(156);
    EIT_data.Timestamp(6,cnt)=temp2(157);
    EIT_data.Timestamp(7,cnt)=typecast(uint8(temp2(158:159)),'uint16'); 
    %         OV(:,cnt) = temp(:,2);
    clear temp;
    cnt = cnt + 1;
end
clear raw_data;