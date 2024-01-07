function [Voltage_data] = FxATR2_VImport_test(EIT_data_path)
meas_num = 256; % 256 measurement
Header_num = 0; % 0 byte header of 1 data
data_num = 2048;

cd(EIT_data_path);
cnt = 1;

file_path = strcat(EIT_data_path,'\','Voltage_0.bin');
fid = fopen(file_path);
raw_data = fread(fid,'double');
scan_num = size(raw_data,1)/meas_num;
Voltage_data = reshape(raw_data,meas_num,scan_num);


fclose(fid);

clear raw_data;