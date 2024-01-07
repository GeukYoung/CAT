function [ vout, dataout ] = ROBIN_Parser_multisetting( filename )
% filename: full file path (e.g.'C:\work\ROBIN\Data\ROBIN_01.TXT')
    
python_call = ['python ROBIN_parse.py' ' ' filename];
system(python_call);

folder_path = filename(1:length(filename)-4);

BioZfilename = [folder_path '/Bioz.csv'];
BioZSettingFilename = [folder_path '/Freq.csv'];

measSetting = csvread(BioZSettingFilename);

vo = dlmread(BioZfilename, ' ');
vi = vo(:,2); vq = vo(:,1);
vout = [vi, vq];

v2 = (floor(measSetting/(2^28)));
v1 = (mod(floor(measSetting/(2^24)),16));
i2 = (mod(floor(measSetting/(2^20)),16));
i1 = (mod(floor(measSetting/(2^16)),16));
freq = mod(measSetting, 256);

% freq = measSetting(:,1);
% i2 = measSetting(:,2);
% i1 = measSetting(:,3);
% v2 = measSetting(:,4);
% v1 = measSetting(:,5);

n = 16;
rfreq = repmat(freq,1,n)';
rfreq = rfreq(:);

ri2 = repmat(i2,1,n)';
ri2 = ri2(:);

ri1 = repmat(i1,1,n)';
ri1 = ri1(:);

rv2 = repmat(v2,1,n)';
rv2 = rv2(:);

rv1 = repmat(v1,1,n)';
rv1 = rv1(:);

setting = [abs(rfreq), abs(ri1), abs(ri2), abs(rv1), abs(rv2)];

offset = -4;

crop_setting = setting(1:end+offset,:);
crop_vout = vout(1-offset:end,:);

dataout = [crop_setting, crop_vout];
csvwrite(strcat(folder_path,'/Bioz_processed.csv'),dataout);
% data format: frequency, I1, I2, V1, V2, Vi-value, Vq-value

delete(BioZfilename);
delete(BioZSettingFilename);

