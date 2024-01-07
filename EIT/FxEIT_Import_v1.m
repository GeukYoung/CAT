function [Data] = FxEIT_Import_v1(Path_main)
cd(Path_main);
dirlist = dir('.');

cnt.eit = 0;
cnt.voltage = 0;
cnt.phase = 0;
cnt.pleth = 0;
for cnt_file = 1:length(dirlist)
    if contains(dirlist(cnt_file).name,'EIT_','IgnoreCase',true)
        if contains(dirlist(cnt_file).name,num2str(cnt.eit))
            cnt.eit = cnt.eit + 1;
            Path(cnt.eit).eit = strcat(Path_main, '\', dirlist(cnt_file).name);
        end
    end

    if contains(dirlist(cnt_file).name,'Voltage','IgnoreCase',true)
        if contains(dirlist(cnt_file).name,num2str(cnt.voltage))
            cnt.voltage = cnt.voltage + 1;
            Path(cnt.voltage).voltage = strcat(Path_main, '\', dirlist(cnt_file).name);
        end
    end
    
    if contains(dirlist(cnt_file).name,'Phase','IgnoreCase',true)
        if contains(dirlist(cnt_file).name,num2str(cnt.phase))
            cnt.phase = cnt.phase + 1;
            Path(cnt.phase).phase = strcat(Path_main, '\', dirlist(cnt_file).name);
        end
    end
    
    if contains(dirlist(cnt_file).name,'Pleth','IgnoreCase',true)
        if contains(dirlist(cnt_file).name,num2str(cnt.pleth))
            cnt.pleth = cnt.pleth + 1;
            Path(cnt.pleth).pleth = strcat(Path_main, '\', dirlist(cnt_file).name);
        end
    end
end
total_EIT = cnt.eit; clear cnt_file;

% check data version
fid = fopen(Path(1).eit);
version = fread(fid,4,'uint8');
version = version(4)*2^3 + version(3)*2^2 + version(2)*2^1 + version(1)*2^0;
fclose(fid);

Data.version = version;
Data.path = Path_main;
Data.mask_0skip = FxEIT_mask(16);
Data.mask_2skip = FxEIT_mask(16,[16 1 3 4]);

switch version
%% v0000 for Old version
    case 0 % 0000
        header.meas_num = 256; % 256 measurement
        header.header_num = 72; % 72 byte header of 1 data
        header.data_num = 1536;
        
        header.amp = 5;
        header.gain = 13;
        header.ci_flag = 53;
        header.as_flag = 54;
        header.time = 55;
        
        cnt_scan = 1;
        for cnt_EIT = 1:total_EIT
            % EIT raw
            fid = fopen(Path(cnt_EIT).eit);
            raw_data = fread(fid,'uint8');
            scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            raw_data = reshape(raw_data,header.meas_num*6+header.header_num,scan_num);
            fclose(fid);
            
            header_data=raw_data(1:header.header_num,:);
            raw_data(1:header.header_num,:) = [];
            
            for cnt_scan2 = 1:scan_num
                temp = raw_data(:,cnt_scan2);
                temp = reshape(temp,6,256);
                temp(1,:)=temp(1,:)-128;
                temp(1:2,:) = temp([2 1],:);
                temp(3,:) = temp(4,:).*256 + temp(3,:); % VM용 L H 변경
                temp(4,:) = temp(6,:).*256 + temp(5,:); % VM용 L H 변경
                temp(5:6,:) = [];
                temp(temp>32767) = temp(temp>32767) - 65536;
                temp = temp';
                Data.EIT_raw(:,cnt_scan) = sqrt(temp(:,3).^2 + temp(:,4).^2);
                temp2 = header_data(:,cnt_scan2);
                Data.Amp(:,cnt_scan)=typecast(uint8(temp2(header.amp:header.amp+7)),'double');
                
                for cnt_gain=1:4
                    Data.Gain(cnt_gain,cnt_scan) = typecast(uint8(temp2(header.gain+8*(cnt_gain-1):header.gain+8*(cnt_gain-1)+7)),'double');
                end
                
                Data.CI_flag(cnt_scan)=temp2(header.ci_flag);
                Data.AS_flag(cnt_scan)=temp2(header.as_flag);
                hms.Y(cnt_scan) = double(typecast(uint8(temp2(header.time:header.time+1)),'uint16'));
                hms.M(cnt_scan) = temp2(header.time+2);
                hms.D(cnt_scan) = temp2(header.time+3);
                hms.H(cnt_scan) = temp2(header.time+4);
                hms.MI(cnt_scan) = temp2(header.time+5);
                hms.S(cnt_scan) = temp2(header.time+6);
                hms.MS(cnt_scan) = double(typecast(uint8(temp2(header.time+7:header.time+8)),'uint16'));
                
                clear temp;
                cnt_scan = cnt_scan + 1;
            end
            clear raw_data;
            
            % voltage
            if cnt.voltage >= cnt_EIT
                fid = fopen(Path(cnt_EIT).voltage);
                raw_data = fread(fid,'double');
                scan_num = size(raw_data,1)/header.meas_num;
                if cnt_EIT == 1
                    Data.EIT_v = reshape(raw_data,header.meas_num,scan_num);
                else
                    Data.EIT_v = [Data.EIT_v reshape(raw_data,header.meas_num,scan_num)];
                end
            end
            
            % phase
            if cnt.voltage >= cnt_EIT
                fid = fopen(Path(cnt_EIT).phase);
                raw_data = fread(fid,'double');
                scan_num = size(raw_data,1)/header.meas_num;
                if cnt_EIT == 1
                    Data.EIT_phase = reshape(raw_data,header.meas_num,scan_num);
                else
                    Data.EIT_phase = [Data.EIT_phase reshape(raw_data,header.meas_num,scan_num)];
                end
            end
        end
        Data.t_hms = datetime(hms.Y, hms.M, hms.D, ...
            hms.H, hms.MI, hms.S, hms.MS);

%% v0010 for AirTom-R2
    case 2 % 0010
        header.meas_num = 256; % 256 measurement
        header.header_num = 168; % 168 byte header of 1 data
        header.data_num = 1536;
        
        header.amp = 5;
        header.gain = 13;
%         header.vsf = 45;
        header.ci_flag = 149;
        header.as_flag = 150;
        header.time = 151;
        
        cnt_scan = 1;
        for cnt_EIT = 1:total_EIT
            % EIT raw
            fid = fopen(Path(cnt_EIT).eit);
            raw_data = fread(fid,'uint8');
            scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            raw_data = reshape(raw_data,header.meas_num*6+header.header_num,scan_num);
            fclose(fid);
            
            header_data=raw_data(1:header.header_num,:);
            raw_data(1:header.header_num,:) = [];
            
            for cnt_scan2 = 1:scan_num
                temp = raw_data(:,cnt_scan2);
                temp = reshape(temp,6,256);
                temp(1,:)=temp(1,:)-128;
                temp(1:2,:) = temp([2 1],:);
                temp(3,:) = temp(4,:).*256 + temp(3,:); % VM용 L H 변경
                temp(4,:) = temp(6,:).*256 + temp(5,:); % VM용 L H 변경
                temp(5:6,:) = [];
                temp(temp>32767) = temp(temp>32767) - 65536;
                temp = temp';
                Data.EIT_raw(:,cnt_scan) = sqrt(temp(:,3).^2 + temp(:,4).^2);
                temp2 = header_data(:,cnt_scan2);
                Data.Amp(:,cnt_scan)=typecast(uint8(temp2(header.amp:header.amp+7)),'double');
                
                for cnt_gain=1:16
                    Data.Gain(cnt_gain,cnt_scan) = typecast(uint8(temp2(header.gain+8*(cnt_gain-1):header.gain+8*(cnt_gain-1)+7)),'double');
                end
                
                Data.CI_flag(cnt_scan)=temp2(header.ci_flag);
                Data.AS_flag(cnt_scan)=temp2(header.as_flag);
                hms.Y(cnt_scan) = double(typecast(uint8(temp2(header.time:header.time+1)),'uint16'));
                hms.M(cnt_scan) = temp2(header.time+2);
                hms.D(cnt_scan) = temp2(header.time+3);
                hms.H(cnt_scan) = temp2(header.time+4);
                hms.MI(cnt_scan) = temp2(header.time+5);
                hms.S(cnt_scan) = temp2(header.time+6);
                hms.MS(cnt_scan) = double(typecast(uint8(temp2(header.time+7:header.time+8)),'uint16'));
                
                clear temp;
                cnt_scan = cnt_scan + 1;
            end
            clear raw_data;
            
            % voltage
            if cnt.voltage >= cnt_EIT
                fid = fopen(Path(cnt_EIT).voltage);
                raw_data = fread(fid,'double');
                scan_num = size(raw_data,1)/header.meas_num;
                if cnt_EIT == 1
                    Data.EIT_v = reshape(raw_data,header.meas_num,scan_num);
                else
                    Data.EIT_v = [Data.EIT_v reshape(raw_data,header.meas_num,scan_num)];
                end
            end
            
            % phase
            if cnt.voltage >= cnt_EIT
                fid = fopen(Path(cnt_EIT).phase);
                raw_data = fread(fid,'double');
                scan_num = size(raw_data,1)/header.meas_num;
                if cnt_EIT == 1
                    Data.EIT_phase = reshape(raw_data,header.meas_num,scan_num);
                else
                    Data.EIT_phase = [Data.EIT_phase reshape(raw_data,header.meas_num,scan_num)];
                end
            end
        end
        Data.t_hms = datetime(hms.Y, hms.M, hms.D, ...
            hms.H, hms.MI, hms.S, hms.MS);
        
%% v0011 for AirTom-PV/HV
    case 3 % 0011
        header.meas_num = 256; % 256 measurement
        header.header_num = 168; % 168 byte header of 1 data
        header.data_num = 1536*2; % add filtered data
        
        header.amp = 5;
        header.gain = 13;
        header.idx_scan = 141;
        header.ci_flag = 149;
        header.as_flag = 150;
        header.time = 151+1;
        
        cnt_scan = 1;
        for cnt_EIT = 1:total_EIT
            % EIT raw
            fid = fopen(Path(cnt_EIT).eit);
            raw_data = fread(fid,'uint8');
            scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            raw_data = reshape(raw_data,(header.header_num + header.data_num),scan_num);
            fclose(fid);
            
            header_data=raw_data(1:header.header_num,:);
            raw_data(1:header.header_num,:) = [];
            
            for cnt_scan2 = 1:scan_num
                temp = raw_data(:,cnt_scan2);
                temp = reshape(temp,6,256*2);
                temp(1,:)=temp(1,:)-128;
                temp(1:2,:) = temp([2 1],:);
                temp(3,:) = temp(4,:).*256 + temp(3,:); % VM용 L H 변경
                temp(4,:) = temp(6,:).*256 + temp(5,:); % VM용 L H 변경
                temp(5:6,:) = [];
                temp(temp>32767) = temp(temp>32767) - 65536;
                temp = temp';
                Data.EIT_raw(:,cnt_scan) = sqrt(temp(1:256,3).^2 + temp(1:256,4).^2);
                Data.EIT_raw_filt(:,cnt_scan) = sqrt(temp(257:end,3).^2 + temp(257:end,4).^2);
                
                temp2 = header_data(:,cnt_scan2);
                Data.Amp(:,cnt_scan)=typecast(uint8(temp2(header.amp:header.amp+7)),'double');
                
                for cnt_gain=1:16
                    Data.Gain(cnt_gain,cnt_scan) = typecast(uint8(temp2(header.gain+8*(cnt_gain-1):header.gain+8*(cnt_gain-1)+7)),'double');
                end
                Data.CI_flag(cnt_scan)=temp2(header.ci_flag);
                Data.AS_flag(cnt_scan)=temp2(header.as_flag);
                epochtime(cnt_scan) = typecast(uint8(temp2(header.time:header.time+7)),'uint64');

                clear temp;
                cnt_scan = cnt_scan + 1;
            end
            clear raw_data;
            
            % voltage
            if cnt.voltage >= cnt_EIT
                fid = fopen(Path(cnt_EIT).voltage);
                raw_data = fread(fid,'double');
                scan_num = size(raw_data,1)/header.meas_num;
                if cnt_EIT == 1
                    Data.EIT_v = reshape(raw_data,header.meas_num,scan_num);
                else
                    Data.EIT_v = [Data.EIT_v reshape(raw_data,header.meas_num,scan_num)];
                end
            end
            
            % phase
            if cnt.voltage >= cnt_EIT
                fid = fopen(Path(cnt_EIT).phase);
                raw_data = fread(fid,'double');
                scan_num = size(raw_data,1)/header.meas_num;
                if cnt_EIT == 1
                    Data.EIT_phase = reshape(raw_data,header.meas_num,scan_num);
                else
                    Data.EIT_phase = [Data.EIT_phase reshape(raw_data,header.meas_num,scan_num)];
                end
            end
            
            % pleth
            if cnt.Pleth >= cnt_EIT
                fid = fopen(Path(cnt_EIT).pleth);
                raw_data = fread(fid,'uint16');
                scan_num = size(raw_data,1)/9;
                if cnt_EIT == 1
                    Data.EIT_pleth = reshape(raw_data,9,scan_num);
                else
                    Data.EIT_pleth = [Data.EIT_phase reshape(raw_data,header.meas_num,scan_num)];
                end
            end
        end
        Data.t_hms = datetime(epochtime,'ConvertFrom','epochtime','TicksPerSecond',1e3);
        
%% v1000 for AirTom-RM
    case 8 % 1000
        header.meas_num = 256; % 256 measurement
        header.header_num = 72; % 72 byte header of 1 data
        header.data_num = 1536;
        
        header.amp = 5;
        header.gain = 13;
        header.vsf = 45;
        header.ci_flag = 53;
        header.as_flag = 54;
        header.time = 55;
        
        cnt_scan = 1;
        for cnt_EIT = 1:total_EIT
            % EIT raw
            fid = fopen(Path(cnt_EIT).eit);
            raw_data = fread(fid,'uint8');
            scan_num = size(raw_data,1)/(header.header_num + header.data_num);
            raw_data = reshape(raw_data,header.meas_num*6+header.header_num,scan_num);
            fclose(fid);
            
            header_data=raw_data(1:header.header_num,:);
            raw_data(1:header.header_num,:) = [];
            
            for cnt_scan2 = 1:scan_num
                temp = raw_data(:,cnt_scan2);
                temp = reshape(temp,6,256);
                temp(1,:)=temp(1,:)-128;
                temp(1:2,:) = temp([2 1],:);
                temp(3,:) = temp(4,:).*256 + temp(3,:); % VM용 L H 변경
                temp(4,:) = temp(6,:).*256 + temp(5,:); % VM용 L H 변경
                temp(5:6,:) = [];
                temp(temp>32767) = temp(temp>32767) - 65536;
                temp = temp';
                Data.EIT_raw(:,cnt_scan) = sqrt(temp(:,3).^2 + temp(:,4).^2);
                temp2 = header_data(:,cnt_scan2);
                Data.Amp(:,cnt_scan)=typecast(uint8(temp2(header.amp:header.amp+7)),'double');
                
                for cnt_gain=1:4
                    Data.Gain(cnt_gain,cnt_scan) = typecast(uint8(temp2(header.gain+8*(cnt_gain-1):header.gain+8*(cnt_gain-1)+7)),'double');
                end
                
                Data.VSF(cnt_scan)=typecast(uint8(temp2(header.vsf:header.vsf+7)),'double');
                Data.CI_flag(cnt_scan)=temp2(header.ci_flag);
                Data.AS_flag(cnt_scan)=temp2(header.as_flag);
                hms.Y(cnt_scan) = double(typecast(uint8(temp2(header.time:header.time+1)),'uint16'));
                hms.M(cnt_scan) = temp2(header.time+2);
                hms.D(cnt_scan) = temp2(header.time+3);
                hms.H(cnt_scan) = temp2(header.time+4);
                hms.MI(cnt_scan) = temp2(header.time+5);
                hms.S(cnt_scan) = temp2(header.time+6);
                hms.MS(cnt_scan) = double(typecast(uint8(temp2(header.time+7:header.time+8)),'uint16'));
                
                clear temp;
                cnt_scan = cnt_scan + 1;
            end
            clear raw_data;
            
            % voltage
            if cnt.voltage >= cnt_EIT
                fid = fopen(Path(cnt_EIT).voltage);
                raw_data = fread(fid,'double');
                scan_num = size(raw_data,1)/header.meas_num;
                if cnt_EIT == 1
                    Data.EIT_v = reshape(raw_data,header.meas_num,scan_num);
                else
                    Data.EIT_v = [Data.EIT_v reshape(raw_data,header.meas_num,scan_num)];
                end
            end
            
%             % phase
%             if cnt.voltage >= cnt_EIT
%                 fid = fopen(Path(cnt_EIT).phase);
%                 raw_data = fread(fid,'double');
%                 scan_num = size(raw_data,1)/header.meas_num;
%                 if cnt_EIT == 1
%                     Data.EIT_phase = reshape(raw_data,header.meas_num,scan_num);
%                 else
%                     Data.EIT_phase = [Data.EIT_phase reshape(raw_data,header.meas_num,scan_num)];
%                 end
%             end
        end
        Data.t_hms = datetime(hms.Y, hms.M, hms.D, ...
            hms.H, hms.MI, hms.S, hms.MS);       
end
end

function [mask] = FxEIT_mask(ch,first_sat)
    temp = zeros(ch,ch);
    if nargin == 1
        first_sat = [ch,1,2];
    end

    if length(first_sat) == 4
        cnt1 = first_sat(1); cnt2 = first_sat(2); cnt3 = first_sat(3); cnt4 = first_sat(4);
        for i = 1:ch
            temp([cnt1,cnt2,cnt3,cnt4],i) = 1;
            cnt1 = cnt1 + 1;
            cnt2 = cnt2 + 1;
            cnt3 = cnt3 + 1;
            cnt4 = cnt4 + 1;
            if cnt1 > ch
                cnt1 = 1;
            end
            if cnt2 > ch
                cnt2 = 1;
            end
            if cnt3 > ch
                cnt3 = 1;
            end
            if cnt4 > ch
                cnt4 = 1;
            end
        end
    else
        cnt1 = first_sat(1); cnt2 = first_sat(2); cnt3 = first_sat(3);
        for i = 1:ch
            temp([cnt1,cnt2,cnt3],i) = 1;
            cnt1 = cnt1 + 1;
            cnt2 = cnt2 + 1;
            cnt3 = cnt3 + 1;
            if cnt1 > ch
                cnt1 = 1;
            end
            if cnt2 > ch
                cnt2 = 1;
            end
            if cnt3 > ch
                cnt3 = 1;
            end
        end
    end

    temp2 = reshape(temp,ch*ch,1);
    [mask, ~] = find(temp2 == 1);
end