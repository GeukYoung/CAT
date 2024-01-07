function [Zout1 Zout2 Zout3 Zout4 result] = Fx_Call_CCScalfile(folder_path,freq_num,ch)
% Function info
% [Zout1 Zout2 Zout3 Zout4] = Fx_Call_CCScalfile(forler_path,freq,ch)
% input files -> 3 nargin
%       [folder path, freq, IMM ch]
%        folder path (str) -> ex) D:\Dropbox\#Lab Work\1. EIT_System\1. 16ch EIT\[FINAL_SYDNEY]EIT_Mark25_20120128\Debug\Calibration\eit1
%        freq (int) -> 1 : 1.125kHz
%                      2 : 11.25KHz
%                      3 : 56.25kHz
%                      4 : 112.5KHz
%        IMM ch (int) -> ex) 1~16
% 
% output files -> 5 matrix out
%       Zout1 : coarse 10 step(26x26) 
%       Zout2 : coarse  1 step(21x21) 
%       Zout3 : coarse 10 step(26x26)
%       Zout4 : coarse  1 step(21x21)
%       result : [CoarseR CoarseC FineR FineC Zout]

freq_index = {'1.125kHz','12.5KHz','62.5kHz','125kHz'};
freq = freq_index{freq_num};

switch freq_num
    case 1 % 1kHz
        file_path = strcat(folder_path,'\OutputImpedance\',freq,'\Log\Rout\Rout',int2str(ch),'Ch.txt');
        fid = fopen(file_path);
        temp=fscanf(fid,'%s',1);
        
        %% coarse
        for i = 1:26
            temp = fscanf(fid,'%*s%g',1);
            Zout1(i,:) = temp;
            clear temp;
        end
        temp = fscanf(fid,'%*s%g%*s%g%*s%g',6);
        
        for i = 1:21
            temp = fscanf(fid,'%*s%g',1);
            Zout2(i,:) = temp;
            clear temp;
        end
        temp = fscanf(fid,'%*s%g%*s%g%*s%g',6);
        result(1) = temp(2);
        
        %% fine
        for i = 1:26
            temp = fscanf(fid,'%*s%g',1);
            Zout3(i,:) = temp;
            clear temp;
        end
        temp = fscanf(fid,'%*s%g%*s%g%*s%g',6);
        
        for i = 1:21
            temp = fscanf(fid,'%*s%g',1);
            Zout4(i,:) = temp;
            clear temp;
        end
        temp = fscanf(fid,'%*s%g%*s%g%*s%g',6);
        result([2 4]) = 128;
        result(3) = temp(2);
        result(5) = temp(1);
        clear temp;
        
    otherwise % above 10kHz
        file_path = strcat(folder_path,'\OutputImpedance\',freq,'\Log\Zout\Zout',int2str(ch),'Ch.txt');
        fid = fopen(file_path);
        temp=fscanf(fid,'%s',52);
        
        %% course RC
%         temp=fscanf(fid,'%s',50);
        for i = 1:26
            temp = fscanf(fid,'%*s%*s%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g',26);
            Zout1(i,:) = temp;
            clear temp;
        end
        temp = fscanf(fid,'%*s%g%*s%g%*s%g',3);
        
        temp=fscanf(fid,'%s',42);
        for i = 1:21
            temp = fscanf(fid,'%*s%*s%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g',21);
            Zout2(i,:) = temp;
            clear temp;
        end
        temp = fscanf(fid,'%*s%g%*s%g%*s%g',3);
        result(1:2) = temp(2:3);
        clear temp;
        
        %% fine RC
        temp=fscanf(fid,'%s',52);
        for i = 1:26
            temp = fscanf(fid,'%*s%*s%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g',26);
            Zout3(i,:) = temp;
            clear temp;
        end
        temp = fscanf(fid,'%*s%g%*s%g%*s%g',3);
        
        temp=fscanf(fid,'%s',42);
        for i = 1:21
            temp = fscanf(fid,'%*s%*s%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g',21);
            Zout4(i,:) = temp;
            clear temp;
        end
        temp = fscanf(fid,'%*s%g%*s%g%*s%g',3);
        result(3:4) = temp(2:3);
        result(5) = temp(1);
        clear temp;
end

% subplot(221); surfl(Zout1); title('Course (10 step)');
% subplot(222); surfl(Zout2); title('Course (1 step)');
% subplot(223); surfl(Zout3); title('Fine (10 step)');
% subplot(224); surfl(Zout4); title('Fine (1 step)');