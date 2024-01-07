clear all;
EITtools_path = 'D:\Dropbox\#Lab Work\2. Recon & simulation\Matlab code\Functions\EITtools.m';
run(EITtools_path);
clc; close all; clear EITtools_path;

Saline_sigma = 0.4;
Ref_path  =   'C:\Users\Jang\Desktop\20161206 MEBIS\20161206 MEBIS\data\1\Conductivity Spectrum';
Folder_path = 'C:\Users\Jang\Desktop\20161206 MEBIS\20161206 MEBIS\data';

for i = 1:7
    Sigma(:,i) = FxMbis_SigmaRecon(Ref_path,strcat(Folder_path,'\',num2str(i),'\Conductivity Spectrum'));
    Sigma(:,i) = Sigma(:,i)
end
Sigma = Sigma*Saline_sigma;
clear all; close all;

%
freqnum_axis = [100, 1000, 5000, 10000, 50000, 100000];
figure;
semilogx(freqnum_axis,Sigma,'--o'); colormap jet; 
xlabel('Frequency(Hz)'); ylabel('Conductivity');
legend('1', '2' , '3' , '4' , '5' , '6' , '7');
% axis([-inf inf 0 0.1]);