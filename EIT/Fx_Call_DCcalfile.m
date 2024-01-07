function [result] = Fx_Call_DCcalfile(folder_path,freq_num)
% Function info
% [result] = Fx_Call_CCScalfile(forler_path,freq,ch)
% input files -> 2 nargin
%       [folder path, freq, IMM ch]
%        folder path (str) -> ex) D:\Dropbox\#Lab Work\1. EIT_System\1. 16ch EIT\[FINAL_SYDNEY]EIT_Mark25_20120128\Debug\Calibration\eit1
%        freq (int) -> 1 : 1.125kHz
%                      2 : 11.25KHz
%                      3 : 56.25kHz
%                      4 : 112.5KHz
% output files -> DC offset cal result

freq_index = {'1.125kHz','12.5KHz','62.5kHz','125kHz'};
freq = freq_index{freq_num};

file_path = strcat(folder_path,'\DCOffset\',freq,'\DCOffset.txt');
result = load(file_path);
