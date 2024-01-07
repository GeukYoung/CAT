function [real, quad] = fxTxtImport(data_path)
temp = load(data_path);
if size(temp,1)>256
    temp(257:end,:) = [];
end
real = temp(:,3);
quad = temp(:,4);
mask = [1,2,16,17,18,19,34,35,36,51,52,53,68,69,70,85,86,87,102,103,104,119,120,121,136,137, ...
138,153,154,155,170,171,172,187,188,189,204,205,206,221,222,223,238,239,240,241,255,256;];
real(mask,:) = [];
quad(mask,:) = [];
end