function [temp2] = FxByte_conv(Data,nByte)
    if rem(size(Data,1),nByte) ~= 0
        Data = Data';
    end
    nData = size(Data,1)/nByte;
    
    for cntData = 1:nData
        temp = zeros(1,size(Data,2));
        for cntByte = 1:nByte
            temp = temp + Data((cntData-1)*nByte+cntByte,:)*256^(cntByte-1);
        end
        temp2(cntData,:) = temp;
    end
