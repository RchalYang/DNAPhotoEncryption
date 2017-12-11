function result=my_bin2dec(bin)
result = 0;
for i = 1:64
    result = result + str2num(bin(i+2))*2^(-i);
end