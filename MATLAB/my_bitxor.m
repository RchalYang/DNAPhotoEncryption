function result=my_bitxor(bin1,bin2)
result = '0.';
for i = 1:64
    result = strcat(result,num2str(bitxor(str2num(bin1(i+2)),str2num(bin2(i+2)))));
end