function result=my_dec2bin(dec,pre)
result = '0.';
for i = 1:pre
    dec = dec*2;
    if dec>1
        result = strcat(result,'1');
        dec = dec-1;
    else
        result = strcat(result,'0');
    end
end