function str=gainlevelwatermark(B)  
wr=size(B,1);
wl=wr;
str='';
    for m=1:wr
        for n=1:wl
           waterbin=dec2bin(B(m,n),8);
           str=strcat(str,waterbin);
        end
    end
end