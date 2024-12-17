function levelimage=gainlevelimage(watermarkstr)

     B=watermarkstr;
     lB=size(B,2);
      for i=1:lB/8
          extract(1,i)=bin2dec(B(1,(i-1)*8+1:(i-1)*8+8));
      end
      extracts=reshape(extract,32,32);
      levelimage=extracts';
end