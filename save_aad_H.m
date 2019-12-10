H = reshape(aad_H.',[1,512*40]).';
H_real = real(H);
H_imag = imag(H);
H1 = [H_real,H_imag];
H2 = reshape(H1.',[],1);
dlmwrite('aad_H_48K.txt',H2.','delimiter',',','precision','%.6f')
fid1=fopen('hello.txt','w');
count=fprintf(fid1,' %f \n',H_imag);