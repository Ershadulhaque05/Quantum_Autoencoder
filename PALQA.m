clc
clear all;
Im = imread('P_gray.jpg');
%Im=Im(1:224,1:224);
[row1 col1]=size(Im);

Img2=Im(1:row1,1:col1);

blocksize=8;

PSNR1=0;
absB1q=0;
X_max_p=0;
BR=0;
Z=zeros(row1, col1);
Y=zeros(row1, col1);
q=8; %Q=8,16,32,70,120

for i=1:blocksize:row1
  for j=1:blocksize:col1

    W1=Img2(i:i+blocksize-1,j:j+blocksize-1);
       
    DCT=dct2(W1);
    
    Q_bit=dec2bin(q,log2(q)+1);
    Q_bit=Q_bit-48;
    Q_bit=flip(Q_bit);
    Q_bit_t=transpose(Q_bit);
    [Q_x Q_Y]= find(Q_bit_t);
    
    Q_bit_po=nnz(dec2bin(Q_x)-48)+1;

    B1q=round(DCT/q);
    
    B1q_sign=sign(B1q);

    
    [D_sig E_Sign F]=find(B1q_sign);

    absB1q= abs(B1q);
    
    [x y z]=find(absB1q);


    
    G=F.*z;
           
    dectoBin=dec2bin(z,8);
    
    dectoBina=uint16(dectoBin);
    dectoBina=dectoBina-48;
    numone=nnz(dectoBina);
    
    rowb=dec2bin(x,log2(blocksize)+1);
    
    X_max=max(x); %check maximum value in the position

    if X_max==8
        X_max_p=X_max_p+X_max;
    else
    end

    urowb=uint16(rowb);
    aurowb1=urowb-48;
    aurowb=aurowb1(:,:);

    aurowb_1=aurowb(:,4);
    N_row= aurowb(:,1:3);
     
    [X_p,Y_P]=size(aurowb_1);

   [C_x C_y C_z]=find(aurowb_1);
    
   %% For x position recovery, LSB based 
    
    
    if sum(aurowb_1)==0
        A=zeros(X_p,Y_P);
        
    else
        A=zeros(X_p,Y_P);
        %[C_x C_y C_z]=find(aurowb_1);
        
        for l=1:1:length(C_z)
            A(C_x(l),C_y(l))=C_z(l);
        end
    end 

    N_row(:,4)=A;

    N_num_str=num2str(N_row);
    N_row_bin_dec=bin2dec(N_num_str);


    erowb=nnz(aurowb(:,1:3));
    
    
    colb=dec2bin(y,log2(blocksize)+1);

    ucolb=uint16(colb);
    aucolb=ucolb-48;

    N_col=num2str(aucolb);

    N_col_bin_dec=bin2dec(N_col);

    H=zeros(8,8);

    for l=1:1:length(G)
        H(N_row_bin_dec(l),N_col_bin_dec(l))=G(l);
    end
    
    aucolb=aucolb(:,:);


    ecolb=nnz(aucolb);

    %col_b_max=max(y);
    %if col_b_max==8
    %    C_max=col_b_max;
   % else
   % end
   


    statebit=erowb+nnz(A)+ecolb;
    
    sbit=nnz(B1q_sign);
    
    auxbit=numone;
    
    
    %BR=BR+(numone+sbit+auxbit)/(1024*1024);
    BR=BR+(Q_bit_po+numone+sbit+statebit+auxbit)/(1024*1024);

    B2=H.*q;
    
   %B2=B1q.*q;

    RI1=idct2(B2);
    Z(i:i+7,j:j+7)=RI1;

    Y(i:i+7,j:j+7)=absB1q;
  end
end



PSNR1=CalculatePSNR(Z,Img2);


numofblockr=row1/blocksize;
numofblockc=col1/blocksize;

rposioflastblock=row1-blocksize;
cposioflastblock=col1-blocksize;

numofcolbit=nnz(dec2bin(rposioflastblock))+nnz(dec2bin(blocksize));
numofrowbit=nnz(dec2bin(cposioflastblock))+nnz(dec2bin(blocksize));
tbitr=numofblockr*numofblockc*(numofcolbit+numofrowbit)/(1024*1024);


bpp=(BR+tbitr)*(1024*1024)/(row1*col1)

PSNR1