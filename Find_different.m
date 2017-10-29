%%整理之后的CLS
clc;
close all;
clear all;

%% 信号参数
N = 200;  %序列长度
n = 1:1:N; 
M  = 100;      %自相关序列的长度


w = 0.3*pi;
V1 = sqrt(2);
phi1 = 0*pi;
sn = V1*cos(w*n + phi1);


%% 独立运行
SNR = 20;
T = 100;     %%独立运行的次数
wAcu = w * ones(1,T);             %真实的信号
for t = 1:T
    xn = awgn(sn,SNR);
    rn1 = r1(xn, N,M);
    RN1 = length(rn1);  %%自相关序列的长度
    
     %% CLS方法
     v_current = [rn1(3:1:RN1)];   %x(n+1)
     v_prev_1 = [rn1(2:1:RN1-1)];  %x(n)
     v_prev_2 = [rn1(1:1:RN1-2)]; %x(n-1)

     argumentc1 = (v_prev_1'*(v_current+v_prev_2))/(2* (v_prev_1)'*(v_prev_1)); 
     if (argumentc1<-1)
        display('-1')
        argumentc1 = -1;
    end
    if (argumentc1>1)
        display('1')
        argumentc1 = 1;
    end
     CLS_gu1(t) = acos(argumentc1);
      %% 第二种CRPHD的方法
     A=real(AR(rn1,RN1));
    B=real(BR(rn1,RN1));
    argumentc = (B+sqrt(B^2+8*A^2))/(4*A); 
    if (argumentc<-1)
        display('-1')
        argumentc = -1;
    end
    if (argumentc>1)
        display('1')
        argumentc = 1;
     end
     CRPHD_gu(t) = acos(argumentc);
     
end
MSE_CLS = 10*log10(mse(wAcu, CLS_gu1))
MSE_SCRPHD = 10*log10(mse(wAcu, CRPHD_gu))

 