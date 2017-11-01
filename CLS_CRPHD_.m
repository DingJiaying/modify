%%整理之后的对三种自相关序列的CLS
clc;
close all;
clear all;

%% 信号参数
N = 100;
n = 0:1:N-1;

M = N/2;

w = 0.3*pi;
V1 = sqrt(2);
phi1 = 0*pi;
sn = V1*cos(w*n + phi1);

%% 独立运行
SNR = 20;
T = 100;
wAcu1 = w * ones(1,T);             %真实的信号
for m = 1:T
    xn = awgn(sn,SNR);
    
    rn1 = r1(xn, N,M);RN1 = length(rn1);
    Rx_k=rn1(2:N-M);
    Rx_k_1=rn1(1:N-M-1);
    Rx_k1=rn1(3:N-M+1);  
    %% 第一种CLS的方法
    A_CLS=real(Rx_k'*(Rx_k_1+Rx_k1));
    B_CLS=2*norm(Rx_k)^2;
    argumentc1 =A_CLS/(B_CLS);
    if (argumentc1<-1)
        display('-1')
        argumentc1 = -1;
    end
    if (argumentc1>1)
        display('1')
        argumentc1 = 1;
     end
     CLS_gu1(m) = acos(argumentc1);
    %% 第二种CRPHD的方法
    A_r=Rx_k'*(Rx_k_1+Rx_k1);
    B_r=norm(Rx_k_1)^2+norm(Rx_k1)^2-2*norm(Rx_k)^2+2*real(Rx_k_1'*Rx_k1);
    argumentc = (B_r+sqrt(B_r^2+8*A_r^2))/(4*A_r);
    if (argumentc<-1)
        display('-1')
        argumentc = -1;
    end
    if (argumentc>1)
        display('1')
        argumentc = 1;
     end
     CRPHD_gu(m) = acos(argumentc); 

    
end
MSE_CLS1 = 10*log10(mse(wAcu1, CLS_gu1))
MSE_CRPHD = 10*log10(mse(wAcu1, CRPHD_gu))
