%%����M�Թ������ܵ�Ӱ��
%%����֮���CRPHD������������
clc;
close all;
clear all;

%% �źŲ���
N = 200;
n = 0:1:N-1;

w = 0.3*pi;
V1 =1;
V2 = 2;
phi1 = 0.1*pi;
phi2 = 0.2*pi;
sn = V1*exp(1j*w*n - 1j*phi1) +  V2*exp(-1j*w*n -1j * phi2);

%% ��������
SNR =20;
T = 100;
M = 100;

index = 1;

    wAcu = w * ones(1,T);             %��ʵ���ź�
    for m = 1:T

        xn = awgn(sn, SNR); 
        rn = r(xn, N,M);
          
         %% CRPHD�Ĺؼ�����
        RN = length(rn);
        A_CRPHD=real(AR(rn,RN));
        B_CRPHD=real(BR(rn,RN));
        argumentc_CRPHD = (B_CRPHD+sqrt(B_CRPHD^2+8*A_CRPHD^2))/(4*A_CRPHD); 
        if (argumentc_CRPHD<-1)
            display('-1')
            argumentc_CRPHD = -1;
        end
        if (argumentc_CRPHD>1)
            display('1')
            argumentc_CRPHD = 1;
         end
         CRPHD_gu(m) = acos(argumentc_CRPHD);
     
    end
    MSE_CRPHD= 10*log10(mse(wAcu, CRPHD_gu))


% plot(M, MSE_CLS3, '-*', 'LineWidth', 2);
% 
% %�����ʽ����
% legend('SNR=20dB');
% % legend('SNR=-10dB','SNR=0dB','SNR=20dB','SNR=40dB','SNR=60dB');%��ָ���������ڵ�ǰ�������ж��������ݵ�ÿһ������ʾһ��ͼ��
% xlabel('\fontname{Times New Roman}M', 'FontWeight','bold');%����Times New Roman���Ӵ�
% ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
% grid on;
% hold off;