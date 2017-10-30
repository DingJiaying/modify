%%测试M对估计性能的影响
%%整理之后的CLS针对自相关序列
clc;
close all;
clear all;

%% 信号参数
N = 200;
n = 0:1:N-1;

w = 0.3*pi;
V1 =1;
V2 = 2;
phi1 = 0.1*pi;
phi2 = 0.2*pi;
sn = V1*exp(1j*w*n - 1j*phi1) +  V2*exp(-1j*w*n -1j * phi2);

SNR =20;

%% 独立运行

T = 100;
M = 100;


    wAcu = w * ones(1,T);             %真实的信号
    for m = 1:T

        xn = awgn(sn, SNR); 
        rn = r(xn, N,M);
         RN = length(rn);
         %% CLS的关键参数
        
        A_CLS=real(AC(rn,RN));
        B_CLS=real(BC(rn,RN));
        argumentc_CLS = A_CLS/(2*B_CLS); 
        if (argumentc_CLS<-1)
            display('-1')
            argumentc_CLS = -1;
        end
        if (argumentc_CLS>1)
            display('1')
            argumentc_CLS = 1;
         end
         CLS_gu(m) = acos(argumentc_CLS);
        
    end

    MSE_CLS = 10*log10(mse(wAcu, CLS_gu))
% plot(M, MSE_CLS3, '-*', 'LineWidth', 2);
% %仿真格式部分
% legend('SNR=20dB');
% % legend('SNR=-10dB','SNR=0dB','SNR=20dB','SNR=40dB','SNR=60dB');%用指定的文字在当前坐标轴中对所给数据的每一部分显示一个图例
% xlabel('\fontname{Times New Roman}M', 'FontWeight','bold');%字体Times New Roman，加粗
% ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
% grid on;
% hold off;