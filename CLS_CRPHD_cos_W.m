%%三种自相关序列的随着SNR的变化
clc;
close all;
clear all;

%% 信号参数
N = 200;

L = 50;  %CLS的参数
M = 100;
l1=1;
l2=floor(N/3); %Close-form expended粗估计

p=2;   
q = 20; %Multiple autocorrelation lags

wv = 0 : pi/20 : pi;   %从pi/25到pi-pi/25每隔pi/25取一个
A = sqrt(2);
phiA = 0*pi;


%% 独立运行
SNR= 30;
T = 1000;

index = 1;
for w0 = wv
    wAcu = w0 * ones(1,T);             %真实的信号
    for t = 1:T     
        n = 1:1:N;
        sn = A*cos(w0*n + phiA);
        xn = awgn(sn,SNR);
        rn = r(xn, N,M);
        RN = length(rn);

        %% CLS-CRPHD方法的参数
         v_current = zeros(L,1); %x(n+1)
         v_prev_1 = zeros(L,1);  %x(n)
         v_prev_2 = zeros(L,1);  %x(n-1)
         Fre_CLS = zeros(RN,1);
         w = zeros(RN,1);
         for kk = RN:-1:(L+2)  
             if kk > (L + 1)
                %% CLS方法
                v_current = [rn(kk:-1:(kk-L+1))];   %x(n+1)
                v_prev_1 = [rn((kk-1):-1:(kk-L))];  %x(n)
                v_prev_2 = [rn((kk-2):-1:(kk-L-1))]; %x(n-1)
                w(kk) = (v_prev_1'*(v_current+v_prev_2))/( (v_prev_1)'*(v_prev_1)); 

                %% CRPHD方法
                CRPHD_A = real(AR(rn(kk:-1:(kk-L+1)), L));
                CRPHD_B = real(BR(rn(kk:-1:(kk-L+1)), L));
                argument_r = (CRPHD_B+sqrt(CRPHD_B^2+8*CRPHD_A^2))/(4*CRPHD_A);  
                if (argument_r>1)
                    argument_r = 1;
                end
                if (argument_r<-1)
                   argument_r = -1;
                end
             end
             Fre_CLS(kk) = acos(0.5 * w(kk));
             Fre_CRPHD(kk) = acos(argument_r);         
         end
         CLS_gu(t) = mean(Fre_CLS(L+2:RN));
         CRPHD_gu(t) = mean(Fre_CRPHD(L+2:RN));
      
         %% Close-form expended
            for ii=1:N-1                  %%产生自相关序列
                Rn3(ii)=r3(xn,N,ii);
            end
            rn3 = [Rn3(l1:l2)];   %%截取自相关序列以减弱噪声的影响
            RN3 = length(rn3);
            A3=real(AC( rn3, RN3));
            B3=real(BC(rn3, RN3));
            argumentc3 =A3/(2*B3);
            if (argumentc3<-1)
                display('-1')
                argumentc3 = -1;
            end
            if (argumentc3>1)
                display('1')
                argumentc3 = 1;
             end
             CEA_gu(t) = acos(argumentc3);
         %% Multiple autocorrelation lags
             for ii=1:N-1                  %%产生自相关序列
                Rn4(ii)=r4(xn,N,ii);
             end
                rn4 = [Rn4(p:q)];   %%截取自相关序列以减弱噪声的影响
                RN4 = length(rn4);

                A4=real(AC( rn4, RN4));
                B4=real(BC(rn4, RN4));
                argumentc4 =A4/(2*B4);
                if (argumentc4<-1)
                    display('-1')
                    argumentc4 = -1;
                end
                if (argumentc4>1)
                    display('1')
                    argumentc4 = 1;
                end
            MAL_gu(t) = acos(argumentc4);
    end    
    MSE_CLS(index)  = 10*log10(mse(wAcu, CLS_gu));
    MSE_CRPHD(index)  = 10*log10(mse(wAcu, CRPHD_gu));  
    MSE_CEA(index) = 10*log10(mse(wAcu, CEA_gu));
    MSE_MAL(index) = 10*log10(mse(wAcu, MAL_gu));
    index = index + 1;
    index
end
%CRLB参考下界
CRLBtemp = 0;
w1 = 0*pi:0.01*pi:0.99*pi;
for k = 0:N-1
%     CRLBtemp = CRLBtemp + (k^2*(sin(w1*k+phi1).^2));
%     CRLBtemp = CRLBtemp + (k*2*pi*sin(w1*k+phi1).^2);
    CRLBtemp = CRLBtemp + (k*k/8*pi*sin(w1*k+phiA).^2);
end

CRLB = 1 ./ (10^(SNR/10)*CRLBtemp);

plot(wv/pi,  MSE_CLS, '*', 'LineWidth', 2);
hold on;
plot(wv/pi,   MSE_CEA, '^', 'LineWidth', 2);
plot(wv/pi,   MSE_MAL, 'd', 'LineWidth', 2);
plot(wv/pi,   MSE_CRPHD, 'o', 'LineWidth', 2);

plot(w1/pi, 10*log10(CRLB),'k-', 'LineWidth', 2);
%仿真格式部分
% set(gca, 'xtick', -10:5:50);
legend('CLS','Close-form expended','Multiple autocorrelation lags','CRPHD','CRLB');%用指定的文字在当前坐标轴中对所给数据的每一部分显示一个图例
xlabel('\fontname{Times New Roman}SNR', 'FontWeight','bold');%字体Times New Roman，加粗
ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
grid on;
hold off;
