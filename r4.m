%�����ַ����õ�������غ���
function res = r4( x, N, k )
res = 0;
for n = k+1:N
    res = res + x(n)*x(n-k);
end
res = res / (N-k);
end