%�����ַ����õ�������غ���
function res = r3( x, N, k )
res = 0;
for n = 1:N-k
    res = res + x(n)*x(n+k);
end
res = res / (N-k);
end