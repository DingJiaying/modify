% 第一种方法对序列的CLS_BC
function res = BC( x, N )
res = 0;
for n = 3:N
    res = res +x(n-1).*conj(x(n-1));
end
end