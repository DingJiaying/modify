% 第一种方法信号序列的CRPHD_BR
function res = BR( x, N )
res = norm(x(N)).^2 - norm(x(N-1)).^2 - norm(x(2)).^2 + norm(x(1)).^2;
for n = 3:N
    res = res + 2*conj(x(n))*x(n-2);
end
end