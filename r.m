% 第一种方法用到的自相关函数
function res = r( x, N, M)
 res = zeros(N-M+1,1);
for k=0:N-M
    for n = 1:M
        res(k+1) = res(k + 1) + x(n)*conj(x(n+k));
    end
    res(k+1) = res(k + 1) / M;
end