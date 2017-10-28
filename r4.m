%第五种方法用到的自相关函数
function res = r4( x, N, k )
res = 0;
for n = k+1:N
    res = res + x(n)*x(n-k);
end
res = res / (N-k);
end