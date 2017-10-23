% 第一种方法信号序列的CRPHD_AR
function res = AR( x, N )
res = 0;
for n = 3:N
   res = res + conj(x(n-1))*(x(n-2)+x(n)); 
end
end                                                                       