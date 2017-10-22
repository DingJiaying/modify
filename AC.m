% 第一种方法对序列的CLS_AC
function res = AC( x, N )
res = 0;
for n = 3:N
   res = res + conj(x(n-1))*(x(n-2)+x(n)); 
end
end