function y = neville(xx,n,x,Q)
% Neville's algorithm as a function (save as "nev.m")
% 
% inputs:
%    n = order of interpolation (n+1 = # of points)
%    x(1),...,x(n+1)    x coords
%    Q(1),...,Q(n+1)    y coords
%    xx=evaluation point for interpolating polynomial p
%
% output:  p(xx)
if isempty(Q)
    y = [];
else
for i = n:-1:1
   for j = 1:i
      Q(j) = (xx-x(j))*Q(j+1) - (xx-x(j+n+1-i))*Q(j);
      Q(j) = Q(j)/(x(j+n+1-i)-x(j));
   end
end
y = Q(1);
end
end