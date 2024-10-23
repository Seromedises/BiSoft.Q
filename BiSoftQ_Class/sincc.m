function y = sincc(x)
% normalized sinc function, sin(pi*x)/(pi*x), no checks on the input
%
% y = sincc(x)
y = sin(pi*x)./(pi*x);
y(x==0) = 1;
end