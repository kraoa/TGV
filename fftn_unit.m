function y = fftn_unit(x)
% Unitary FFT operator
s = size(x);
n = sqrt(prod(s));
y = (1/n)*fftn(x);
end