function y = ifftn_unit(x)
% Unitary FFT operator
s = size(x);
n = sqrt(prod(s));
y = n*ifftn(x);
end