function g = eval_tukeywin(k,a)
% I have only tested this for k = [-64:63]';
% a that best eliminates the gibbs ringing at 2mm is 
% a = 0.9
N = 128;%length(k);
n = k + N/2;
g = zeros(size(k));
idx = (0 <= n) & (n <= a*(N-1)/2);
g(idx) =  (1/2)*(1 + cos(pi*( (2*n(idx))/(a*(N-1)) - 1 )));

idx = ((N-1)*(1-a/2) <= n) & (n <= N-1);
g(idx) =  (1/2)*(1 + cos(pi*( (2*n(idx))/(a*(N-1)) -2/a + 1 )));

idx = (a*(N-1)/2 <= n) & (n <= (N-1)*(1-a/2));
g(idx) =  1;

% % h should exactly be equal to g if k is a regular +/- grid 
% % (eg. k = [-64:63]';).
% h = tukeywin(N,a);
% if norm(h-g)/N > 1e-14
%     error('something is wrong in tukey window generation');
% end
end