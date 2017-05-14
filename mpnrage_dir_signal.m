function [s1 s2] = mpnrage_dir_signal(TR,T1,alpha,beta,TI,TD,num_ti,W,varargin)
num_args = 8;

if(TR<0)
    error('TR must be greater than 0');
end

% if(T1<0)
%     error('T1 must be greater than 0');
% end

if(numel(alpha)~=2)
    error('must provide two flip angles for alpha');
end
% if( min(alpha) < 0 )
%     error('alpha cannot be negative');
% end

if(numel(beta)~=2)
    error('must provide two flip angles for beta');
end
% if( min(beta)<0 )
%     error('beta cannot be negative');
% end

if(numel(TI)~=2)
    error('must provide two values for TI');
end

if(numel(TD)~=2)
    error('must provide two values for TD');
end

if(numel(num_ti)~=2)
    error('must provide two values for num_ti');
end

if(numel(W)~=2)
    error('must provide two values for W');
end

N = num_ti(1);
M = num_ti(2);

E1 = exp(-TR/T1);

EI = exp(-TI/T1);

ED = exp(-TD/T1);

a = cosd(alpha).*E1;
Mss = (1-E1)./(1-a);

s = cosd(alpha).*cosd(beta).*ED.*EI;

AN = a(1)^(num_ti(1)-1);
AM = a(2)^(num_ti(2)-1);
den = 1-s(1)*s(2)*AN*AM;

if( nargin == 10 )
    m = varargin{2};
else
    m=0:(M-1);
end

am = a(2).^m;

if(M>0)
    
w = ones( num_ti(2), 1 )*W(2);
temp=(numel(w)-0)-(0:1:numel(w)-1);
idx = find( temp < W(2) );
w(idx) = temp(idx);
w = w';
if( nargin== 10 )
    w = w(varargin{2}+1);
end


amvs = (1-(a(2).^w)).*(a(2).^(m))./(1-a(2))./w;
    
term(:,1) = Mss(2).*(1-amvs);

term(:,2) = Mss(2)*(1-AM)*s(2)*s(1)*AN*amvs/den;

term(:,3) = (1-ED(2))*EI(2)*cosd(beta(2))*s(1)*AN*amvs/den;

term(:,4) = (1-EI(2))*amvs/den;

term(:,5) = Mss(1)*(1-AN)*ED(1)*cosd(alpha(1))*EI(2)*cosd(beta(2))*amvs/den;

term(:,6) = (1-ED(1))*EI(2)*cosd(beta(2))*amvs/den;

term(:,7) = AN*(1-EI(1))*ED(1)*cosd(alpha(1))*EI(2)*cosd(beta(2))*amvs/den;

s2 = sum(term,2);
else
    s2 = [];
end

if( nargin == 10 )
    n = varargin{1};
else
    n=0:(N-1);
end

an = a(1).^n;

w = ones( num_ti(1), 1 )*W(1);
temp=(numel(w)-0)-(0:1:numel(w)-1);
idx = find( temp < W(1) );
w(idx) = temp(idx);
w = w';
if( nargin== 10 )
    w = w(varargin{1}+1);
end

anvs = (1-(a(1).^w)).*(a(1).^(n))./(1-a(1))./w;

clear term
term(:,1) = Mss(1).*(1-anvs);

term(:,2) = Mss(1)*(1-AN)*s(1)*s(2)*AM*anvs/den;

term(:,3) = (1-ED(1))*EI(1)*cosd(beta(1))*s(2)*AM*anvs/den;

term(:,4) = (1-EI(1))*anvs/den;

term(:,5) = Mss(2)*(1-AM)*ED(2)*cosd(alpha(2))*EI(1)*cosd(beta(1))*anvs/den;

term(:,6) = (1-ED(2))*EI(1)*cosd(beta(1))*anvs/den;

term(:,7) = AM*(1-EI(2))*ED(2)*cosd(alpha(2))*EI(1)*cosd(beta(1))*anvs/den;

s1 = sum(term,2);


























