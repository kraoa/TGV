function DATA = gen_sos_ir_1(k,num_t,coil_sens)
% generate stack of stars vfa data based on julia's dataset t1 and m0 maps.
% this accepts k-space data as input.
% addpath(genpath('./NUFFT'));

% num_t = 300;
TR = 7e-3;
TD = [500 500]*1e-3+TR;
TI = [12 12]*1e-3;
beta = [155 155];
fa = [4 4];
% The temporal sampling period is tr.
t = [0:num_t-1]'*TR;
num_t = numel(t);
N = num_t;
M = num_t;
num_coil = size(coil_sens,3);
% [t1_fs,m0_fs] = calc_t1_sos_8c('fs');
load('/home/aroor/matlab/code/junta/karthikaroor/fmri/matlab/mri/t1_fs_steve')
load('/home/aroor/matlab/code/junta/karthikaroor/fmri/matlab/mri/rho_fs_steve')
load('/home/aroor/matlab/code/junta/karthikaroor/fmri/matlab/mri/mask_steve');
t1_fs(165,93) = t1_fs(165,92);
rho_fs(165,93) = rho_fs(165,92);
mask(165,93) = mask(165,92);
t1_fs = t1_fs.*1e-3;
t1_fs = t1_fs.*mask;

% m0_fs = ones(size(t1_fs)).*mask;
m0_fs = rho_fs.*mask;
warning('change kappa later');
kappa = ones(size(t1_fs));
idx = (t1_fs ~= 0) & (m0_fs ~= 0);
% s := source image.
s = m0_fs;

sz_im = size(t1_fs);

s_arr = zeros([sz_im num_t]);

for i = 1:sz_im(1)
    for j = 1:sz_im(2)
        if mask(i,j) == 1
            [f,g] = mpnrage_2fa(TR,t1_fs(i,j),fa,beta,TI,TD,[N,M],kappa(i,j));
            s_arr(i,j,:) = f;
        end
    end
end

% multiply by rho
for i_t = 1:num_t
    s_arr(:,:,i_t) = s_arr(:,:,i_t).*m0_fs;
end

% add noise
% compute the std dev of noise
tmp_t1 = 1000e-3;
tmp_E1 = exp(-TR/tmp_t1);
tmp_a1 = cosd(1.*fa(1))*tmp_E1;
tmp_m_inf = median(m0_fs(m0_fs>0)).*(1-tmp_E1)./(1-tmp_a1).*sind(1.*fa(1));
snr = 30;
sig = tmp_m_inf/snr;

% s_arr = s_arr + sig*randn(size(s_arr));

% s_arr1 = calc_ma(s_arr,10);

% [k,w] = gen_kspace;
% [k,w] = gen_kspace(num_t);
w = abs(k);


% 
% load('/scratch/aroor_sos/coil_sens.mat');
% for i_coil = 1:num_coil
%     coil_sens(:,:,i_coil) = coil_sens(:,:,i_coil).*mask;
% end

warning('tbd add noise to kspace samples. actually this is done but i need to double check');
tmp_FT = NUFFT(k(:,:,1), 1, 1, 0, [240,240], 2);
tmp = tmp_FT*(s_arr(:,:,1).*coil_sens(:,:,1));
X = zeros([size(tmp) num_coil num_t]);
for i_t = 1:num_t
    s = s_arr(:,:,i_t);
    FT = NUFFT(k(:,:,i_t), 1, 1, 0, [240,240], 2);
    for i_coil = 1:num_coil
        x = s.*coil_sens(:,:,i_coil);
        X(:,:,i_coil,i_t) = FT*x;
%         noise = sig*randn(size(s));
% %         X(:,:,i_coil,i_t) = FT*(x+noise);
%         noise = sig*randn(size(X(:,:,i_coil,i_t)));
%         X(:,:,i_coil,i_t) = X(:,:,i_coil,i_t) + noise;
%         x1 = FT'*(w.*X);
    end
end
% add noise in k-space
rng(1);
X = X + sig*randn(size(X));

DATA.X = X;
DATA.k = k;
DATA.w = w;
% save('/scratch/aroor_sos/DATA_ir_300t','DATA','-v7.3');

% % r = temp_gen_sos_ir;
% r = temp_gen_sos_ir(num_t);
% 
% % DEL
% for i_t = 1:num_t
% FT = NUFFT(k(:,:,i_t), 1, 1, 0, [256,256], 2);
% y(:,:,i_t) = IFT(FT,w(:,:,i_t),X(:,:,:,i_t),coil_sens,[256 256 num_coil num_t]);
% end
% % DEL
% % adjust for the scaling.
% for i_t = 1:num_t
% y(:,:,i_t) = r*y(:,:,i_t);
% end
% % DEL
% i_t = 300;
% imagesc(abs(y(:,:,i_t))-abs(s_arr(:,:,i_t)),[-1e-6 1e-6])
end

function y = IFT(K,D,x,coil_sens,sz)
% x := #non-cartesian x #coils.
% y := 3 spatial.
% sz := 3d x #coil x #parametric.

num_coil = sz(3);
% y := 3 spatial x #coils.
y = zeros([sz(1:3)]);
for i_coil = 1:num_coil
   tmp = x(:,:,i_coil);
   im_n = D.*tmp;
   % Gridding from non-cart to cart.
   im_c = K'*(im_n);
   y(:,:,i_coil) = conj(coil_sens(:,:,i_coil)).*(im_c);
end
y = sum(y,3);
end

function s = calc_vfa(fa,t1,m0)
% Calculate the vfa curve at flip angles fa.
% [s,si,sl] = calc_vfa(fa);
% tr = 3.8e-3;
tr = 6e-3;
t1 = t1*1e-3;
e = exp(-tr./t1);
s = m0.*(1-e)*sind(fa)./(1-e*cosd(fa));
end

function [k,w] = gen_kspace(num_t)
% k - num_rad x num_proj x num_t
% eg. k - 526 x 402 x num_t
dir_in = '/scratch/aroor_sos/VOL080811/';
load([dir_in 'kspace_traj']);
% choose 1 projection to be used as a baseline
% we rotate all other projections
kb = [kxx(:,1) kyy(:,1)]'*0.5/128;

% the fully sampled data has 402 projections for image size of 256
num_proj = 402;
R = @(th) [cosd(th) -sind(th);sind(th) cosd(th)];

k = zeros(size(kxx,1),num_proj,num_t);
for i_t = 1:num_t
    for i_p = 1:num_proj
        th = unifrnd(0,180);
        tmp = R(th)*kb;
%         plotv(tmp*128/0.5);
%         xlim([-128 128]);
%         ylim([-128 128]);
        tmp = tmp';
        k(:,i_p,i_t) = tmp(:,1) + 1i*tmp(:,2);
    end
end

w = abs(k);
end

function [k,w] = gen_kspace_old
% numSpokes = 402;
% numSamplesOnSpoke = 526;
% angles = (1:numSpokes)/numSpokes*pi;
% line = linspace(-0.5,0.5,numSamplesOnSpoke);
% [arg1,arg2] = meshgrid(angles,line);
% coordx = arg2.*cos(arg1);
% coordy = arg2.*sin(arg1);
% k = coordx + 1i*coordy;
% 
% w = (abs(line)*numSamplesOnSpoke/numSpokes/2)';
% w = repmat(w, [1, numSpokes]);

dir_in = '/scratch/aroor_sos/VOL080811/';
load([dir_in 'kspace_traj']);
kx = kxx*0.5/128;
ky = kyy*0.5/128;

k = kx + 1i*ky;
w = abs(k);
% k = k(:);
% w = w(:);

end

function r = temp_gen_sos_ir(num_t)
% r := ratio of the two norms.
% [t1_fs,m0_fs] = calc_t1_sos_8c('fs');
load('/scratch/aroor_sos/t1_fs')
load('/scratch/aroor_sos/m0_fs')
s = m0_fs;
[k,w] = gen_kspace(num_t);
FT = NUFFT(k, 1, 1, 0, [256,256], 2);
load('/scratch/aroor_sos/coil_sens.mat')
for i_coil = 1:num_coil
    x = s.*coil_sens(:,:,i_coil);
    X = FT*x;
    x1(:,:,i_coil) = (FT'*(w.*X)).*conj(coil_sens(:,:,i_coil));
end
sh = sum(x1,3);
% determine the scaling.
r = norm(s(:))/norm(sh(:));
end
