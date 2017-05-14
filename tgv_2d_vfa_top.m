function tgv_2d_vfa_top(idx_param)
% this is for steve's dataset in which is the result averaged from 8 scans.

% here we generate the kspace locations and data on the fly.
% we solve the optimization problem,
% || Ax - y ||^2 + TGV(x)

% the time series follows the mpnrage signal model.

addpath(genpath('./NUFFT'));
addpath(genpath('./GRIDDING'));
idx_param = str2num(idx_param);

% num_proj := number of projections (fully sampled data has pi/2*256 projections).
% sz_bin := size of the bin
[TV_TSWeight, num_proj, sz_bin] = read_param(idx_param)
t = [1:386]';
% num_t := number of inversion time points.
num_t = length(t);

reduction = 2^(-8);     % usually there is no need to change this
alpha = TV_TSWeight
maxits = 500;
num_coil = 1;
[A, kr, coil_sens, im_k_arr] = read_kspace_vfa_2(t,num_proj,sz_bin);
% n(1:2) := spatial dimension.
% n(3) := number of coils.
% n(4) := number of time points.
sz = [240 240 num_coil num_t-(sz_bin-1)];%num_tr/4

% TGV iterations
tic
im_cs = tgv_vfa(A, coil_sens, im_k_arr, kr, 2*alpha, alpha, maxits, reduction, sz, num_proj,idx_param);
toc

end
