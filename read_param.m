function [TV_TSWeight, num_proj, sz_bin] = read_param(idx_param)
dat = dlmread('param_cs_3d_vfa_rad.txt');
% st := start of each block of parameters.
st = find(dat == 666666);

idx = find(dat(st + 1) == idx_param);
% skip the 666666 and idx_param.
idx = st(idx)+1;

TV_TSWeight = dat(idx+1);

num_proj = dat(idx+2);

sz_bin = dat(idx+3);
end