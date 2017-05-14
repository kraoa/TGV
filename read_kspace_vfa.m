function [A, kr_arr, coil_sens, im_k_arr] = read_kspace_vfa(fa,num_proj,DATA)
% Here we load k-space locations and data for each coil and tr separately. 
% This is slower but the code is simpler since I've tested 
% read_kspace_rad_vfa_vd_coil.
dir_in_cs = '/scratch/aroor_sos/';
% n_r := number of samples along each radial line i.e. each spoke.
n_r = 526;
% n_th := number of spokes.
% For fully sampled the #spokes is pi/2*128^2 = 25736.
n_th = 332;
% num_coil := number of coils.
num_coil = 8;
% num_tr := number of flip angles.
num_tr = length(fa);

load('/home/aroor/matlab/code/junta/karthikaroor/fmri/matlab/mri/mask_sos_8c');

fid = fopen([dir_in_cs 'coil_sens_r'],'r');
coil_sens_r = fread(fid,'double');
fclose(fid);
fid = fopen([dir_in_cs 'coil_sens_i'],'r');
coil_sens_i = fread(fid,'double');
fclose(fid);
coil_sens = coil_sens_r + 1i*coil_sens_i;
coil_sens = reshape(coil_sens,[256 256 20 8]);
coil_sens = coil_sens(:,:,11,:);
coil_sens = squeeze(coil_sens);
% % we can instead load coil as follows
% load /scratch/aroor_sos/coil_sens

for i_coil = 1:num_coil
    coil_sens(:,:,i_coil) = coil_sens(:,:,i_coil).*mask;
end

for i_tr = 1:num_tr
   im_k = [];
   for i_coil = 1:num_coil
      [kx,ky,kr,im_k(:,i_coil),~] = read_kspace_rad_vfa_vd_coil(fa(i_tr),i_coil,num_proj,DATA);
   end
%    kx_arr{i_tr} = kx;
%    ky_arr{i_tr} = ky;
%    kz_arr{i_tr} = kz;
   k = kx + 1i*ky;
   A{i_tr} = NUFFT(k, 1, 1, 0, [256,256], 2);
   kr_arr{i_tr} = kr;
   im_k_arr{i_tr} = im_k;
%    A{i_tr} = make_grid_opt_2d([256 256],kx,ky);
%    A{i_tr} = make_grid_opt_2d_test([256 256],kx,ky);
%    im_k_arr{i_tr} = lpf_tukeywin(im_k,kx,ky,kz);
%    offset_arr{i_tr} = calc_offset(kx,ky,kz);
%    ker_arr{i_tr} = calc_ker(kx,ky,kz);
end
% DEL
% make the signal in image-space to be a linear function across time.
% for i_tr = 1:num_tr
%     im_k_arr{i_tr} = i_tr*im_k_arr{1};
% end
% load '/scratch/aroor_sos/A_66'
% load '/scratch/aroor_sos/A_332'
% load '/scratch/aroor_sos/A_kb'
end

function im_k_lpf = lpf_tukeywin(im_k,kx,ky,kz)
% Here we apply a low pass filter to the input k-space data.

% This is only valid for image size of 128.
num_coil = 32;
k = [-64:63]';
a = 0.9;
gx = eval_tukeywin(kx,a);
gy = eval_tukeywin(ky,a);
gz = eval_tukeywin(kz,a);
for i_coil = 1:num_coil
im_k_lpf(:,i_coil) = im_k(:,i_coil).*gx.*gy.*gz;
end
end


function im_k_lpf = lpf(im_k,kx,ky,kz)
% Here we apply a low pass filter to the input k-space data.

% This is only valid for image size of 128.
num_coil = 32;
kc = 1;
k = [-64:63]';
w = max(abs(k)) - kc;
% g = cos(pi*(abs(k)-kc)/(2*w));
% ng = norm(g);
gx = cos(pi*(abs(kx)-kc)/(2*w));
gy = cos(pi*(abs(ky)-kc)/(2*w));
gz = cos(pi*(abs(kz)-kc)/(2*w));
for i_coil = 1:num_coil
im_k_lpf(:,i_coil) = im_k(:,i_coil).*gx.*gy.*gz;
end
end