function [A, kr_arr, coil_sens, im_k_arr] = read_kspace_vfa_2(fa,num_proj,sz_bin)
% Here we load k-space locations and data for each coil and tr separately. 
% This is slower but the code is simpler since I've tested 
% read_kspace_rad_vfa_vd_coil.
dir_in_cs = '/study/mrphys/karthik/in/dat_t1_rho_b1_steve/';
% n_r := number of samples along each radial line i.e. each spoke.
n_r = 526;
% n_th := number of spokes.
% For fully sampled the #spokes is pi/2*128^2 = 25736.
warning('fix this #fs projs');
n_th = 332;
% num_coil := number of coils.
num_coil = 1;
% num_tr := number of flip angles.
num_tr = length(fa);

load('/home/aroor/matlab/code/junta/karthikaroor/fmri/matlab/mri/mask_steve');

warning('change coil sense later');
% % load coil sensitivity
% load([dir_in_cs 'coil_sens'])
% for i_coil = 1:num_coil
%     coil_sens(:,:,i_coil) = coil_sens(:,:,i_coil).*mask;
% end
coil_sens = ones(240,240,num_coil);
% i_t = 2;
% generate kspace locations
% plot(imag(k1(:,1:2,i_t)),real(k1(:,1:2,i_t)),'LineWidth',2)
k1 = pr_order5_test(num_proj,num_tr);
% % use same kspace samples at all time points
% for i = 1:300
%     k1(:,:,i) = k1(:,:,1);
% end
% generate kspace data
DATA = gen_sos_ir_1(k1,num_tr,coil_sens);

for i_tr = 1:num_tr-(sz_bin-1)
   i_tr;
   k{i_tr} = [];
   kr_arr{i_tr} = [];
   im_k_arr{i_tr} = [];
   
   for j = 0:(sz_bin-1)
       im_k = [];
       for i_coil = 1:num_coil
          [kx,ky,kr,im_k(:,i_coil),~] = read_kspace_rad_vfa_vd_coil(fa(i_tr+j),i_coil,num_proj,DATA);
       end
       k{i_tr} = [k{i_tr};kx + 1i*ky];
       kr_arr{i_tr} = [kr_arr{i_tr};kr];
       im_k_arr{i_tr} = [im_k_arr{i_tr};im_k];
   end
      
   A{i_tr} = NUFFT_main(k{i_tr}, 1, 1, 0, [240,240], 2);
   
end
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
