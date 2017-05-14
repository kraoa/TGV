function calc_cs_top(DATA)
% DATA := raw k-space data off the scanner.
dir_in = '/scratch/aroor_sos/VOL080811/';
% load([dir_in 'raw_data']);
load([dir_in 'kspace_traj']);
% Here we make process the data to compute the coil sens.
acc = pi/2*256/(332);
% acc = 100;
% num_proj := number of projections (fully sampled data has pi/2*128^2 projections).
num_proj = round(pi/2*256/332);
% i_tr := the tr used to compute the coil sens.
i_tr = 5;
n = [256 256 8];
% Read the k-space data and gridding variables (but not the coil sens since we haven't computed it yet).
% [ker, offset, kr, ~, im_k_arr] = read_kspace_vfa_single_tr(i_tr,num_proj);


% kx = kxx(:);
% ky = kyy(:);
% idx = (abs(kx)<=124) & (abs(ky)<=124);
% kx = kx(idx)*0.5/128;
% ky = ky(idx)*0.5/128;
kx = kxx*0.5/128;
ky = kyy*0.5/128;
kr = sqrt(kx.^2+ky.^2);
k = kx + 1i*ky;

% A = make_grid_opt_2d([256 256],kx(:),ky(:));
A = NUFFT(k, 1, 1, 0, [256,256], 2);

% x = ifft(DATA,[],3);
x = DATA;
x = x(:,:,:,:,i_tr);
for i_z = 1:size(DATA,3)
x1 = x(:,:,i_z,:);
x1 = squeeze(x1);
im_k_arr = [];
for i_coil = 1:8
%     tmp = x1(:,:,i_coil);
%     tmp = tmp(:);
%     im_k_arr(:,i_coil) = tmp(idx);
    im_k_arr(:,:,i_coil) = x1(:,:,i_coil);
end
im_zf(:,:,i_z,:) = ig(im_k_arr,n,A,kr);
end
coil_sens = calc_cs(im_zf,1,8);

% Store the coil sens.
dir_out = '/scratch/aroor_sos/';
coil_sens_r = real(coil_sens);
coil_sens_i = imag(coil_sens);

fid = fopen([dir_out 'coil_sens_r'],'w');
fwrite(fid,coil_sens_r,'double');
fclose(fid);
fid = fopen([dir_out 'coil_sens_i'],'w');
fwrite(fid,coil_sens_i,'double');
fclose(fid);


% Check if the coil sens is indeed correct.
dir_in = '/scratch/aroor_sos/';
fid = fopen([dir_in 'coil_sens_r'],'r');
coil_sens_r = fread(fid,'double');
fclose(fid);
fid = fopen([dir_in 'coil_sens_i'],'r');
coil_sens_i = fread(fid,'double');
fclose(fid);
coil_sens = coil_sens_r + 1i*coil_sens_i;
coil_sens = reshape(coil_sens,[256 256 20 8]);


for i_coil = 1:8
for i_z = 11
    tmp = (coil_sens(:,:,i_z,i_coil));
    imagesc(abs(tmp(:,:,:)))
%     colorbar
%     impixelinfo
    title(num2str(i_coil))
    pause(1)
end
end



end


function im_zf = ig(x,sz,A,D)
im_zf = zeros(sz(1),sz(2),sz(3));
% Inverse grid with density compensation.
for i_coil = 1:sz(3)
    % data on the cartesian grid.
    tmp_d = D;
%     tmp_d = tmp_d(:);
%     tmp_x = x;
%     tmp_x = tmp_x(:,i_coil);
    tmp_x = x(:,:,i_coil);
    tmp = A'*(tmp_x(:).*tmp_d(:));
%     tmp = (reshape(tmp,sz(1),sz(2)));
%     im_zf(:,:,i_coil) = ifftn_unit(tmp);
    im_zf(:,:,i_coil) = tmp;
end
end
