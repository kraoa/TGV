function [kx,ky,kr,img_k_arr,dat_nc] = read_kspace_rad_vfa_vd_coil(vd,coil,num_proj,DATA)
% Read radial kspace locations and data (multiple coils and flip angles).
% num_proj := number of projections to be used.

kx = real(DATA.k(:,:,vd));
ky = imag(DATA.k(:,:,vd));
kr = DATA.w(:,:,vd);

idx = [1:num_proj];

kx = kx(:,idx);
ky = ky(:,idx);
kr = kr(:,idx);
img_k_arr = DATA.X(:,idx,coil,vd);
dat_nc = img_k_arr.*kr;

% Vectorize the output.
kx = double(kx(:));
ky = double(ky(:));
kr = double(kr(:));
img_k_arr = double(img_k_arr(:));
dat_nc = double(dat_nc(:));

end