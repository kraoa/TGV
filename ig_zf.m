function im_zf = ig_zf(x,sz,coil_sens,ker,offset,D)
% Zero filled inverse gridding.
% sz := 3d x coil x time.
% D := kr.
im_zf = zeros(sz(1),sz(2),sz(3),sz(5));
% Inverse grid with density compensation.
for i_tr = 1:sz(5)
   for i_coil = 1:sz(4)
      % data on the cartesian grid.
      tmp_d = D{i_tr}.^2;
      tmp_d = tmp_d(:);
      tmp_x = x{i_tr};
      tmp_x = tmp_x(:,i_coil);
      tmp = n2c_crop_fast1(ker{i_tr},offset{i_tr},tmp_x(:).*tmp_d(:));
      tmp = fftshift(reshape(tmp,sz(1),sz(2),sz(3)));
      im_zf(:,:,:,i_tr) = im_zf(:,:,:,i_tr) + ifftn_unit(tmp).*conj(coil_sens(:,:,:,i_coil,1));
   end
end
end