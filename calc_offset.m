function offset = calc_offset(kx_arr,ky_arr,kz_arr)
% Calculate for each non-cartesian point, the offset 
% (relative to the first location of the cartesian array 
% (which corresponds to (kx,ky,kz) = (-64,-64,-64)).
% offset := the index described above.
% offset = calc_offset(kx,ky,kz);

% w := width of the interpolation kernel.
w = 3;
offset = int32(zeros(length(kx_arr),343));
parfor i = 1:length(kx_arr)
   kx = kx_arr(i);
   ky = ky_arr(i);
   kz = kz_arr(i);

   % Find a local nbd around each non-cartesian k-space location.
   % cx := index of the cartesian array. Since (-64,-64,-64) is the first
   % location of this array, we need to do -(-64) = +64 indicated below.
   cx_min = floor(kx) - w + 64;
   cx_max = floor(kx) + w + 64;
   cx_box = [cx_min:cx_max]';

   cy_min = floor(ky) - w + 64;
   cy_max = floor(ky) + w + 64;
   cy_box = [cy_min:cy_max]';

   cz_min = floor(kz) - w + 64;
   cz_max = floor(kz) + w + 64;
   cz_box = [cz_min:cz_max]';

   % Increment in x location is an increment of 1 linear index.
   % Increment in y location is an increment of 128 (1 column) linear index.
   % Increment in z location is an increment of 128^2 (1 matrix) linear index.

   % cx_box = cx_box*1;
   cy_box = cy_box*128;
   cz_box = cz_box*128*128;

   % x changes the fastest followed by y and z.
   % Equivalent statement in C : offset_xyz = (cz << 14) + (cy << 7) + cx
   offset_xy = bsxfun(@plus,cx_box,cy_box');
   offset_xy = offset_xy(:);
   offset_xyz = bsxfun(@plus,offset_xy,cz_box');
   offset(i,:) = offset_xyz(:)';
end
offset = int32(offset);
offset = offset';
end
% fid = fopen('X:\aroor_temp\offset_1000','w');
% fwrite(fid,offset,'int32');
% fclose(fid);
%------ C code to be emulated -------%
% for (int i = 0; i < n; i++)
% {
%    cx_min = floor(*kx) - w + 64;
%    cx_max = floor(*kx) + w + 64;
%    cy_min = floor(*ky) - w + 64;
%    cy_max = floor(*ky) + w + 64;
%    cz_min = floor(*kz) - w + 64;
%    cz_max = floor(*kz) + w + 64;
% 
%    // Loop over the local nbd.
%    for (cz = cz_min; cz <= cz_max; cz++)
%    {
%       for (cy = cy_min; cy <= cy_max; cy++)
%       {
%          for (cx = cx_min; cx <= cx_max; cx++)
%          {
%             out_r = out_origin_r + (cz << 14) + (cy << 7) + cx;
%             out_i = out_origin_i + (cz << 14) + (cy << 7) + cx;
%          }
%       }
%    }
%    // Since we move to a new non-cartesian location, increment the data and
%    // location pointer.
%    kx++;
%    ky++;
%    kz++;
% }
