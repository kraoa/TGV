function ker = calc_ker(kx_arr,ky_arr,kz_arr)
% Calculate for each non-cartesian point, the offset 
% (relative to the first location of the cartesian array 
% (which corresponds to (kx,ky,kz) = (-64,-64,-64)).
% offset := the index described above.
% ker = calc_ker(kx,ky,kz);

% w := width of the interpolation kernel.
w = 3;
ker = single(zeros(length(kx_arr),343));
parfor i = 1:length(kx_arr)
   kx = kx_arr(i);
   ky = ky_arr(i);
   kz = kz_arr(i);

   % Find a local nbd around each non-cartesian k-space location.
   % cx := index of the cartesian array. Since (-64,-64,-64) is the first
   % location of this array, we need to do -(-64) = +64 indicated below.
   cx_min = floor(kx) - w;
   cx_max = floor(kx) + w;
   cx_box = [cx_min:cx_max]';

   cy_min = floor(ky) - w;
   cy_max = floor(ky) + w;
   cy_box = [cy_min:cy_max]';

   cz_min = floor(kz) - w;
   cz_max = floor(kz) + w;
   cz_box = [cz_min:cz_max]';

   dx = kx - (cx_box);
   dy = ky - (cy_box);
   dz = kz - (cz_box);

   ker_x = wsinc(dx);
   ker_y = wsinc(dy);
   ker_z = wsinc(dz);
   % x changes the fastest followed by y and z.
   ker_xyz = kron(ker_z',ker_x*ker_y.');
   ker(i,:) = ker_xyz(:)';
end

% ker = single(ker);

% fid = fopen('X:\aroor_temp\ker_1000','w');
% fwrite(fid,ker,'float');
% fclose(fid);

ker = ker';
end

function y = wsinc(x)
% Windowed sinc with a Gaussian window.
w = exp(-(x/3).^2);
y = sinc(x).*w.*(abs(x)<3);
end

function y = kb(x)
% Kaiser Bessel window.
n = 128;
g = 1*n;
a = g/n;
w = 6;
b = pi*sqrt(w^2/a^2*(a-0.5)^2 - 0.8);
c = g/w*besseli(0,b*sqrt(1-(2*a*x/w).^2));
i = abs(x) > w/(2);
c(i) = 0;
y = c;
end

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
