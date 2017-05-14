function u = tgv_vfa(K, coil_sens, f, kr, alpha0, alpha1, maxits, reduction, sz, num_proj,idx_param)
% Function implements tgv in the parametric dimension.
% load /scratch/aroor_sos/u_fs
check_it = 1;

alpha00 = alpha0;
alpha10 = alpha1;
alpha01 = alpha0*reduction;
alpha11 = alpha1*reduction;

% sz = [128 128 128 32 num_t]; % 3d x #coil x #parametric
num_coil = sz(3);
num_tr = sz(4);
sz_st = [sz(1:2) sz(4)];
for i_tr = 1:num_tr
    num_nc(i_tr) = size(f{i_tr},1);
end

% DEL
% for i_tr = 1:10
% kr{i_tr} = kr{i_tr}/2;
% end

for i_tr = 1:num_tr
    for i_coil = 1:num_coil
%         f{i_tr}(:,i_coil) = f{i_tr}(:,i_coil) + 2e-5*randn(size(f{i_tr}(:,i_coil)));
        f{i_tr}(:,i_coil) = f{i_tr}(:,i_coil).*sqrt(kr{i_tr});
    end
end

% v = zeros(size(f));
for i_tr = 1:num_tr
   v{i_tr} = 0;
end
p = zeros(sz_st);
q = zeros(sz_st);

u_zf = zeros(sz_st);
for i_tr = 1:num_tr
    for i_coil = 1:num_coil
%         f_dc{i_tr}(:,i_coil) = f{i_tr}(:,i_coil).*(kr{i_tr});
        f_dc{i_tr}(:,i_coil) = f{i_tr}(:,i_coil);
    end
    u_zf(:,:,i_tr) = IFT(K{i_tr},kr{i_tr},f_dc{i_tr},coil_sens,sz);
end
save(['/scratch/aroor_tgv_ivptm/u_zf_' num2str(idx_param)],'u_zf');

for i_tr = 1:num_tr
    tmp_zf{i_tr} = FT(K{i_tr},kr{i_tr},u_zf(:,:,i_tr),coil_sens,sz,num_nc(i_tr));
    tmp_scl(i_tr) = (tmp_zf{i_tr}(:)'*f{i_tr}(:))/(tmp_zf{i_tr}(:)'*tmp_zf{i_tr}(:));
%     u_zf(:,:,i_tr) = u_zf(:,:,i_tr)*tmp;
%     tmp_zf{i_tr} = FT(K{i_tr},kr{i_tr},u_zf(:,:,i_tr),coil_sens,sz,num_nc(i_tr));
end

tmp_scl = mean(tmp_scl);

for i_tr = 1:num_tr
    u_zf(:,:,i_tr) = u_zf(:,:,i_tr)*tmp_scl;
end

% u_zf = u_zf*pi/num_proj;
u = u_zf;
xi = zeros(sz_st);
u_ = u;
xi_ = xi;

% e = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multiindices of the spatial derivatives

% derivatives
% | uxx uxy |
% | uyx uyy |

% multiindices
% | 1 3 |
% | 3 2 |

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% L = sqrt(64);
L = 1;
tau_p = 1/4;
tau_d = 1/4;

for k=1:1000
%     k
%     if alpha00 ~= 0
%       alpha0 = exp(k/maxits*log(alpha01) + (maxits-k)/maxits*log(alpha00)); alpha0 = max(alpha0,2e-8);
%       alpha1 = exp(k/maxits*log(alpha11) + (maxits-k)/maxits*log(alpha10)); alpha1 = max(alpha1,1e-8);
%       alpha0_arr(k+1) = alpha0;
%       alpha1_arr(k+1) = alpha1;
%     end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SAVE VARIABLES
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  uold = u; xiold = xi;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % DUAL UPDATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  clear tmp_r
  % operator
  parfor i_tr = 1:num_tr
      Ku0 = FT(K{i_tr},kr{i_tr},u_(:,:,i_tr),coil_sens,sz,num_nc(i_tr));
      Ku_ = Ku0;
      r =  Ku_ - f{i_tr};
      tmp_r{i_tr} = r;
      v{i_tr} = (v{i_tr} + tau_d*r)/(1+tau_d*L);
   end
  
  % gradient
  du__p = dp(u_);
  
  p = p + tau_d*(du__p - xi_);
  
  % projection
  absp = sqrt(abs(p).^2);
  denom = max(1,absp/alpha1);
  p = p./denom;
  
  % symmetrized gradient
  dxi__p = dp(xi_);
  
  q = q + tau_d*dxi__p;
  
  % projection
  absq = sqrt(abs(q).^2);
  denom = max(1,absq/alpha0);
  q = q./denom;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PRIMAL UPDATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  % dual operator
  ww = zeros(sz_st);
  parfor i_tr = 1:num_tr
     ww(:,:,i_tr) = IFT(K{i_tr},kr{i_tr},v{i_tr},coil_sens,sz);
  end
  
  % divergence
  divp = -idp(p);
  
  u = u + tau_p*(divp - ww);
  
  % divergence
  divq = -idp(q);
  
  xi = xi + tau_p*(p + divq);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % AUXILIARY UPDATE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  u_ = 2*u - uold;
  xi_ = 2*xi - xiold;
  
%   if mod(k,check_it) == 0
%     primal_dual_energy = 0;
%     dual_energy = sum(sum(abs(u).^2))/2.0;
%     e = [e dual_energy];
%     fprintf('TGV2-L2-2D-PD: it = %4d, alpha0 = %f, rmse = %f\n', k, alpha0, 0);    
%     imshow(abs(u),[]);
%     drawnow;
%   end
  tmp_r = cell2mat(tmp_r);
  e(k) = norm(tmp_r(:));
  display(num2str(k));
  display(num2str(e(k)));
save(['/scratch/aroor_tgv_ivptm/u_' num2str(idx_param)],'u');
end
end

function y = dp(x)
y = x(:,:,[2:end,end]) - x;
end

function y = idp(x) 
y = x(:,:,[1,1:end-1]) - x;
y(:,:,1) = -x(:,:,1);
y(:,:,end) = x(:,:,end-1);
% y = cat(4,x(:,:,:,1:end-1),zeros(:,:,:)) - cat(4,zeros(:,:,:),x(:,:,:,1:end-1));
end

function y = FT(K,D,x,coil_sens,sz,num_nc)
num_coil = sz(3);
y = zeros(num_nc,num_coil);
for i_coil = 1:num_coil
   % im_c := the k-space data for the input b and a correspoding coil
   % sensitivity for different i_coil.
   y(:,i_coil) = NUFFT_mtimes(K,(coil_sens(:,:,i_coil).*x));
end
for i_coil = 1:num_coil
    y(:,i_coil) = y(:,i_coil).*sqrt(D);
end
end

function y = IFT(K,D,x,coil_sens,sz)
% x := #non-cartesian x #coils.
% y := 3 spatial.
% sz := 3d x #coil x #parametric.

num_coil = sz(3);
% y := 3 spatial x #coils.
y = zeros([sz(1:3)]);
for i_coil = 1:num_coil
   im_n = sqrt(D).*x(:,i_coil);
   % Gridding from non-cart to cart.
   Kt = NUFFT_ctranspose(K);
   im_c = NUFFT_mtimes(Kt,(im_n));
   y(:,:,i_coil) = conj(coil_sens(:,:,i_coil)).*(im_c);
end
y = sum(y,3);
end

function y = FT1(x,coil_sens,ker,offset,sz,num_nc)
num_coil = sz(4);
y = zeros(num_nc,num_coil);
for i_coil = 1:num_coil
   % im_c := the k-space data for the input b and a correspoding coil
   % sensitivity for different i_coil.
   im_c = fftn_unit(coil_sens(:,:,:,i_coil).*x);
   % The center of k-space needs to be at the center of the 3d array for
   % the c-gridding code to work.
   im_c = ifftshift(im_c);
   % Gridding from cart to non-cart.
   y(:,i_coil) =  c2n_crop_fast1(ker,offset,im_c);
end
y = y/20;
end

function y = IFT1(x,coil_sens,ker,offset,sz)
% x := #non-cartesian x #coils.
% y := 3 spatial.
% sz := 3d x #coil x #parametric.

num_coil = sz(4);
% y := 3 spatial x #coils.
y = zeros([sz(1:4)]);
for i_coil = 1:num_coil
   im_n = x(:,i_coil);
   % Gridding from non-cart to cart.
   im_c = n2c_crop_fast1(ker,offset,im_n);
   im_c = reshape(im_c,sz(1:3));
   % Get the center of k-space to the top left corner of the array.
   im_c = fftshift(im_c);
   y(:,:,:,i_coil) = conj(coil_sens(:,:,:,i_coil)).*ifftn_unit(im_c);
end
y = sum(y,4);
y = y/20;
end
