function coil_sens = calc_cs(img_arr,num_vd,num_coil)
% Calculate coil sensitivities.
s = size(img_arr);

den = sqrt(sum(abs(img_arr).^2,4));
mask_brain = 1;
% mask_brain = abs(den) > 5;

for vd=1:num_vd
    for coil=1:num_coil
%         mask_brain(:,:,:,vd);
        coil_sens(:,:,:,coil,vd) = img_arr(:,:,:,coil,vd).*mask_brain./den(:,:,:,1,vd);
    end
end

% if( sum( isinf(coil_sens(:) ) ) )
%     error('infinite coil sensitivity values detected\n');
% end
% 
% if( sum( isnan(coil_sens(:) ) ) )
%     error('isnan coil sensitivity values detected\n');
% end
% 
% coil_sens(isinf(coil_sens)) = 0;
% coil_sens(isnan(coil_sens)) = 0;

end


% den = zeros([s(1:3),s(5)]);
% for vd=1:num_vd
%     for coil=1:num_coil
%         den(:,:,:,vd) = den(:,:,:,vd) + abs(img_arr(:,:,:,coil,vd)).^2;
%     end
%     den(:,:,:,vd) = sqrt(den(:,:,:,vd));
%     mask_brain(:,:,:,vd) = 1;
% %     mask_brain(:,:,:,vd) = abs(den(:,:,:,vd)) > 5;
% end