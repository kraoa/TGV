% clc
% clear
% close all

% addpath('/home/aroor/matlab/code/junta/karthikaroor/fmri/matlab/mri/GRIDDING');
function k = pr_order5_test(subproj,interleaved_angles)

% % subproj := number of projections at each time point
% subproj = 300;

% % interleaved_angles := number of time points
% interleaved_angles = 300;

interleaves = interleaved_angles;

num_views = interleaves*subproj;

num_slices=1;
num_views=interleaved_angles*subproj;
num_partitions=subproj;

interleave_order_host = bitreverse(interleaved_angles);
sub_order = bitreverse(subproj);

for pos_slice=0:(num_slices-1) 
  for pos=0:(num_views-1)
     inter_index=mod(pos,interleaved_angles);
     inter_pos=floor( pos / interleaved_angles );
     sub_pos = inter_index + inter_pos;

     aov(pos + pos_slice*num_views + 1) = interleaved_angles*sub_order(mod(sub_pos,num_partitions)+1) ...
         + interleave_order_host(inter_index+1);
  end
end

AOV=reshape(aov,[interleaved_angles subproj])';

theta=linspace(0,pi,numel(AOV)+1);
theta(end)=[];
idx = find(theta>pi);
theta(idx) = theta(idx) + 0.5*(theta(2)-theta(1));


angles_1 = theta(AOV+1);

opxres=256;
for mf=1:size(angles_1,1)
    [KX(:,mf,:) KY(:,mf,:)] = get_traj(opxres,angles_1(mf,:));
end

k = KX + 1i*KY;
k = k*0.5/128;
return;


ColOrd = get(gca,'ColorOrder');
figure
w=2
for cf=1:size(angles_1,2)
    st=cf;
    en=st+w-1;
    if en > size(KX,3)
        en = size(KX,3)
    end
    hold off
    for pos=st:en
        if en < size(KX,3)
            kx = reshape( KX(:,:,pos), size(KX,1), [] );
            ky = reshape( KY(:,:,pos), size(KX,1), [] );
            plot( kx, ky,'Color',ColOrd(mod(pos,size(ColOrd,1))+1,:) );axis([-opxres/2 opxres/2 -opxres/2 opxres/2]);pause(0.5)
            hold on
        end
    end
end

return;

xx = cosd(theta(:,:,1));
yy = sind(theta(:,:,1));

XX = cosd(theta(:,:,1)+180);
YY = sind(theta(:,:,1)+180);


echo = 1:num_echoes;
figure, 
for tr=1:interleaved_trs
    plot( xx(echo,tr), yy(echo,tr),'bx','LineWidth',10 );daspect([1 1 1]);
    axis([-1 1 -1 1]);
    hold on
    plot( XX(echo,tr), YY(echo,tr),'rx','LineWidth',10 );daspect([1 1 1]);    
    pause
end

return
st=1;
en=1;

st_ileave=1;
en_ileave=st_ileave+interleaved_angles-1;

t1=AOV(st:en,st_ileave:en_ileave);
t2=AOV2(st:en,st_ileave:en_ileave);


sort(t1(:))'
sort(t2(:))'

S(:,1) = sort(t1(:));
S(:,2) = sort(t2(:));


D=diff(S);
D(end+1,:)=0;

[S D]

mean_vals = mean(D)

std_vals = std(D,0,1)
end