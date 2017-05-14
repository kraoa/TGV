function B = bitreverse(N)

N_bits = ceil(log2(N));

B_tmp = zeros(1, N);

for j = 0 : N-1
    tmp = j;
    for i = 0 : N_bits
        B_tmp(j+1) = B_tmp(j+1) + 2^(N_bits - i - 1) * mod(tmp, 2);
        tmp = floor(tmp/2);
    end
end
%B_tmp
B(1) = B_tmp(1);
curr_high = 0;

for j = 2 : N
    tmp = N+1;
    diff_elt = 1;
    for i = 2 : N
        if B_tmp(i) > curr_high
            diff = B_tmp(i) - curr_high;
            if diff < tmp
               diff_elt = i;
               tmp = diff;
            end
        end
    end
        curr_high = B_tmp(diff_elt);
        B(diff_elt) = j-1;
   % end
end
