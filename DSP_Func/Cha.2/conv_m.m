function [y,ny] = conv_m(x,nx,h,nh)
% Modified convolution routine for signal processing
% --------------------------------------------------
% [y,ny] = conv_m(x,nx,h,nh)
% [y,ny] = convolution result
% [x,nx] = first signal
% [h,nh] = second signal
% when 'n' in the convolution formula equal to 0, the index of the corresponding value is nyb.  
nyb = nx(1)+nh(1); nye = nx(length(x)) + nh(length(h));
ny = (nyb:nye); y = conv(x,h);
end

