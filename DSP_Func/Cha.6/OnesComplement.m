function [y] = OnesComplement(x,B)
% y = OnesComplement(x,B)
% ---------------
% Decimal equivalent of Sign-Magnitude format integer to b-bit Ones¡¯-Complement format conversion
% 
% x:integer between -2^(b-1) < x < 2^(b-1) (sign-magnitude)
% y:integer between 0 <= y <= 2^b-1 (1's-complement)
%
if any((x <= -2^(B-1) | (x >= 2^(B-1))))
    error('Number must satisfy -2^(B-1) < x < 2^(B-1)')
end
s = sign(x); % sign fo x (-1 if x<0, 0 if x=0, 1 if x>0)
sb = (s < 0); % sign-bit (0 if x>=0, 1 if x<0)
y = (1-sb).*x + sb.*(2^B+x-1);
end

