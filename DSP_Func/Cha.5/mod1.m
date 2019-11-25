function m = mod1(n,N)
% Computes m = (n mod N) index, 
% note that the arguments in MATLAB begin in 1, not 0. 
% ----------------------------
% m = mod(n,N)
m = rem(n,N); m = m+N; m = rem(m,N);
end

