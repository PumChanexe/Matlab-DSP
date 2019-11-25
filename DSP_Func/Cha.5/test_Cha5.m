clc; clear; close;
%% dfs function test
% xn = [0,1,2,3]; N = 4; Xk = dfs(xn,N);

%% Example 5.2
% L = 5; N = 20; k = (-N/2:N/2); % Sq wave parameters
% xn = [ones(1,L), zeros(1,N-L)]; % Sq wave x(n)
% Xk = dfs(xn,N); % DFS
% magXk = abs([Xk(N/2+1:N) Xk(1:N/2+1)]); % DFS magnitude
% stem(k,magXk); axis([-N/2,N/2,-0.5,5.5])
% xlabel('k'); ylabel('Xtilde(k)')
% title('DFS of SQ. wave: L=5, N=20')

%% Example 5.5 discuss sample number N's effect
% N = 5; k = 0:1:N-1; 
% wk = 2*pi*k/N; zk = exp(1i*wk); 
% Xk = (zk)./(zk - 0.7);
% xn = real(idfs(Xk,N)); 
% xtilde = xn'*ones(1,8); 
% xtilde = (xtilde(:))'; % 5*8 --> 40*1 --> 1*40
% subplot(2,2,1); stem(0:39,xtilde); axis([0,40,-0.1,1.5]) 
% xlabel('n'); ylabel('xtilde(n)'); title('N=5'); 

%% Example 5.6
% x = [1,1,1,1,0,0,0,0]; N = 8; X = dft(x,N); 
% magX = abs(X); phaX = angle(X)*180/pi; 

%% Example 5.8
% n = (0:1:99); x = cos(0.48*pi*n)+cos(0.52*pi*n);
% --------------------------------------------------------------------------
% n1 = (0:1:9); y1 = x(1:1:10);
% subplot(2,1,1); stem(n1,y1); title('signal x(n), 0 <= n <= 9'); xlabel('n')
% Y1 = dft(y1,10); magY1 = abs(Y1(1:1:6));
% k1 = 0:1:5; w1 = 2*pi/10*k1;
% subplot(2,1,2);stem(w1/pi,magY1);title('Samples of DTFT Magnitude');
% xlabel('frequency in pi units')
% ------------------------------------------------------------------------
% n2 = (0:1:99); y2 = [x(1:1:10) zeros(1,90)];
% subplot(2,1,1) ;stem(n2,y2) ;title('signal x(n), 0 <= n <= 9 + 90 zeros');
% xlabel('n');
% Y2 = dft(y2,100); magY2 = abs(Y2(1:1:51));
% k2 = 0:1:50; w2 = 2*pi/100*k2;
% subplot(2,1,2); plot(w2/pi,magY2); title('DTFT Magnitude');
% xlabel('frequency in pi units');
% ------------------------------------------------------------------------
% subplot(2,1,1); stem(n,x);
% title('signal x(n), 0 <= n <= 99'); xlabel('n')
% X = dft(x,100); magX = abs(X(1:1:51));
% k = 0:1:50; w = 2*pi/100*k;
% subplot(2,1,2); plot(w/pi,magX); title('DTFT Magnitude');
% xlabel('frequency in pi units')

%% Example 5.9
% n = 0:10; x = 10*(0.8).^n; y = x(mod(-n,11)+1); % note the argument.
% subplot(2,1,1); stem(n,x); title('Original sequence')
% xlabel('n'); ylabel('x(n)');
% subplot(2,1,2); stem(n,y); title('Circularly folded sequence')
% xlabel('n'); ylabel('x(-n mod 10)');
% 
% X = dft(x,11); Y = dft(y,11);
% subplot(2,2,1); stem(n,real(X));
% title('Real{DFT[x(n)]}'); xlabel('k');
% subplot(2,2,2); stem(n,imag(X));
% title('Imag{DFT[x(n)]}'); xlabel('k');
% subplot(2,2,3); stem(n,real(Y));
% title('Real{DFT[x(-n)11]}'); xlabel('k');
% subplot(2,2,4); stem(n,imag(Y));
% title('Imag{DFT[x(-n)11]}'); xlabel('k');

%% Example 5.10
% n = 0:10; x = 10*(0.8).^n;
% [xec,xoc] = circevod(x);
% subplot(2,1,1); stem(n,xec); title('Circular-even component')
% xlabel('n'); ylabel('xec(n)'); axis([-0.5,10.5,-1,11])
% subplot(2,1,2); stem(n,xoc); title('Circular-odd component')
% xlabel('n'); ylabel('xoc(n)'); axis([-0.5,10.5,-4,4])

%% Example 5.12
% n = 0:10; x = 10*(0.8).^n; y = cirshftt(x,6,15);
% n = 0:14; x = [x zeros(1,4)];
% subplot(2,1,1); stem(n,x); title('Original Sequence')
% xlabel('n'); ylabel('x(n)');
% subplot(2,1,2); stem(n,y);
% title('Circularly shifted sequence, N=15')
% xlabel('n'); ylabel('x((n-6) mod 15)')

%% Example 5.14
% x1 = [1,2,2]; x2 = [1,2,3,4]; y = circonvt(x1, x2, 4);

%% Example 5.19
% n = 0:9; x = n+1; h = [1,0,-1]; N = 6; 
% y = ovrlpsav(x,h,N);

%% study the execution time of fft function
% Nmax = 2048; fft_time = zeros(1,Nmax); 
% for n = 1:1:Nmax 
%     x = rand(1,n);
%     t = clock; fft(x); fft_time(n) = etime(clock,t);
% end
% n = (1:1:Nmax); 
% plot(n,fft_time,'.'); xlabel('N'); ylabel('Time in Sec'); title('FFT execution times')

%% Example 5.23 demonstrate the effectiveness of the high-speed convolution
% conv_time = zeros(1,150); fft_time = zeros(1,150);
% for L = 1:150
%     tc = 0; tf=0;
%     N = 2*L-1; 
%     nu = ceil(log10(N)/log10(2)); N = 2^nu;
%     for I=1:100
%         h = randn(1,L); x = rand(1,L);
%         t0 = clock; y1 = conv(h,x); t1=etime(clock,t0); tc = tc+t1;
%         t0 = clock; y2 = ifft(fft(h,N).*fft(x,N)); t2=etime(clock,t0); tf = tf+t2;
%     end
%     conv_time(L)=tc/100; fft_time(L)=tf/100;
% end
% n = 1:150; subplot(1,1,1);
% plot(n(25:150),conv_time(25:150),n(25:150),fft_time(25:150))
% ------------------------ not work --------------------------------




 













