clc; clear;close;
%% Example 3.3
% w = (0:1:500)*pi/500; % [0, pi] axis divided into 501 points.
% X = exp(1i*w) ./ (exp(1i*w) - 0.5*ones(1,501));
% magX = abs(X); angX = angle(X); realX = real(X); imagX = imag(X);
% 
% subplot(2,2,1); plot(w/pi,magX); grid
% xlabel('frequency in pi units'); title('Magnitude Part'); ylabel('Magnitude')
% subplot(2,2,3); plot(w/pi,angX); grid
% xlabel('frequency in pi units'); title('Angle Part'); ylabel('Radians')
% subplot(2,2,2); plot(w/pi,realX); grid
% xlabel('frequency in pi units'); title('Real Part'); ylabel('Real')
% subplot(2,2,4); plot(w/pi,imagX); grid
% xlabel('frequency in pi units'); title('Imaginary Part'); ylabel('Imaginary')

%% Example 3.4 -- matrix-vector multiplication
% n = -1:3; x = 1:5; k = 0:500; w = (pi/500)*k;
% X = x * (exp(-1i*pi/500)) .^ (n'*k); % (500 & k):value between [0,pi]
% magX = abs(X); angX = angle(X);
% realX = real(X); imagX = imag(X);
% % 
% subplot(2,2,1); plot(k/500,magX);grid
% xlabel('frequency in pi units'); title('Magnitude Part')
% subplot(2,2,3); plot(k/500,angX/pi);grid
% xlabel('frequency in pi units'); title('Angle Part')
% subplot(2,2,2); plot(k/500,realX);grid
% xlabel('frequency in pi units'); title('Real Part')
% subplot(2,2,4); plot(k/500,imagX);grid
% xlabel('frequency in pi units'); title('Imaginary Part')

%% Example 3.5 -- x(n) is complexed valued
% n = 0:10; x = (0.9*exp(1i*pi/3)).^n; 
% k = -200:200; w = (pi/100)*k; % frequency between [-2*pi,2*pi]
% X = x * (exp(-1i*pi/100)) .^ (n'*k); 
% magX = abs(X); angX = angle(X);
% 
% subplot(2,1,1); plot(w/pi,magX);grid
% xlabel('frequency in units of pi'); ylabel('|X|'); title('Magnitude Part')
% subplot(2,1,2); plot(w/pi,angX/pi);grid
% xlabel('frequency in units of pi'); ylabel('radians/pi'); title('Angle Part')

%% Example 3.6
% n = -5:5; x = (-0.9).^n;
% k = -200:200; w = (pi/100)*k; 
% X = x * (exp(-1i*pi/100)) .^ (n'*k);
% magX = abs(X); angX = angle(X);
% 
% subplot(2,1,1); plot(w/pi,magX);grid; axis([-2,2,0,15])
% xlabel('frequency in units of pi'); ylabel('|X|'); title('Magnitude Part')
% subplot(2,1,2); plot(w/pi,angX/pi);grid; axis([-2,2,-1,1])
% xlabel('frequency in units of pi'); ylabel('radians/pi'); title('Angle Part')

%% Example 3.7 verify the linearity property(3.5)
% x1 = rand(1,11); x2 = rand(1,11); n = 0:10;
% alpha = 2; beta = 3; k = 0:500; w = (pi/500)*k; 
% X1 = x1 * (exp(-1i*pi/500)).^(n'*k); % DTFT of x1
% X2 = x2 * (exp(-1i*pi/500)).^(n'*k); % DTFT of x2
% x = alpha*x1 + beta*x2; % linear combination of x1 & x2
% X = x * (exp(-1i*pi/500)).^(n'*k); % DTFT of x
% % verification
% X_check = alpha*X1 + beta*X2; % Linear Combination of X1 & X2
% error = max(abs(X-X_check)); % Difference

%% Example 3.8 verify the shift property(3.6)
% x = rand(1,11); n = 0:10;
% k = 0:500; w = (pi/500)*k;
% X = x * (exp(-1i*pi/500)).^(n'*k); % DTFT of x
% % signal shifted by two samples
% y = x; m = n+2;
% Y = y * (exp(-1i*pi/500)).^(m'*k);
% % verification
% Y_check = (exp(-1i*2).^w).*X; % multiplication by exp(-j2w)
% error = max(abs(Y-Y_check)); % Difference

%% Example 3.9 verify the frequency shift property(3.7)
% n = 0:100; x = cos(pi*n/2);
% k = -100:100; w = (pi/100)*k; % frequency between [-pi,+pi]
% X = x * (exp(-1i*pi/100)).^(n'*k); % DTFT of x
% %
% y = exp(1i*pi*n/4).*x; % DTFT of y
% Y = y * (exp(-1i*pi/100)).^(n'*k); % DTFT of y
% % Graphical verification
% subplot(2,2,1); plot(w/pi,abs(X)); grid; axis([-1,1,0,60])
% xlabel('frequency in pi units'); ylabel('|X|'); title('Magnitude of X')
% subplot(2,2,2); plot(w/pi,angle(X)/pi); grid; axis([-1,1,-1,1])
% xlabel('frequency in pi units'); ylabel('radiands/pi'); title('Angle of X')
% subplot(2,2,3); plot(w/pi,abs(Y)); grid; axis([-1,1,0,60])
% xlabel('frequency in pi units'); ylabel('|Y|'); title('Magnitude of Y')
% subplot(2,2,4); plot(w/pi,angle(Y)/pi); grid; axis([-1,1,-1,1])
% xlabel('frequency in pi units'); ylabel('radians/pi'); title('Angle of Y')

%% Example 3.10 verify the conjugation property (3.8)
% n = -5:10; x = rand(1,length(n)) + 1i*rand(1,length(n));
% k = -100:100; w = (pi/100)*k; % frequency between [-pi,+pi]
% X = x * (exp(-1i*pi/100)).^(n'*k); % DTFT of x
% % conjugation property
% y = conj(x); % signal conjugation
% Y = y * (exp(-1i*pi/100)).^(n'*k); % DTFT of y
% % verification
% Y_check = conj(fliplr(X)); % conj(X(-w))
% error = max(abs(Y-Y_check)); % Difference

%% Example 3.11 verify the folding property (3.9)
% n = -5:10; x = rand(1,length(n));
% k = -100:100; w = (pi/100)*k; % frequency between [-pi,+pi]
% X = x * (exp(-1i*pi/100)).^(n'*k); % DTFT of x
% % folding property
% y = fliplr(x); m = -fliplr(n); % signal folding
% Y = y * (exp(-1i*pi/100)).^(m'*k);
% % verification
% Y_check = fliplr(X); % X(-w)
% error = max(abs(Y-Y_check)); % Difference

%% Example 3.12 verify the symmetry property (3.10)
% n = -5:10; x = sin(pi*n/2);
% k = -100:100; w = (pi/100)*k; % frequency between [-pi,+pi]
% X = x * (exp(-1i*pi/100)).^(n'*k); % DTFT of x
% % signal decomposition
% [xe,xo,m] = evenodd(x,n); % even and odd parts
% XE = xe * (exp(-1i*pi/100)).^(m'*k); % DTFT of xe
% XO = xo * (exp(-1i*pi/100)).^(m'*k); % DTFT of xo
% % verification
% XR = real(X); % real part of X
% XI = imag(X); % imag part of X
% error1 = max(abs(XE-XR)); % Difference
% error2 = max(abs(XO-1i*XI));
% % graphical verification
% subplot(2,2,1); plot(w/pi,XR); grid; axis([-1,1,-2,2])
% xlabel('frequency in pi units'); ylabel('Re(X)'); title('Real part of X');
% subplot(2,2,2); plot(w/pi,XI); grid; axis([-1,1,-10,10])
% xlabel('frequency in pi units'); ylabel('Im(X)'); title('Imaginary part of X');
% subplot(2,2,3); plot(w/pi,real(XE)); grid; axis([-1,1,-2,2])
% xlabel('frequency in pi units'); ylabel('XE'); title('Transform of even part');
% subplot(2,2,4); plot(w/pi,imag(XO)); grid; axis([-1,1,-10,10])
% xlabel('frequency in pi units'); ylabel('XO'); title('Transform of odd part');

%% Example 3.13 frequency response
% w = (0:1:500)*pi/500; % [0, pi] axis divided into 501 points.
% H = exp(1i*w) ./ (exp(1i*w) - 0.9*ones(1,501));
% magH = abs(H); angH = angle(H);
% subplot(2,1,1); plot(w/pi,magH); grid;
% xlabel('frequency in pi units'); ylabel('|H|'); title('Magnitude Response');
% subplot(2,1,2); plot(w/pi,angH/pi); grid
% xlabel('frequency in pi units'); ylabel('Phase in pi Radians'); title('Phase Response');

%% Example 3.15 frequency response by difference equation
% subplot(1,1,1)
% b = 1; a = [1,-0.8];
% n=(0:100);x = cos(0.05*pi*n);
% y = filter(b,a,x);
% subplot(2,1,1); stem(n,x);
% xlabel('n'); ylabel('x(n)'); title('Input sequence')
% subplot(2,1,2); stem(n,y);
% xlabel('n'); ylabel('y(n)'); title('Output sequence')

%% Example 3.16 matrix-vector multiplication form
% b = [0.0181, 0.0543, 0.0543, 0.0181]; % filter coefficient array b
% a = [1.0000, -1.7600, 1.1829, -0.2781]; % filter coefficient array a
% m = 0:length(b)-1; l = 0:length(a)-1; % index arrays m and l
% K = 500; k = 0:1:K; % index array k for frequencies
% w = pi*k/K; % [0, pi] axis divided into 501 points.
% num = b * exp(-1i*m'*w); % Numerator calculations. (b*m'*w):(1,4)*(4,1)*(1,501) --> 1*501, implement the sum. 
% den = a * exp(-1i*l'*w); % Denominator calculations
% H = num ./ den; % Frequency response
% magH = abs(H); angH = angle(H); % mag and phase responses
% subplot(2,1,1); plot(w/pi,magH); grid; axis([0,1,0,1])
% xlabel('frequency in pi units'); ylabel('|H|'); title('Magnitude Response');
% subplot(2,1,2); plot(w/pi,angH/pi); grid
% xlabel('frequency in pi units'); ylabel('Phase in pi Radians'); title('Phase Response');

%% Example 3.18 approximate the analog signal's analysis
% % analog signal 
% Dt = 0.00005; t = -0.005:Dt:0.005; xa = exp(-1000*abs(t)); 
% % Continuous-time Fourier Transform(approximate)
% Wmax = 2*pi*2000; K = 500; k = 0:1:500; W = k*Wmax/K; 
% Xa = xa * exp(-1i*t'*W) * Dt; % implement (3.35). (t*t'*W):(1,Dt)*(Dt,1)*(1,501) --> 1*501, implement the sum.
% Xa = real(Xa);
% W = [-fliplr(W), W(2:501)]; % Omega from -Wmax to Wmax
% Xa = [fliplr(Xa), Xa(2:501)]; % Xa over -Wmax to Wmax interval
% subplot(2,1,1); plot(t*1000,xa);
% xlabel('t in msec.'); ylabel('xa(t)');title('Analog Signal')
% subplot(2,1,2); plot(W/(2*pi*1000),Xa*1000);
% xlabel('Frequency in KHz'); ylabel('Xa(jW)*1000'); title('Continuous-time Fourier Transform')

%% Example 3.19 
%  % Analog Signal
% Dt = 0.00005; t = -0.005:Dt:0.005; xa = exp(-1000*abs(t));
% % Discrete-time Signal
% Ts = 0.0002; n = -25:1:25; x = exp(-1000*abs(n*Ts));
% % Discrete-time Fourier transform
% K = 500; k = 0:1:K; w = pi*k/K;
% X = x * exp(-1i*n'*w); 
% X = real(X);
% w = [-fliplr(w), w(2:K+1)]; X = [fliplr(X), X(2:K+1)];
% subplot(2,1,1);plot(t*1000,xa);
% xlabel('t in msec.'); ylabel('x1(n)'); title('Discrete Signal'); hold on
% stem(n*Ts*1000,x); gtext('Ts=0.2 msec'); hold off
% subplot(2,1,2);plot(w/pi,X);
% xlabel('Frequency in pi units'); ylabel('X1(w)'); title('Discrete-time Fourier Transform')

%% Example 3.21 reconstruct xa(t)
% % Discrete-time Signal x1(n)
% Ts = 0.0002; Fs = 1/Ts; n = -25:1:25; nTs = n*Ts; x = exp(-1000*abs(nTs));
% % Analog Signal reconstruction
% Dt = 0.00005; t = -0.005:Dt:0.005;
% xa = x * sinc(Fs*(ones(length(n),1)*t-nTs'*ones(1,length(t))));
% % check
% error = max(abs(xa - exp(-1000*abs(t))))

%% Example 3.22
% % Discrete-time Signal x2(n)
% Ts = 0.001; Fs = 1/Ts; n = -5:1:5; nTs = n*Ts; x = exp(-1000*abs(nTs));
% % Analog Signal reconstruction
% Dt = 0.00005; t = -0.005:Dt:0.005;
% xa = x * sinc(Fs*(ones(length(n),1)*t-nTs'*ones(1,length(t))));
% % check
% error = max(abs(xa - exp(-1000*abs(t))))

%% Example 3.23, stairs function --> ZOH, plot --> FOH
% % a) Discrete-time Signal x1(n) : Ts = 0.0002
% Ts = 0.0002; n = -25:1:25; nTs = n*Ts; x = exp(-1000*abs(nTs)); 
% % Plots
% subplot(2,1,1); stairs(nTs*1000,x);
% xlabel('t in msec.'); ylabel('xa(t)')
% title('Reconstructed Signal from x1(n) using zero-order-hold'); hold on
% stem(n*Ts*1000,x); hold off

% % b) Discrete-time Signal x2(n) : Ts = 0.001
% Ts = 0.001; n = -5:1:5; nTs = n*Ts; x = exp(-1000*abs(nTs));
% % Plots
% subplot(2,1,2); plot(nTs*1000,x);
% xlabel('t in msec.'); ylabel('xa(t)')
% title('Reconstructed Signal from x2(n) using zero-order-hold'); hold on
% stem(n*Ts*1000,x); hold off

%% Example 3.24 spline function --> cubic spline functions
% a) Discrete-time Signal x1(n): Ts = 0.0002
Ts = 0.0002; n = -25:1:25; nTs = n*Ts; x = exp(-1000*abs(nTs));
% Analog Signal reconstruction
Dt = 0.00005; t = -0.005:Dt:0.005; xa = spline(nTs,x,t);
% check
error = max(abs(xa - exp(-1000*abs(t))));

% b) Discrete-time Signal x2(n): Ts = 0.001
Ts = 0.001; n = -5:1:5; nTs = n*Ts; x = exp(-1000*abs(nTs));
% Analog Signal reconstruction
Dt = 0.00005; t = -0.005:Dt:0.005; xa = spline(nTs,x,t);
% check
error = max(abs(xa - exp(-1000*abs(t))));



















