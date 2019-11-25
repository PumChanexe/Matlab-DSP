clc; clear; close;
%% Example 4.4
% x1 = [2,3,4]; x2 = [3,4,5,6]; x3 = conv(x1,x2);

%% Example 4.5
% x1 = [1,2,3]; n1 = [-1:1]; x2 = [2,4,3,5]; n2 = [-2:1];
% [x3,n3] = conv_m(x1,n1,x2,n2)
%% Problem 4.10 modify the deconv function to obtain sample index
%% Example 4.6
% b = [0,0,0,0.25,-0.5,0.0625]; a = [1,-1,0.75,-0.25,0.0625];
% [delta,n]=impseq(0,0,6);
% stem(n,delta);
% x_check = filter(b,a,delta); % check sequence, filter output is x(n)
% x_ori = ((n-2).*(1/2).^(n-2).*cos(pi*(n-2)/3)).*stepseq(2,0,7); % original sequence

%% Example 4.8, residuez function
% b = [0,1]; a = [3,-4,1]; [R,p,C] = residuez(b,a)

%% Example 4.9 
% b = 1; a = poly([0.9,0.9,-0.9]);
% [R,p,C] = residuez(b,a);
% [delta,n] = impseq(0,0,7); x_check = filter(b,a,delta); % filter implement the inverse z trans
% x_answer = (0.75)*(0.9).^n + (0.5)*n.*(0.9).^n + (0.25)*(-0.9).^n; % answer sequence
% error = x_check - x_answer;

%% Example 4.10
% b = [1,0.4*sqrt(2)]; a = [1,-0.8*sqrt(2),0.64]; [R,p,C] = residuez(b,a);
% Mp = (abs(p))'; % pole magnitudes
% Ap = (angle(p))'/pi; % pole angles in pi units
% [delta,n] = impseq(0,0,6);
% x_check = filter(b,a,delta);
% x_answer = ((0.8).^n).*(cos(pi*n/4)+2*sin(pi*n/4));
% error = x_check - x_answer;

%% Example 4.11
% b = [1, 0]; a = [1, -0.9]; 
% zplane(b,a);
% -----------------
% [H,w] = freqz(b,a,100); 
% magH = abs(H); phaH = angle(H);
% --------------
% [H,w] = freqz(b,a,200,'whole');
% magH = abs(H(1:101)); phaH = angle(H(1:101));
% -----------------------
% subplot(2,1,1); plot(w/pi,magH); grid;
% xlabel('frequency in pi units'); ylabel('Magnitude');
% title('Magnitude Response');
% subplot(2,1,2);plot(w/pi,phaH/pi);grid
% xlabel('frequency in pi units'); ylabel('Phase in pi units');
% title('Phase Response')

%% Example 4.12, transfer function, difference equation, impulse response

%% Example 4.13, causal LTI system

%% Example 4.14 ** not solve
% n = (0:7); b = [2,-9/4,1/2]; a = poly([1/2,1,1/4]); 
% x = (1/4).^n; 
% % xic = [1,-2];
% format long; 
% y1 = filter(b,a,x);
% y2 = (1/3)*(1/4).^n+(1/2).^n+(2/3)*ones(1,8); % MATLAB Check

%% Example 4.15 
% b = [1,1,1]/3; a = [1,-0.95,0.9025];
% Y = [-2,-3]; X = [1,1]; xic = filtic(b,a,Y,X);
% bxplus = [1,-0.5]; axplus = [1,-1,1]; % X(z) transform coeff.
% ayplus = conv(a,axplus); % Denominator of Yplus(z)
% byplus = conv(b,bxplus)+conv(xic,axplus); % Numerator of Yplus(z)
% [R,p,C] = residuez(byplus,ayplus);
% Mp = abs(p); Ap = angle(p)/pi; % Polar form














