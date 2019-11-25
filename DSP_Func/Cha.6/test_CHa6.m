clc; clear; close;
%% Example 6.1 DIRECT form to CASCADE form
% b = [1 -3 11 -27 18]; a = [16 12 2 -4 -1];
% [b0,B,A] = dir2cas(b,a);

%% Example 6.2 DIRECT form to PARALLEL form
% b = [1 -3 11 -27 18]; a = [16 12 2 -4 -1];
% [C,B,A] = dir2par(b,a);

%% Example 6.3 combine form 
% C0=0; B1=[2 4;3 1]; A1=[1 1 0.9; 1 0.4 -0.4];
% B2=[0.5 0.7;1.5 2.5;0.8 1]; A2=[1 -1 0.8;1 0.5 0.5;1 0 -0.5];
% [b1,a1]=par2dir(C0,B1,A1);
% [b2,a2]=par2dir(C0,B2,A2);
% b=conv(b1,b2); % Overall direct form numerator
% a=conv(a1,a2); % Overall direct form denominator
% [b0,Bc,Ac]=dir2cas(b,a); % Overall cascade form
% [C0,Bp,Ap]=dir2par(b,a); % Overall parallel form

%% Example 6.5 cascade with linear-phase form
% b = [1,0,0,0,16+1/16,0,0,0,1]; 
% broots = roots(b);
% B1 = real(poly([broots(1),broots(2),broots(5),broots(6)]));
% B2 = real(poly([broots(3),broots(4),broots(7),broots(8)]));

%% Example 6.6 frequency sample form
% h = [1,2,3,2,1]/9; [C,B,A] = dir2fs(h); 

%% Example 6.7
% M = 32; alpha = (M-1)/2;
% magHk = [1,1,1,0.5,zeros(1,25),0.5,1,1];
% k1 = 0:15; k2 = 16:M-1;
% angHk = [-alpha*(2*pi)/M*k1, alpha*(2*pi)/M*(M-k2)]; % (6.11)
% H = magHk.*exp(1i*angHk); h = real(ifft(H,M)); 
% [C,B,A] = dir2fs(h);

%% Example 6.8 FIR direct form to lattice form conversion
% b = [2, 13/12, 5/4, 2/3]; 
% K = dir2latc(b);
% [delta,n] = impseq(0,0,3); format long; hdirect = filter(b,1,delta); 
% hlattice = latcfilt(K,delta); 
% error = (hdirect - hlattice);

%% Example 6.10 IIR filter to lattice-ladder structure
% b = [1,2,2,1]; a = [1,13/24,5/8,1/3]; [K,C] = dir2ladr(b,a);
% [x,n] = impseq(0,0,7); format long; 
% hdirect = filter(b,a,x); 
% hladder = ladrfilt(K,C,x);
% error = hdirect - hladder;

%% Example 6.11 
x = -7:7;
y = OnesComplement(x,4);






