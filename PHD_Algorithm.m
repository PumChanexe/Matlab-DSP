close all;
clc;
clear variables;
% Pisarenko Harmonic Decomposition
% rx:autocorrelation function of a random process
% M:The number of complex exponentials
% K:The length of FFT
Nx = 25;
n = (0 : Nx-1); M = 2; SNR = 30;
x = exp(1i*2*pi*0.5*n) + exp(1i*(2*pi*0.52*n + pi/4)); 
xn = awgn(x, SNR); % �ź���Ӱ�����
%%%% plot signal and noise %%%%
% figure(1)
% subplot(211);
% plot(n, x, 'b-'); 
% axis([0,Nx-1,-15,15]); title('�ź�x');
% xlabel('�źų���'); ylabel('x(n)'); grid
% subplot(212),
% plot(n, xn, 'r-');
% axis([0,Nx-1,-4, 4]); title('�ź�x + ������');
% xlabel('��������'); ylabel('w(n)'); grid
%% PHD Algorithm
rx = xcorr(xn, 'biased'); % ��ƫ����غ���ֵ��ȡ
% rx = [6; 1.92705 + 1i*4.58522; -3.42705 + 1i*3.49541]; M = 2;
Rx = toeplitz(rx(Nx : Nx+M)); % Create toeplitz matrix
Rx = transpose(Rx); % ʹ��toeplitz������ԭ��ת��һ�Σ��������Ƶ�ʷ���
[V,D] = eig(Rx); % ��������ֵ�ĶԽǾ��� D �;��� V�������Ƕ�Ӧ��������������ʹ�� A*V = V*D��
sigma = min(diag(D)); % diag����ȡ����ֵ����ĶԽ�Ԫ�أ�min���ҵ���Сֵ����Ϊ����������
index = find(diag(D) == sigma); % �ҵ���Ӧ��С����ֵ���������±꣩
vmin = V(:,index); % ��С����ֵ��Ӧ��������������
z = roots(vmin); omega = angle(z); % �������ʽvmin(z)�ĸ������ķ��Ǽ�Ϊ��ָ���źŵ�Ƶ��
j = 0 : length(vmin) - 1; % ��⹦�ʵ����Է����������
e = exp(-1i * omega * j); ev = e * V(: ,index+1 : M+1); % �����ź�����ʸ����������Ƶ���ϵĸ���Ҷ�任�ı��ʽ

A = abs(ev).^2; % ���������������źţ��ĸ���Ҷ�任��ƽ���������ź�Ƶ���ϵ�ֵ
A = rot90(A); A = fliplr(A); % fliplr��������������ҷ�ת
C = sum(D,2) - sigma; 
B = C(index + 1 : M+1); B = flipud(B); % sum(D,2)����ÿһ���ܺ͵���������flipud��������������·�ת
P = linsolve(A,B); % �ź�Ƶ�ʶ�Ӧ�Ĺ���
%%%% ����������ʸ������FFT�任��ȡģֵƽ��֮�����������õ�Ƶ�ʹ��ƺ��� %%%%
K =256; k = 0 : K-1; w = 2*pi*k/K;
F = fft(vmin, K); FEF = 1./(abs(F).^2 + eps);
FEF_dB = 10*log10(FEF);
%%%% plot %%%%
figure(2)
subplot(211)
plot(w, FEF);
set(gca,'XTick',(0 : 0.5*pi : 2*pi)); set(gca,'xtickLabel',{'0', '0.5��', '��', '1.5��', '2��'}) % ���ÿ̶ȱ�ǩ
xlabel('Ƶ��/rad'); ylabel('���ܶ�'); grid
title('Ƶ�ʹ��ƺ���'); 
subplot(212)
plot(w, FEF_dB);
set(gca,'XTick',(0 : 0.5*pi : 2*pi)); set(gca,'xtickLabel',{'0', '0.5��', '��', '1.5��', '2��'}) % ���ÿ̶ȱ�ǩ
xlabel('Ƶ��/rad'); ylabel('���ܶ�/dB'); grid
title('Ƶ�ʹ��ƺ���/dB'); 









