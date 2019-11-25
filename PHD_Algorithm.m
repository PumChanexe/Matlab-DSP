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
xn = awgn(x, SNR); % 信号添加白噪声
%%%% plot signal and noise %%%%
% figure(1)
% subplot(211);
% plot(n, x, 'b-'); 
% axis([0,Nx-1,-15,15]); title('信号x');
% xlabel('信号长度'); ylabel('x(n)'); grid
% subplot(212),
% plot(n, xn, 'r-');
% axis([0,Nx-1,-4, 4]); title('信号x + 白噪声');
% xlabel('噪声长度'); ylabel('w(n)'); grid
%% PHD Algorithm
rx = xcorr(xn, 'biased'); % 有偏自相关函数值求取
% rx = [6; 1.92705 + 1i*4.58522; -3.42705 + 1i*3.49541]; M = 2;
Rx = toeplitz(rx(Nx : Nx+M)); % Create toeplitz matrix
Rx = transpose(Rx); % 使用toeplitz函数的原因，转置一次，仿真出现频率反号
[V,D] = eig(Rx); % 返回特征值的对角矩阵 D 和矩阵 V，其列是对应的右特征向量，使得 A*V = V*D。
sigma = min(diag(D)); % diag：获取特征值矩阵的对角元素；min：找到最小值，即为白噪声功率
index = find(diag(D) == sigma); % 找到对应最小特征值的索引（下标）
vmin = V(:,index); % 最小特征值对应的噪声特征向量
z = roots(vmin); omega = angle(z); % 计算多项式vmin(z)的根，根的辐角即为复指数信号的频率
j = 0 : length(vmin) - 1; % 求解功率的线性方程组的索引
e = exp(-1i * omega * j); ev = e * V(: ,index+1 : M+1); % 计算信号特征矢量在在两个频率上的傅里叶变换的表达式

A = abs(ev).^2; % 计算特征向量（信号）的傅里叶变换的平方幅度在信号频率上的值
A = rot90(A); A = fliplr(A); % fliplr：将数组从左向右翻转
C = sum(D,2) - sigma; 
B = C(index + 1 : M+1); B = flipud(B); % sum(D,2)包含每一行总和的列向量；flipud：将数组从上向下翻转
P = linsolve(A,B); % 信号频率对应的功率
%%%% 对噪声特征矢量进行FFT变换，取模值平方之后做倒数，得到频率估计函数 %%%%
K =256; k = 0 : K-1; w = 2*pi*k/K;
F = fft(vmin, K); FEF = 1./(abs(F).^2 + eps);
FEF_dB = 10*log10(FEF);
%%%% plot %%%%
figure(2)
subplot(211)
plot(w, FEF);
set(gca,'XTick',(0 : 0.5*pi : 2*pi)); set(gca,'xtickLabel',{'0', '0.5π', 'π', '1.5π', '2π'}) % 设置刻度标签
xlabel('频率/rad'); ylabel('谱密度'); grid
title('频率估计函数'); 
subplot(212)
plot(w, FEF_dB);
set(gca,'XTick',(0 : 0.5*pi : 2*pi)); set(gca,'xtickLabel',{'0', '0.5π', 'π', '1.5π', '2π'}) % 设置刻度标签
xlabel('频率/rad'); ylabel('谱密度/dB'); grid
title('频率估计函数/dB'); 









