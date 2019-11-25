% exp3.9.5 : MUSIC
Nx = 25; n = 0 : Nx - 1; % length of data
M = 2; SNR = 30;% 
N = 12; % Size of data autocorrelation matrix
x = exp(1i*2*pi*0.5*n) + exp(1i*(2*pi*0.52*n + pi/4)); 
xn = awgn(x, SNR); 
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
%% MUSIC Algorithm
rx = xcorr(xn, 'coeff'); % Autocorrelation function of data
Rx = toeplitz(rx(Nx : Nx+N-1)); % Autocorrelation matrix of data
Rx = transpose(Rx);
[V, D] = eig(Rx); % Engendecomposition of data autocorrelation matrix
%%%% 求解频率估计函数z变换的根，根的辐角即为信号频率 %%%%
d = 0; % Routing the Polynomial D(z)
for k = 1 : N-M
    v = V(:, k);
    v1 = flipud(v);
    d = d + conv(v, conj(v1));
end
roots_d = roots(d);
roots_d1 = roots_d(abs(roots_d) < 1);
[not_used, index] = sort(1 - abs(roots_d1)); % sort:升序，找到最接近单位圆的零点的索引
sorted_roots = roots_d1(index);
omega = angle(sorted_roots(1:M)); % 辐角即复指数信号频率
lambda_all = diag(D); sigma_all = lambda_all(1:N-M);
sigma = sum(sigma_all)/(N-M); % noise variance
%%%% 求信号功率 %%%%
j = 0 : N - 1;
e = exp(-1i * omega * j); 
ev = e * V(: ,N-M+1 : N);
A = abs(ev).^2; 
A = rot90(A); A = fliplr(A); % fliplr：将数组从左向右翻转
C = lambda_all - sigma; 
B = C(N-M+1:N); 
B = flipud(B); % sum(D,2)包含每一行总和的列向量；flipud：将数组从上向下翻转
P = linsolve(A,B);
%%%% 频率估计函数 %%%%
[y, i] = sort(diag(D));
Px = 0; K =1024; k = 0 : K-1; w = 2*pi*k/K;
for j = 1:N-M
    Px = Px + abs(fft(V(:,i(j)), K)).^2;
end
Px = flipud(Px); FEF = 1./Px; 
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








