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
% axis([0,Nx-1,-15,15]); title('�ź�x');
% xlabel('�źų���'); ylabel('x(n)'); grid
% subplot(212),
% plot(n, xn, 'r-');
% axis([0,Nx-1,-4, 4]); title('�ź�x + ������');
% xlabel('��������'); ylabel('w(n)'); grid
%% MUSIC Algorithm
rx = xcorr(xn, 'coeff'); % Autocorrelation function of data
Rx = toeplitz(rx(Nx : Nx+N-1)); % Autocorrelation matrix of data
Rx = transpose(Rx);
[V, D] = eig(Rx); % Engendecomposition of data autocorrelation matrix
%%%% ���Ƶ�ʹ��ƺ���z�任�ĸ������ķ��Ǽ�Ϊ�ź�Ƶ�� %%%%
d = 0; % Routing the Polynomial D(z)
for k = 1 : N-M
    v = V(:, k);
    v1 = flipud(v);
    d = d + conv(v, conj(v1));
end
roots_d = roots(d);
roots_d1 = roots_d(abs(roots_d) < 1);
[not_used, index] = sort(1 - abs(roots_d1)); % sort:�����ҵ���ӽ���λԲ����������
sorted_roots = roots_d1(index);
omega = angle(sorted_roots(1:M)); % ���Ǽ���ָ���ź�Ƶ��
lambda_all = diag(D); sigma_all = lambda_all(1:N-M);
sigma = sum(sigma_all)/(N-M); % noise variance
%%%% ���źŹ��� %%%%
j = 0 : N - 1;
e = exp(-1i * omega * j); 
ev = e * V(: ,N-M+1 : N);
A = abs(ev).^2; 
A = rot90(A); A = fliplr(A); % fliplr��������������ҷ�ת
C = lambda_all - sigma; 
B = C(N-M+1:N); 
B = flipud(B); % sum(D,2)����ÿһ���ܺ͵���������flipud��������������·�ת
P = linsolve(A,B);
%%%% Ƶ�ʹ��ƺ��� %%%%
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
set(gca,'XTick',(0 : 0.5*pi : 2*pi)); set(gca,'xtickLabel',{'0', '0.5��', '��', '1.5��', '2��'}) % ���ÿ̶ȱ�ǩ
xlabel('Ƶ��/rad'); ylabel('���ܶ�'); grid
title('Ƶ�ʹ��ƺ���'); 
subplot(212)
plot(w, FEF_dB); 
set(gca,'XTick',(0 : 0.5*pi : 2*pi)); set(gca,'xtickLabel',{'0', '0.5��', '��', '1.5��', '2��'}) % ���ÿ̶ȱ�ǩ
xlabel('Ƶ��/rad'); ylabel('���ܶ�/dB'); grid
title('Ƶ�ʹ��ƺ���/dB'); 








