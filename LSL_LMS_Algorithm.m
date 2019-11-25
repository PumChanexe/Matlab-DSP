close all;
clc;
clear variables;
%%%%%%%%%%%%%%%%% Construct input and desired output %%%%%%%%%%%%%%%%%
%%%% x:input signal; d:expect signal; %%%%
a1 = 1.558; a2 = -0.81; N = 1000; L = 2;
v = randn(1, N);
x = zeros(1, N); x(1) = v(1); x(2) = a1*x(1) + v(2);
d = zeros(1, N); d(1) = x(1); d(2) = x(2);
for n = 3:N
    x(n) = a1*x(n-1) + a2 * x(n-2) + v(n); % only the first column is assigned.
    d(n) = x(n);
end
%%%% plot the signal and noise %%%%
figure(1)
n = 1 : N;
subplot(211);
plot(n, x, 'b-'); 
axis([0,1000,-15,15]); title('基于二阶自回归模型模型产生的信号x');
xlabel('信号长度'); ylabel('x(n)'); grid
subplot(212),
plot(n, v, 'r-');
axis([0,1000,-4, 4]); title('高斯白噪声序列');
xlabel('噪声长度'); ylabel('w(n)'); grid
%% LSL Algorithm
%%%%%%%%%%%%%%%%% LSL Algorithm %%%%%%%%%%%%%%%%%%%%
%%%% Initiation %%%%
ef = zeros(L+1,N); eb = zeros(L+1,N); % 1, 2, 3下标，对应2阶滤波器,所有的前向、后向预测误差初始化都是零。
delta = zeros(L+1,N); % 初始化所有的delta参数都为零。
%%%% D = 0.8;
D = [0.1 1.0 10.0];
% epsilonf = D*ones(L+1,1); epsilonb = D*ones(L+1,1);
kf = zeros(L+1,N); kb = zeros(L+1,N); 
gama = ones(L+1,1); % 创建gama参数，所有初始化为1
%%%% Iteration %%%%
for i = 1:length(D)
    epsilonf = D(i)*ones(L+1,1); epsilonb = D(i)*ones(L+1,1);
    for n = 2:N % n是时间，按时间迭代计算
        eb(1,n) = x(n); ef(1,n) = x(n); % 前向预测误差和后向预测误差的第一个系数，在每个时刻都等于输入x(n)
        epsilonf(1,n) = epsilonf(1,n-1) + (x(n))^2;
        epsilonb(1,n) = epsilonf(1,n);
        gama(1,n) = 1; % 第一个gama系数，在所有时刻都为1
        for m = 1:L
            delta(m+1,n) = delta(m+1,n-1) + eb(m,n-1)*ef(m,n)/gama(m,n-1);
            ef(m+1,n) = ef(m,n) - delta(m+1,n)*eb(m,n-1)/epsilonb(m,n-1);
            eb(m+1,n) = eb(m,n-1) - delta(m+1,n)*ef(m,n)/epsilonf(m,n);
            epsilonf(m+1,n) = epsilonf(m,n) - (delta(m+1,n))^2/epsilonb(m,n-1);
            epsilonb(m+1,n) = epsilonb(m,n-1) - (delta(m+1,n))^2/epsilonf(m,n);
            gama(m+1,n-1) = gama(m,n-1) - (eb(m,n-1))^2/epsilonb(m,n-1);
            kf(m+1,n) = delta(m+1,n)/epsilonf(m,n);
            kb(m+1,n) = delta(m+1,n)/epsilonb(m,n-1);
        end
    end
    k1 = kf(2,:); k1 = circshift(k1, 999); k2 = kf(3,:); 
    k3 = kb(2,:); k4 = kb(3,:);
    LSL_w1(:, :, i) = k3 - k1.*k4;
    LSL_w2(:, :, i) = k4;
end
%%%% Plot the coefficients %%%%
figure(2)
n = 1:N;
subplot(211)
plot(n,LSL_w1(:, :, 1), 'r-', n,LSL_w1(:, :, 2), 'b-', n,LSL_w1(:, :, 3), 'k-');
axis([0,1000,-2, 2]); xlabel('信号长度n'); ylabel('权值系数w1'); grid
legend('D=0.1','D=1.0','D=10.0'); title('LSL算法下不同D值对应的w1收敛曲线');
subplot(212)
plot(n,LSL_w2(:, :, 1), 'r-', n,LSL_w2(:, :, 2), 'b-', n,LSL_w2(:, :, 3), 'k-');
axis([0,1000,-2, 2]); xlabel('信号长度n'); ylabel('权值系数w2'); grid
legend('D=0.1','D=1.0','D=10.0'); title('LMS算法下不同D值对应的w2收敛曲线');
%% LMS Algorithm
%%%%%%%%%%%%%%%%% LMS Algorithm %%%%%%%%%%%%%%%%%
%%%% mu = 0.005;
mu3 = [0.001 0.002 0.005];
w = zeros(L,1); theta = zeros(L,1);
e = zeros(size(x)); q = x(:); 
w1w2 = zeros(length(w),N);
%%%% find optimal weight vector %%%%
for i = 1:length(mu3)
    mu_x = mu3(i);
    for k = 3:N
        theta = q(k-1:-1:k-L);
        e(k) = d(k) - w'*theta;
        w = w + 2*mu_x*e(k)*theta;
        w1w2(:,k) = w(:);
    end
    LMS_w1(:,:,i) = w1w2(1, :);
    LMS_w2(:,:,i) = w1w2(2, :);
end
%%%% Plot the coefficients %%%%
figure(3)
subplot(211)
n = 1:N;
plot(n, LMS_w1(:, :, 1), 'r-', n, LMS_w1(:, :, 2), 'b-', n, LMS_w1(:, :, 3), 'k-');
axis([0,1000,-2, 2]); xlabel('信号长度n'); ylabel('权值系数w1'); grid
legend('μ=0.001','μ=0.002','μ=0.004'); title('LMS算法下不同μ值对应的w1收敛曲线');
subplot(212)
plot(n, LMS_w2(:, :, 1), 'r-', n, LMS_w2(:, :, 2), 'b-', n, LMS_w2(:, :, 3), 'k-');
axis([0,1000,-2, 2]); xlabel('信号长度n'); ylabel('权值系数w2'); grid
legend('μ=0.001','μ=0.002','μ=0.004'); title('LMS算法下不同μ值对应的w2收敛曲线');
%% LMS vs LSL
%%%%%%%%%%%%%%%%% LMS vs. LSL on convergense %%%%%%%%%%%%%%%%%
%%%% LSL : D=1.0, LMS : mu = 0.005 %%%%
figure(4)
subplot(211)
n = 1:N;
plot(n, LMS_w1(:, :, 3), 'r-', n, LSL_w1(:, :, 2), 'b-');
axis([0,1000,-2, 2]); xlabel('信号长度n'); ylabel('权值系数w1'); grid
legend('μ=0.005','D=1.0'); title('LMS vs LSL');
subplot(212)
plot(n, LMS_w2(:, :, 3), 'r-', n, LSL_w2(:, :, 2), 'b-');
axis([0,1000,-2, 2]); xlabel('信号长度n'); ylabel('权值系数w2'); grid
legend('μ=0.001','D=1.0'); title('LMS vs LSL');