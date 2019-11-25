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
axis([0,1000,-15,15]); title('���ڶ����Իع�ģ��ģ�Ͳ������ź�x');
xlabel('�źų���'); ylabel('x(n)'); grid
subplot(212),
plot(n, v, 'r-');
axis([0,1000,-4, 4]); title('��˹����������');
xlabel('��������'); ylabel('w(n)'); grid
%% LSL Algorithm
%%%%%%%%%%%%%%%%% LSL Algorithm %%%%%%%%%%%%%%%%%%%%
%%%% Initiation %%%%
ef = zeros(L+1,N); eb = zeros(L+1,N); % 1, 2, 3�±꣬��Ӧ2���˲���,���е�ǰ�򡢺���Ԥ������ʼ�������㡣
delta = zeros(L+1,N); % ��ʼ�����е�delta������Ϊ�㡣
%%%% D = 0.8;
D = [0.1 1.0 10.0];
% epsilonf = D*ones(L+1,1); epsilonb = D*ones(L+1,1);
kf = zeros(L+1,N); kb = zeros(L+1,N); 
gama = ones(L+1,1); % ����gama���������г�ʼ��Ϊ1
%%%% Iteration %%%%
for i = 1:length(D)
    epsilonf = D(i)*ones(L+1,1); epsilonb = D(i)*ones(L+1,1);
    for n = 2:N % n��ʱ�䣬��ʱ���������
        eb(1,n) = x(n); ef(1,n) = x(n); % ǰ��Ԥ�����ͺ���Ԥ�����ĵ�һ��ϵ������ÿ��ʱ�̶���������x(n)
        epsilonf(1,n) = epsilonf(1,n-1) + (x(n))^2;
        epsilonb(1,n) = epsilonf(1,n);
        gama(1,n) = 1; % ��һ��gamaϵ����������ʱ�̶�Ϊ1
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
axis([0,1000,-2, 2]); xlabel('�źų���n'); ylabel('Ȩֵϵ��w1'); grid
legend('D=0.1','D=1.0','D=10.0'); title('LSL�㷨�²�ͬDֵ��Ӧ��w1��������');
subplot(212)
plot(n,LSL_w2(:, :, 1), 'r-', n,LSL_w2(:, :, 2), 'b-', n,LSL_w2(:, :, 3), 'k-');
axis([0,1000,-2, 2]); xlabel('�źų���n'); ylabel('Ȩֵϵ��w2'); grid
legend('D=0.1','D=1.0','D=10.0'); title('LMS�㷨�²�ͬDֵ��Ӧ��w2��������');
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
axis([0,1000,-2, 2]); xlabel('�źų���n'); ylabel('Ȩֵϵ��w1'); grid
legend('��=0.001','��=0.002','��=0.004'); title('LMS�㷨�²�ͬ��ֵ��Ӧ��w1��������');
subplot(212)
plot(n, LMS_w2(:, :, 1), 'r-', n, LMS_w2(:, :, 2), 'b-', n, LMS_w2(:, :, 3), 'k-');
axis([0,1000,-2, 2]); xlabel('�źų���n'); ylabel('Ȩֵϵ��w2'); grid
legend('��=0.001','��=0.002','��=0.004'); title('LMS�㷨�²�ͬ��ֵ��Ӧ��w2��������');
%% LMS vs LSL
%%%%%%%%%%%%%%%%% LMS vs. LSL on convergense %%%%%%%%%%%%%%%%%
%%%% LSL : D=1.0, LMS : mu = 0.005 %%%%
figure(4)
subplot(211)
n = 1:N;
plot(n, LMS_w1(:, :, 3), 'r-', n, LSL_w1(:, :, 2), 'b-');
axis([0,1000,-2, 2]); xlabel('�źų���n'); ylabel('Ȩֵϵ��w1'); grid
legend('��=0.005','D=1.0'); title('LMS vs LSL');
subplot(212)
plot(n, LMS_w2(:, :, 3), 'r-', n, LSL_w2(:, :, 2), 'b-');
axis([0,1000,-2, 2]); xlabel('�źų���n'); ylabel('Ȩֵϵ��w2'); grid
legend('��=0.001','D=1.0'); title('LMS vs LSL');