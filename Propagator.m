% Propagator
% P = 1/(a'QQ'a)
% Q = A2A1'pinv(A1A1')
clc; clear; close all;

%% Initial signal
% Propagator only work on many snapshots
T = 20;
snsrMap = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12];
tThetas = [5; 10];
[X, tA] = arrayModel(snsrMap, T, tThetas);
tS = pinv(tA' * tA) * tA' * X;

%% Initial scan grid
pThetas = -64:0.5:63.5;
pThetas = pThetas .';
pA = exp(1i * pi * (snsrMap - 1)* sind(pThetas).');
pA = pA ./ norm(pA(:, 1));

%% DBF
pDBF = zeros(size(pThetas));
for i = 1: length(pThetas)
    pDBF(i) = pinv(pA(:, i)' * pA(:, i)) * pA(:, i)' * X(:, 1);
end

%% Propagator Method
R = X * X' ./ T;
tN = length(tThetas);
A1 = R(1: tN, :);
A2 = R(tN + 1:end, :);
% A2 = P * A1;
P = A2 * A1' * pinv(A1 * A1');
Q = zeros(length(snsrMap), length(snsrMap) - tN);
Q(1: size(P', 1), 1: size(P', 2)) = P';
Q(tN + 1: end, :) = -eye(length(snsrMap) - tN);

mP = zeros(size(pThetas));
for i = 1: length(pThetas)
    mP(i) = real(1/ (pA(:, i)' * Q * Q' * pA(:, i)));
end

%% Figure;
figure;
plot(pThetas, 20 * log10(abs(pDBF) / max(abs(pDBF))));
hold on;
plot(pThetas, 20 * log10(real(mP) / max(real(mP))));
for i = 1: length(tThetas)
    xline(tThetas(i), '--g');
end
legend("DBF", "Propagator", "Golden")
grid on;
