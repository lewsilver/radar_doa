% FOCUSS
% W = diag(S)
% T = A*W
% S = inv(T' * T) * T' * X
clc; clear; close all;

%% Initial signal
% FOCUSS can work on one snapshot
T = 1;
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

%% FOCUSS
maxIter = 20;
pS = pinv(pA' * pA) * pA' * X;
for i = 1: maxIter
    W = diag(pS);
    T = pA * W;
    pS = pinv(T' * T) * T' * X;
    figure(1);
    plot(abs(pS));
end

%% Figure;
figure;
plot(pThetas, 20 * log10(abs(pDBF) / max(abs(pDBF))));
hold on;
plot(pThetas, 20 * log10(abs(pS) / max(abs(pS))));
for i = 1: length(tThetas)
    xline(tThetas(i), '--g');
end
legend("DBF", "FOCUSS", "Golden")
grid on;