% IAA
% P = a'R^-1X/(a'R^-1a)
clc; clear; close all;

%% Initial signal
% IAA can work on one or more snapshots, this work just show one snapshot
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

%% IAA
R = X * X' ./ T;
maxIter = 16;
P = diag(abs(pDBF).^2);
for i = 1: maxIter
    Q = R - pA * P * pA';
    invQ = pinv(Q);
    aS = zeros(size(pThetas));
    for j = 1: length(pThetas)
        aS(j) = pA(:, j)' * invQ * X / (pA(:, j)' * invQ * pA(:, j));
    end
    figure(2);
    plot(abs(aS));
    P = diag(abs(aS) .^2);
end

%% Figure;
figure;
plot(pThetas, 20 * log10(abs(pDBF) / max(abs(pDBF))));
hold on;
plot(pThetas, 20 * log10(real(P) / max(real(P))));
for i = 1: length(tThetas)
    xline(tThetas(i), '--g');
end
legend("DBF", "IAA", "Golden")
grid on;
