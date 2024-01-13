% ISTA
clc; clear; close all;

%% Initial signal
% ISTA can work on one snapshot
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

%% ISTA
soft = @(y, T) (max(abs(y) - T, 0)) .* sign(y);
[V, D] = eig(pA' * pA);
alpha = max(diag(D));
lambda = 30;
T = lambda / (2 * alpha);
maxIter = 1000;
pS = zeros(size(pThetas));
for i = 1: maxIter
    S = pS +  (1 / alpha) * pA' * (X - pA * pS);
    pS = soft(S, T);
    figure(1);
    clf;
    plot(abs(S));
    hold on;
    plot(abs(pS));
end

%% Figure;
figure;
plot(pThetas, 20 * log10(abs(pDBF) / max(abs(pDBF))));
hold on;
plot(pThetas, 20 * log10(abs(S) / max(abs(S))));
for i = 1: length(tThetas)
    xline(tThetas(i), '--g');
end
legend("DBF", "ISTA", "Golden")
grid on;