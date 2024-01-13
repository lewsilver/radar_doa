% OMP
clc; clear; close all;

%% Initial signal
% OMP can work on one snapshot
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

%% OMP
maxIter = length(tThetas);
index = zeros(1, maxIter);
oP = zeros(length(snsrMap), maxIter);
oX = X;
for i = 1: maxIter
    pS = pA' * oX;
    [V, I] = max(abs(pS));
    index(i) = I;
    oP(:, i) = pA(:, I);
    oX = oX - oP * pinv(oP' * oP) * oP' * oX;
    figure(1);
    clf;
    plot(abs(pS));
end

%% Figure;
figure;
plot(pThetas, 20 * log10(abs(pDBF) / max(abs(pDBF))));
hold on;
scatter(pThetas(index), zeros(length(index), 1), '*')
for i = 1: length(tThetas)
    xline(tThetas(i), '--g');
end
legend("DBF", "OMP", "Golden")
grid on;