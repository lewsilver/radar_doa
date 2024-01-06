% CAPON
% P = 1/(a'R^-1a)
clc; clear; close all;

%% Initial signal
% Capon only work on many snapshots
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

%% Capon
% Capon result is the power, so just get real of the result 
R = X * X' ./ T;
invR = pinv(R);
cP = zeros(size(pThetas));
for i = 1: length(pThetas)
    cP(i) = 1 / (pA(:, i)' * invR * pA(:, i));
end

%% Figure;
figure;
plot(pThetas, 20 * log10(abs(pDBF) / max(abs(pDBF))));
hold on;
plot(pThetas, 20 * log10(real(cP) / max(real(cP))));
for i = 1: length(tThetas)
    xline(tThetas(i), '--g');
end
legend("DBF", "Capon", "Golden")
grid on;
