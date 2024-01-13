% ESPRIT
clc; clear; close all;

%% Initial signal
% ESPRIT only work on many snapshots
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

%%ESPRIT
Ms = length(snsrMap) - 1;
X1 = X(1: Ms, :);
X2 = X(length(snsrMap) - Ms + 1: end, :);
XE = [X1; X2];
R = XE * XE' ./ T;
[U, S, V] = svd(R);
R = R - S(2 * Ms, 2 * Ms) * eye(2 * Ms);
[U, S, V] = svd(R);
Us = U(:, 1: length(tThetas));
Us1 = Us(1: Ms, :);
Us2 = Us(Ms + 1: end, :);
M = pinv(Us1) * Us2;
[V, D] = eig(M);
aThetas = asind(angle(diag(D).') / pi);

%% Figure;
figure;
plot(pThetas, 20 * log10(abs(pDBF) / max(abs(pDBF))));
hold on;
scatter(aThetas, zeros(length(aThetas), 1), '*')
for i = 1: length(tThetas)
    xline(tThetas(i), '--g');
end
legend("DBF", "ESPRIT", "Golden")
grid on;
