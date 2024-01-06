% DML
% P = tr(A(A'*A)^-1A'R)
clc; clear; close all;

%% Initial signal
% DML work on one or more snapshots
% This work just work for only two targets
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

%% DML
R = X * X' ./ T;
maxPower = 0;
aThetas = [];
for i = 1: length(pThetas)
    eA = zeros(length(snsrMap), 2);
    eA(:, 1) = pA(:, i);
    for j = 1: length(pThetas)
        eA(:, 2) = pA(:, j);
        ePower = trace(eA * pinv(eA' * eA) * eA' * R);
        if ePower > maxPower
            maxPower = ePower;
            aThetas = [i, j];
        end
    end
end

%% Figure;
figure;
plot(pThetas, 20 * log10(abs(pDBF) / max(abs(pDBF))));
hold on;
scatter(pThetas(aThetas), zeros(length(aThetas), 1), '*')
for i = 1: length(tThetas)
    xline(tThetas(i), '--g');
end
legend("DBF", "DML", "Golden")
grid on;