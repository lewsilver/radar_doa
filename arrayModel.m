function [X, tA] = arrayModel(snsrMap, T, target)
    % X = As + e
    % snsrMap: Sensor array, d = lambda / 2
    % T      : Snapshot
    % target : Targets at angle thetas

    tA = exp(1i * pi * (snsrMap - 1) * sind(target).');
    tA = tA ./ norm(tA(:, 1));

    pd = makedist('Normal','mu',30,'sigma',5);
    s = @() random(pd, size(target));

    noiseSNR = 15;
    X = zeros(length(snsrMap), T);
    for i = 1: T
        X(:, i) = awgn(tA * s(), noiseSNR);
    end

end

