function x_test_emd

rng(100);

f1 = [[100, 40, 22]; [211, 20, 2]; [32, 190, 150]; [2, 100, 100]];
f2 = [[0, 0, 0]; [50, 100, 80]; [255, 255, 255]];
w1 = [0.4; 0.3; 0.2; 0.1];
w2 = [0.5; 0.3; 0.2];

[x fval] = emd(f1, f2, w1, w2, @gdf);

fval

Ns = 20;

fvs = [];

for Ns = 20:20:100
    
    f1 = randn(100,10);
    f2 = randn(Ns,10);
    w1 = 0.01 * ones(100,1);
    w2 = (1/Ns) * ones(Ns,1);
    
    [x fval] = emd(f1, f2, w1, w2, @gdf);
    
    fvs = [fvs fval];
    
end

fvs

end

function [E] = gdf(V1, V2)
%
% GDF   Ground distance between two vectors
%    [E] = GDF(F1, F2) is the ground distance between two feature vectors.
%
%    Example:
%    -------
%        v1 = [100, 40, 22];
%        v2 = [50, 100, 80];
%        ...
%        [e] = gdf(v1, v2);
%        ...
%
%    This file and its content belong to Ulas Yilmaz.
%    You are welcome to use it for non-commercial purposes, such as
%    student projects, research and personal interest. However,
%    you are not allowed to use it for commercial purposes, without
%    an explicit written and signed license agreement with Ulas Yilmaz.
%    Berlin University of Technology, Germany 2006.
%    http://www.cv.tu-berlin.de/~ulas/RaRF
%
E = norm(V1 - V2, 2);
end
