function fval = emd2(gdm)

% same as edm, but with gdm passed

% ground distance matrix
f = gdm; %(F1, F2, Func);
% number of feature vectors
[m n] = size(f);

W1 = ones(m,1)./m;
W2 = ones(n,1)./n;

% inequality constraints
A1 = zeros(m, m * n);
A2 = zeros(n, m * n);
for i = 1:m
    for j = 1:n
        k = j + (i - 1) * n;
        A1(i, k) = 1;
        A2(j, k) = 1;
    end
end
A = [A1; A2];
b = [W1; W2];
% equality constraints
Aeq = ones(m + n, m * n);
beq = ones(m + n, 1) * min(sum(W1), sum(W2));
% lower bound
lb = zeros(1, m * n);
% linear programming
options = optimoptions('linprog','Display','off');
[x, fval] = linprog(f, A, b, Aeq, beq, lb, [], options);
fval = fval / sum(x);
end

%%%

function [f] = gdm(F1, F2, Func)
%
% GDM   Ground distance matrix between two signatures
%    [F] = GDM(F1, F2, FUNC) is the ground distance matrix between
%    two signatures whose feature vectors are given in F1 and F2.
%    FUNC is a function which computes the ground distance between
%    two feature vectors.
%
%    Example:
%    -------
%        f1 = [[100, 40, 22]; [211, 20, 2]; [32, 190, 150]; [2, 100, 100]];
%        f2 = [[0, 0, 0]; [50, 100, 80]; [255, 255, 255]];
%        ...
%        [f] = gdm(f1, f2, @gmf);
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
% number and length of feature vectors
[m a] = size(F1);
[n a] = size(F2);
% ground distance matrix
for i = 1:m
    for j = 1:n
        f(i, j) = Func(F1(i, 1:a), F2(j, 1:a));
    end
end
% gdm in column-vector form
f = f';
f = f(:);
end
