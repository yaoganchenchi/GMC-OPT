
function [xg, wg] = gauss_quad_new(ng)
% gaussLegendreQuad - Compute Gauss-Legendre quadrature nodes and weights
% using the Golub-Welsch algorithm.
%
% Inputs:
%   n - Quadrature order (number of points)
% Outputs:
%   xg - Gauss-Legendre nodes (in [-1, 1])
%   wg - Corresponding weights

    % Form symmetric tridiagonal Jacobi matrix
    i = 1:ng-1;
    a = zeros(ng,1);                            % diagonal entries (all 0)
    b = i ./ sqrt(4*i.^2 - 1);                % sub-diagonal entries

    % Construct Jacobi matrix
    J = diag(a) + diag(b,1) + diag(b,-1);

    % Compute eigenvalues (nodes) and eigenvectors
    [V, D] = eig(J);

    % Nodes are eigenvalues of J
    xg = diag(D);

    % Weights are square of first row of eigenvectors, scaled by 2
    wg = 2 * (V(1,:).^2)';

    % Sort the nodes and weights
    [xg, idx] = sort(xg);
    wg = wg(idx);
end


