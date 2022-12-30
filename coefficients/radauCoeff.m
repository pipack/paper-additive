function [Bi, Be] = radauCoeff(q, star)
    if(nargin == 1)
        star = false;
    end

    Bi = zeros(q);
    Be = zeros(q);
    z  = nodes(q);
    Bi(:,2:end) = transpose(quadW(z(2:end), z(1) * ones(q,1), z));
    if(star)
        Be = transpose(quadW(z, z(q) * ones(q,1), z + 2));
    else
        Be(:,2:end) = transpose(quadW(z(2:end), z(q) * ones(q,1), z + 2));
    end
end

function w = quadW(x, a, b)
%QUADW Returns quadrature weights for the nodes x where
%
%   \int^b(i)_a(i) f(x) dx \approx \sum_{j=1}^n w(j,i) f(x(i))
%
% == Parameters ===============================================================
%   x (vector) - nodes
%   a (vector) - left endpoints
%   b (vector) - right endpoints
% =============================================================================

V = transpose(fliplr(vander(x)));
p = (1:length(x))';
b = ((b(:)').^p - (a(:)').^p ) ./ p;
w = V \ b;
end

function z = nodes(q)
%QUADW Returns nodes for radau method
% == Parameters ===============================================================
%   q (vector) - right endpoints
% =============================================================================

if(q > 8) % use chebfun if q > 8
    if(isempty(which('radaupts')))
        error('q > 8 requires Chebfun (https://www.chebfun.org/)');
    end
    z = [-1; -flip(radaupts(q-1))];
    return;
end

nodes = {
    [1]'
    [-0.333333333333333,  1]'
    [-0.689897948556636,  0.289897948556636, 1]'
    [-0.822824080974592, -0.181066271118531,  0.575318923521694,  1]'
    [-0.885791607770965, -0.446313972723752,  0.167180864737834,  0.720480271312439,  1]'
    [-0.920380285897062, -0.603973164252784, -0.124050379505228,  0.390928546707272,  0.802929828402347,  1]'
    [-0.941367145680430, -0.703842800663031, -0.326030619437691,  0.117343037543100,  0.538467724060109,  0.853891342639482,  1]'
}; 
z = [-1; nodes{q-1}];

end