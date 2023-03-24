function H = buildH_mass(poses, factor, order)
%BUILDH_MASS function to create multiple trafo matrices at once
%    Inputs:
%       - poses: multiple poses as a list (nx6)
%       - factor: scaling factor for the translations (scalar, default: 1)
%       - order: rotation order for R matrix ('XXX')
%   Outputs:
%       - H: trafo matrices (4x4xn)
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    if nargin < 3
        order = 'XYZ';
    end
    if nargin < 2
        factor = 1;
    end

    M = size(poses, 1);
    H = zeros(4,4,M);
    
    co = reshape(cos(poses(:,4)), 1,1,M);
    so = reshape(sin(poses(:,4)), 1,1,M);
    cp = reshape(cos(poses(:,5)), 1,1,M);
    sp = reshape(sin(poses(:,5)), 1,1,M);
    ck = reshape(cos(poses(:,6)), 1,1,M);
    sk = reshape(sin(poses(:,6)), 1,1,M);
    
    if strcmp(order, 'XYZ')
        % rotation order (R(z) * R(y) * R(x))
        H(1:3,1:3,:) = [ck .* cp, ck .* so .* sp - co .* sk, sk .* so + ck .* co .* sp;
                        cp .* sk, ck .* co + sk .* so .* sp, co .* sk .* sp - ck .* so;
                        -sp,      cp .* so,                  co .* cp];
    elseif strcmp(order, 'YZX')
        % rotation order (R(x) * R(z) * R(y))
        H(1:3,1:3,:) = [ck .* cp,                  -sk,      ck .* sp;
                        so .* sp + co .* cp .* sk, ck .* co, co .* sk .* sp - cp .* so;
                        cp .* sk .* so - co .* sp, ck .* so, co .* cp + sk .* so .* sp];
    elseif strcmp(order, 'YXZ')
        % rotation order (R(z) * R(x) * R(y))
        H(1:3,1:3,:) = [ck .* cp - sk .* so .* sp, -co .* sk, ck .* sp + cp .* sk .* so;
                        cp .* sk + ck .* so .* sp,  ck .* co, sk .* sp - ck .* cp .* so;
                        -co .* sp,                 so,        co .* cp];
    end
    
    H(1:3,4,:) = permute(poses(:, 1:3) * factor, [2,3,1]);
    H(4,4,:) = ones(1,1,M);

end

