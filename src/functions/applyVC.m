function vec_corr = applyVC(vec, omegas, apply, abs_line)
% APPLYVC applies the vertical correction for velodyne scanners
% expects:
% 	- vec: vector of z-coordinates [nx1]
%	- omegas: elevation angles corresponding to points [nx1]
%	- apply: apply (true) or remove (false) vertical correction
%	- abs_line (optional): ring number of points [nx1] (only for high performance necessary)
% returns:
%	- vec_corr: corrected vector of z-coordinates [nx1]
%
% Copyright (c) 2023 Dominik Ernst under MIT License

    % find corresponding scan line
    if nargin < 4
        abs_line = determineScanLine(omegas);
    end

    % vertical corrections - symmetrical
    vc = [0.7, 2.2, 3.7, 5.1, 6.6, 8.1, 9.7, 11.2]' * 1e-3;

    if apply
        factor = -1;
    else
        factor = 1;
    end
    vec_corr = vec + vc(abs_line) .* factor .* sign(omegas);
end

