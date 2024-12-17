%Integrand 2*pi*f*r*dr for Hankel transform.
%
% f      Function values
% r      Radial positions
% w      Integrand 2*pi*f*r*dr
%
%Principle:    Each function value f(r) is attributed to
%              the ring surface between the neighbouring
%              mid-points.
%
%              The integration interval is [0,r(end)].
%

%     Marcel Leutenegger � June 2006
%
function w = frdr(f, r)
    % 确保 r 是列向量
    r = r(:);
    w = numel(r);
    
    % 根据 r 的数量计算不同情况
    switch w
        case 0
            % 如果 r 为空，则 w 也为空
            w = [];
        case 1
            % 对于单个半径位置，计算公式为： w = pi * r^2 * f
            w = pi * r^2 * f;
        case 2
            % 对于两个半径位置，使用涉及差异和和的公式
            r1 = r(1);
            r2 = r(2);
            term1 = pi / 4 * sqrt(0.5) * (r1 + r2)^2;
            term2 = 4 * r2^2 - (r1 + r2)^2;
            w = (term1 + term2) * f(:);
        otherwise
            % 对于多个半径位置，基于给定公式计算
            term1 = pi / 4 * sqrt(0.5) * (r(1) + r(2)).^2;
            term2 = (r(3:end) - r(1:end-2)) .* (r(1:end-2) + 2 * r(2:end-1) + r(3:end));
            term3 = 4 * r(end)^2 - (r(end-1) + r(end)).^2;
            w = pi / 4 * [term1; term2; term3] .* f(:);
    end
end