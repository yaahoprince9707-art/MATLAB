function Q = calc_Area(c,xst,xen)
    % calculates area of polynomial 
    % c: coefficient of polynomial descending order, xst: start and end

    c = c(:); % Make c a column vector
    c = fliplr(c'); % Flip c (Least squares give a0 + a1 x + a2 x^2 ... we need exactly opposite;
    cint = polyint(c); % integrate the polynomial, contains +C
    Q = polyval(cint,xst)-polyval(cint,xen); % substitute the limits

end