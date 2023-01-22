function [gradfx] = findiff_grad(f, x, h, type)

gradfx = zeros(size(x));

fx = f(x);

switch type
    case 'fw'
        for i=1:length(x)
            xh = x;
            xh(i) = xh(i) + h;
            gradfx(i) = (f(xh) - fx)/ h;
        end
    case 'c'
        for i=1:length(x)
            xh_plus = x;
            xh_minus = x;
            xh_plus(i) = xh_plus(i) + h;
            xh_minus(i) = xh_minus(i) - h;
            gradfx(i) = (f(xh_plus) - f(xh_minus))/(2 * h);
        end
    otherwise
        for i=1:length(x)
            xh = x;
            xh(i) = xh(i) + h;
            gradfx(i) = (f(xh) - fx)/h;
        end
end
end