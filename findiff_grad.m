function [gradfx] = findiff_grad(x, h, type)
    len=length(x);
    F=@(x) (1:len)'.*(x.^2);
    switch type
        case 'fw'
            gradfx = (F(x+h*ones(len,1)) - F(x))/ h;
        case 'c'
            gradfx = (F(x+h*ones(len,1)) - F(x-h*ones(len,1)))/(2 * h);
        otherwise
            gradfx = (F(x+h*ones(len,1)) - F(x))/ h;
    end
end