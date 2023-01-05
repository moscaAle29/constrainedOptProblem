function [gradf]=findiff_grad(f,x,h,type)
dim=length(x);
base=eye(dim);
gradf=zeros(dim,1);
switch type
    case 'fw'
        for i=1:dim
            gradf(i)=(f(x+h*base(:,i))-f(x))/h;
        end
    otherwise
        for i=1:length(x)
            gradf(i)=(f(x+h*base(:,i))-f(x-h*base(:,i)))/(2*h);
        end
end