function y = f(x)
    sum=0;
    len=length(x);
    for i= 1:len
        sum=sum+i*x(i)^2;
    end
    y=sum;
end

