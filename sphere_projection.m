function [xhat]=sphere_projection(x,c,r)
    if norm(x-c)<=r
        xhat= x;
    else 
        xhat=c+r*((x-c)/norm(x-c));
    end
end