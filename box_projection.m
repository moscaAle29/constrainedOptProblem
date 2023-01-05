function xhat=box_projection(x,mins,maxs)
    xhat=max(min(x,maxs),mins);
end