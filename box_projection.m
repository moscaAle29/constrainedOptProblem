function xhat=box_projection(x,feasible_set)
    mins = feasible_set(:,1);
    maxs = feasible_set(:,2);
    xhat = max(min(x, maxs), mins);
end