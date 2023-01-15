function xhat=box_projection(x,feasible_set)
    xhat=max(min(x,feasible_set(2)),feasible_set(1));
end