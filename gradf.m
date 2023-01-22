function gradf_value= gradf(x)
    len_x=length(x);
    gradf_value=zeros(len_x,1);
    for i=1:len_x
        gradf_value(i)=2*i*x(i);

    end
end

