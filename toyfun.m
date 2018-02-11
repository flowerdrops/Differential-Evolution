function res=toyfun(ind)
    load('d_eku.mat');
    res=sum(sum(abs(d_eku-13*ind)));
end