function s = stencil_n(n)
    all = all_stencils();
    s = [];
    if n == 2
        s = all{1};
    end
    if n == 3.5
        s = all{2};
    end
    if n == 4
        s = all{3};
    end
    if n == 4.5
        s = all{4};
    end
    if n == 5.5
        s = all{5};
    end
end
    