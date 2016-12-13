function print_unionfind(q, N)
    u = [];
    p = unionfind('find', q, 1);
    for l=1:N,
        u(l) = unionfind('find', q, l);
    end
    v = unique(u);
    for l=1:length(v),
        fprintf(1,'{ ');
        fprintf(1,'%d ', find(u==v(l)));
        fprintf(1,'} ');
    end
