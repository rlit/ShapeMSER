function [shape_, edges] = TessellateShape( shape, R )

f = fastmarchmex('init', int32(shape.TRIV-1), double(shape.X(:)), double(shape.Y(:)), double(shape.Z(:)));

D = zeros(size(R));

for k=0:max(R(:)),

    r = double(R==k);
    
    % Negative map
    source = repmat(Inf, [size(shape.X) 1]);
    source(r==1) = 0;
    dn = fastmarchmex('march', f, double(source));
    dn(dn>=9999999) = Inf;

    % Positive map
    source = repmat(Inf, [size(shape.X) 1]);
    source(r==0) = 0;
    dp = fastmarchmex('march', f, double(source));
    dp(dp>=9999999) = Inf;

    D(:,k+1) = -(dp-dn);
    
end

fastmarchmex('deinit', f);

[shape_, edges] = voronoi_tessellation(shape, D);