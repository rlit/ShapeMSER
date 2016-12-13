function F = SoftVQ(desc,vocab,sigma)

D = squared_dist(vocab, desc);
w = exp(-0.5*D/sigma^2);
ws = sum(w,1);             % L1-normalization
ws(ws <= 0) = 1;
F  = bsxfun(@rdivide, w, ws);