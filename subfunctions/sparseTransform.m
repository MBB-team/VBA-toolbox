function sx = sparseTransform(x,P)
sx = (x.^2).*(2./(1+exp(-P*x))-1);