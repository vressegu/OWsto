init;

beta = 1:0.1:6;
bound = ((beta.^2-1)/8).^(1./(beta-3));

plot(beta,bound)
