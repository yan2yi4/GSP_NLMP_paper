function [V] = FLOM( p, alpha,gam)
    C = (2^(p+1)*gamma((p+1)/2)*gamma(-p/alpha))/(alpha*(pi^0.5)*gamma(-p/2));
    if(p<0)
        C = (2^(p)*gamma((p+1)/2)*gamma(1-p/alpha))/((pi^0.5)*gamma(1-p/2));
    end
    V = C*gam^(p/alpha);
end
