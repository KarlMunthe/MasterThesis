function [F1, F2, F3, F4] = F(rho, m, n, E, gamma)

    %Make variables that are going to "flux" in x direction

    p = (gamma - 1)*(E - 0.5*(m.^2 + n.^2)./rho);

    F1 = m; 
    F2 = m.^2./rho + p;
    F3 = m.*n./rho;
    F4 = (E + p).*m./rho;
       
end