function [G1, G2, G3, G4] = G(rho, m, n, E, gamma)

    %Make variables that are going to "flux" in y direction


    p = (gamma - 1)*(E - 0.5*(m.^2 + n.^2)./rho);
    
    G1 = n;
    G2 = m.*n./rho;
    G3 = n.^2./rho + p;
    G4 = (E + p).*n./rho;
       
end