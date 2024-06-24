function y = get_CYP(CYP, i_P, j_P, i, j, nx, ny, nxperm, nyperm)
q = i_P+nx*(j_P-1);
r = i+(nxperm)*(j-1);
y = CYP(q,r);
end