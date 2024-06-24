function y = get_CYY(CYY, i_1, j_1, i_2, j_2, nxperm, nyperm)
q = i_1+(nxperm)*(j_1-1);
r = i_2+(nxperm)*(j_2-1);
y = CYY(q,r);
end