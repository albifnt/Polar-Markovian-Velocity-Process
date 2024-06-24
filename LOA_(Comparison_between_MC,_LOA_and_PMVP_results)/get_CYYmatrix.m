function CYY = get_CYYmatrix(Grid, Y_var, Y_corr)

nxperm = Grid.nxperm;
nyperm = Grid.nyperm;
xx = Grid.xx_perm;
yy = Grid.yy_perm;

x_grid = Grid.xx(1,:);
y_grid = Grid.yy(:,1);

CYY=zeros(nxperm*nyperm,nxperm*nyperm);
for q=1:nxperm*nyperm
    for r=1:nxperm*nyperm
        [i_q,j_q]=get_ij(q,nxperm,nyperm);
        [i_r,j_r]=get_ij(r,nxperm,nyperm);
%         dx = (i_q-i_r)/(nxperm);
%         dy = (j_q-j_r)/(nyperm);
        
        dx = ( xx(1,i_q) - xx(1,i_r) );
        dy = ( yy(j_q,1) - yy(j_r,1) );
        
%         TO SLOW
%         [d,id1_x] = min( abs( x_grid - xx(1,i_q) ) );
%         [d,id1_y] = min( abs( y_grid - yy(j_q,1) ) );
%         
%         [d,id2_x] = min( abs( x_grid - xx(1,i_r) ) );
%         [d,id2_y] = min( abs( y_grid - yy(j_r,1) ) );
        
        id1_x = fix(i_q/2);
        if (id1_x == 0)
            id1_x = 1;
        end
        id1_y = fix(j_q/2);
        if (id1_y == 0)
           id1_y = 1;
        end
        
        id2_x = fix(i_r/2);
        if (id2_x == 0)
            id2_x = 1;
        end
        id2_y = fix(j_r/2);
        if (id2_y == 0)
           id2_y = 1;
        end
        
        
        CYY(q,r)= sqrt(Y_var(id1_y, id1_x))*sqrt(Y_var(id2_y, id2_x))*exp(-sqrt((dx/Y_corr)^2+(dy/Y_corr)^2)); %Y_var*exp(-sqrt((dx/Y_corr)^2+(dy/Y_corr)^2));
    end
end
end