function pj( J )

Ncell = 3;

for r = 1 : length(J)
    sprintf( '%d      %.6e   %.6e   %.6e   %.6e   %.6e   %.6e   %.6e\n', r, J(r,1), J(r,2), J(r,3), J(r,4), J(r,5), J(r,6), J(r,7) )
end