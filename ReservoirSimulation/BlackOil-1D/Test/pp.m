function pp( U )

if length(U) == 3
    sprintf( '%.6e   %.6e   %.6e', U(1), U(2), U(3) )
else
    sprintf( '%.6e   %.6e   %.6e   %.6e', U(1), U(2), U(3), U(4) )
end