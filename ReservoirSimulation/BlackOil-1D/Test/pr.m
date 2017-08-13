function pr( Res )

sprintf( '%s      %.6e   %.6e   %.6e\n', 'OIL', Res(1), Res(3), Res(5) )
sprintf( '%s      %.6e   %.6e   %.6e\n', 'WAT', Res(2), Res(4), Res(6) )

if length(Res) > 6
    sprintf( '%s      %.6e\n', 'INJ', Res(7) )
end