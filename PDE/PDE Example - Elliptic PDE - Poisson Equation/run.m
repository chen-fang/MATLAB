nx = 32; ny = 32;
dx = 1 / (nx-1);
dy = 1 / (ny-1);

% diagonal entries
diag_block = eye(ny-1) * (-2/dx^2-2/dy^2);
% upper diagonal entries
diag_block = diag_block + diag( ones(ny-2,1)/dy^2, 1 );
% lower diagonal entries
diag_block = diag_block + diag( ones(ny-2,1)/dy^2, -1 );

% construct diagonal block by repeating diag_block
Matrix = kron( eye(nx-1), diag_block );
% construct upper diagonal block
Matrix = Matrix + diag( ones((nx-2)*(ny-1), 1 ), ny-1 ) * 1/dx^2;
% construct lower diagonal block
Matrix = Matrix + diag( ones((nx-2)*(ny-1), 1 ), 1-ny ) * 1/dx^2;

% check sparsity pattern of matrix
figure;
spy( Matrix );

% meshgrid
x = [1:nx-1] * dx;
y = [1:ny-1] * dy;
[Y,X] = meshgrid( x, y );

% RHS
F = sin(X) .* cos(Y);
f = reshape( F', (nx-1)*(ny-1), 1 );

figure;
surf( X, Y, F );
xlabel('X');
ylabel('Y');

% solve
u = Matrix \ f;
U = reshape( u, ny-1, nx-1 )';

figure;
surf( X, Y, U );
xlabel('X');
ylabel('Y');