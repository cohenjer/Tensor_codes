% Implementation of the vec(x) operator
% obs: For a M x N matrix, the vec(x) operator
% stacks its N columns into a MN x 1 column vector
function a= vec(A);
a=A(:);