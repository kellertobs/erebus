% SecondInv    EREBUS subroutine to calculate tensor magnitude
%
% [aII]  =  SecondInv(a)
%
%   Function calculates second invariant tensor magnitude for stresses or
%   strain-rates.
%
%   modified  20170427  Tobias Keller
%   modified  20200515  Tobias Keller


function  [AII]  =  SecondInv(A)

AII  =  sqrt(1/2.*(A(:,1).^2 + A(:,2).^2 + 2.*A(:,3).^2));

end