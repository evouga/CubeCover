function e = levicivita(I,Dim)
% levicivita
%   e = levicivita(I,Dim)
%
%   Returns the levi civiata of the input indices. When Dim == 1, each
%   column is seen as a vector of indices. When Dim == 2, each row is seen
%   as a vector of indices. Dim defaults to 1 when unspecified.
%
%   levi civiata is defined as follows
%       | -1 if indices are an odd permutation. i.e. [1,3,2], [3,2,1],...
%   e = | +1 if indices are an even permutation. i.e [1,2,3], [3,1,2],...
%       | 0  if any of the indices are equal. i.e. [1,1,1], [3,2,2],...
%
%   Example
%       [I{1:3}] = ndgrid(1:3,1:3,1:3);     % All 3D subscripts
%       I = [I{1}(:) I{2}(:) I{3}(:)];      % 3D subscripts, i,j,k
%       e = levicivita(I,2)                 % Levi civiata for 3D
%       
%       levicivita([1 2 3 4].')             % Levi civiata for 4D
%   
%   See also
%       levicivita permutationparity
%

%% Author Information
% Pierce Brady
% Smart Systems Integration Group - SSIG
% Cork Institute of Technology, Ireland
%

%% Assign defaults
if nargin<2 || isempty(Dim), Dim = 1; end

%% Main body
e = (-1).^permutationparity(I,Dim);                 % Get parity
i = ~logical(prod(diff(sort(I,Dim),[],Dim),Dim));   % Find any vectors with duplicate indices
e(i) = 0;                                           % Replace duplicate indoces with 0;

end

function [p] = permutationparity(P,Dim)
% permutationparity
%   p = permutationparity(P,Dim)
%
%   Returns the parity of a permutation. A value of 0 indicates an even 
%   permutation, while a value of 1 indicates an odd permutation.
%
%   The parity of permutations with duplicate values cannot be determined
%   with this routine.
%
%   `Dim` specifies which dimension the permutations are located on. 
%   If Dim == 1 each column of the matrix `P` is treated as a permutation,
%   If Dim == 2 each row of the matrix `P` is treated as a permutation.
%   If Dim is empty or unspecified, if defaults to treating each column as 
%   a different permutation.
%   
%   Example
%       P = [1 2 3 5 4];            % An odd permutation
%       pP = permutationparity(P)   % Get parity
%       pP = 1
%
%       P = [1 2 3 5 4; 1 2 4 5 3];     
%       pP = permutationparity(P,2)
%       pP = 1 
%            0
%
%       P = [1 2 3; 3 1 2; 2 3 1];  % Rows are odd, columns are even
%       pP = permutationparity(P,1)
%       pP = 0 0 0
%
%       pP = permutationparity(P,2)
%       pP = 1
%            1
%            1
%
%   See also
%       permutationparity perms randperm
%

%% Reference
%   'Pretty Permutation Parity Proof'
%   Alan Macdonald
%   Department of Mathematics
%   Luther College,
%   Decorah, IA 52101, U.S.A.
%   macdonal@luther.edu
%   June 16, 2004
%   
%   Let p be a permutation of [1 2 ... n]. Let G(p) be the number of times 
%   in p that a number is greater than a number to its right. For example,
%   G([2 4 1 3]) = 3. Note that G[1 2 ... n] = 0. 
%   A transposition of numbers in adjacent positions changes G by ±1. Every 
%   transposition can be expressed as a product of an odd number of such 
%   transpositions. Therefore every transposition changes the parity of G. 
%   Thus the number of transpositions used to obtain p is always even or 
%   always odd, according as G(p) is even or odd.
%

%% Author Information
% Pierce Brady
% Smart Systems Integration Group - SSIG
% Cork Institute of Technology, Ireland.
%

%%
[nRows,nCols] = size(P);    % Size of input matrix

%% Assign defaults
if nargin<2 || isempty(Dim)
    if nRows == 1,
        Dim = 2;
    elseif nCols == 1, 
        Dim = 1; 
    else
        Dim = 1;
    end
end

%% Main body
p = 0;  % Initiate counter

% Loop through each column, checking if the numbers in the current column are
% greater than the numbers in the remaining columns. Sum the number of times
% this occurs. This will sum to an even number for an even permutation and
% an odd number for an odd permutation.
if Dim == 1
    for i = 1:(nRows)
        p = sum(repmat(P(i,:),nRows-i,1) > P((i+1):end,:),1) + p;
    end
elseif Dim == 2
    for i = 1:(nCols)
        p = sum(repmat(P(:,i),1,nCols-i) > P(:,(i+1):end),2) + p;
    end
end
p = mod(p,2);   % Return 0 for even numbers and 1 for odd numbers
end

