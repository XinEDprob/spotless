function Ynew = new_tensor_vec(Y, index, newval)
% change value an entry of tensor Y according to index vector index to
% newvalue

% index: index vector
% Y: tensor
% newvalue: a new value

if length(size(Y)) ~= length(index) , error('unable to convert: wrong dimension'); end

Ynew = Y;
m = length(index);

if m == 2
    Ynew(index(1), index(2)) = newval;
elseif m == 4
    Ynew(index(1), index(2), index(3), index(4)) = newval;
elseif m == 6
    Ynew(index(1), index(2), index(3), index(4), index(5), index(6)) = newval;
elseif m == 8
    Ynew(index(1), index(2), index(3), index(4), index(5), index(6), index(7), index(8)) = newval;
elseif m == 10
    Ynew(index(1), index(2), index(3), index(4), index(5), index(6), index(7), index(8), index(9), index(10)) = newval;
   
end
end