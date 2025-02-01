function y=mss_v2sT(x, n , m)
% function y=mss_v2sT(x, n, m)
%
% re-arrange the rows of x into a symmetrical tensor
% m is the order of the targeted tensor
% n is the dimension of the tensor


x=x(:);

if length(x)~=(nchoosek(m+n, n) - nchoosek(m+n-1, m-1)) , error('unable to convert: wrong dimension'); end

ytemp = zeros(n * ones(1,m));
%y = zeros(m, m, m, m);
% we may improve this part later


temp_matrix = ones(1,m);


a = 1:n;
if m == 2
    allind = allcomb(a, a);
elseif m ==4
    allind = allcomb(a, a, a, a);
elseif m == 6
    allind = allcomb(a, a, a, a, a, a);
elseif m == 8
    allind = allcomb(a, a, a, a, a, a, a, a);
elseif m == 10
    allind = allcomb(a, a, a, a, a, a, a, a, a, a);   
end 


size_allind = size(allind);
for i = 1:size_allind(1)
    indicator = 1;
    for j = 1:size(temp_matrix,1)
        if ismember(allind(i,:), perms(temp_matrix(j,:)), 'row')
            ytemp = new_tensor_vec(ytemp, allind(i,:), j);
            indicator = 0;
            break;
        end
    end
    if indicator == 1
        temp_matrix = [temp_matrix; allind(i,:)];
        ytemp = new_tensor_vec(ytemp, allind(i,:), size(temp_matrix,1));
    end
end


% temp_matrix = ones(1,4);
% for i1 = 1:m
%     for i2 = 1:m
%         for i3 = 1:m
%             for i4 = 1:m
%                 indicator = 1;
%                 for j = 1:size(temp_matrix,1)
%                     if ismember([i1, i2, i3, i4], perms(temp_matrix(j,:)), 'row')
%                         ytemp(i1,i2,i3,i4) = j;
%                         %y(i1, i2, i3, i4) = x(j);
%                         indicator = 0;
%                         break;
%                     end
%                 end
%                 if indicator == 1
%                     temp_matrix = [temp_matrix; [i1, i2, i3, i4]];
%                     ytemp(i1,i2,i3, i4) = size(temp_matrix,1); 
%                     %y(i1, i2, i3, i4) = x(j);
%                 end
%             end
%         end
%     end
% end

y = ytemp;