function [sq1, sq2] = sharedVertsVectorized(t1, t2)
 %   [~,sq1]=intersect(t1, t2);
 %   [~,sq2]=intersect(t2, t1);
    
%    sq1 = find(ismember(t1,t2));
%    sq2 = find(ismember(t2,t1));

%    sq1 = fastintersect(t1, t2); 
%    sq2 = fastintersect(t2, t1);
   

%     ia = find(tf);
%     ib = ib(ia);
%     [c,iau] = unique([a(ia) b([])],order);
%     ia = ia(iau);
%     ib = ib(iau);

%     sq1 = [];
%     sq2 = [];
%     for i = 1:4
%         for j = 1:4
%             if t1(i) == t2(j)
%                 sq1 = [sq1 i];
%                 sq2 = [sq2 j];
%             end
%         end
%     end
    
    %% vectorized version starts here
    assert(size(t1,2)==4)
    assert(size(t2,2)==4)
    N = size(t1,1);
    equalities = t1==reshape(t2,N,1,4);
    [ii,jjkk] = find(reshape(equalities,[],16));
    [sortedii,perm] = sort(ii);
    
    % Expects inputs to have 3 matches per row. This must be guaranteed for any type of vectorization to work.
    assert(norm(reshape(sortedii,3,[]) - repmat(1:size(t1,1),3,1))==0);
        
    sortedjjkk = jjkk(perm);
    [jj, kk] = ind2sub([4 4],sortedjjkk);
    sq1 = reshape(jj,3,N)';
    sq2 = reshape(kk,3,N)';
    
    % randomized correctness verification
%     i = randi(N);
%     assert(all([t1(i,sq1(i,:))==t2(i,sq2(i,:))]));
end

% function idx = fastintersect(t1, t2)
%    P=zeros(min(min(t1),min(t2)), max(max(t1),max(t2)));
%    P(t1) = 1;
%    idx = find(logical(P(t2)));
% end
