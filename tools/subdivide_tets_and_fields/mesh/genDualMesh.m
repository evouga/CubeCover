% 
% This code is part of a utility called mint, 
% 
% which aspires to be useful for designing frame fields.  
% 
% See the readme.md
%
% written by josh vekhter in 2020.  Some Permissive License, tbd.  
%
%
%%% TODO: Come back and fix this, no reason to do it the slow way.

function [internalEdges, boundaryTris2] = genDualMesh(tets, verts)
data= processXTtetra_paul(tets,verts);
    
ntets = size(tets,1);
                  
internalEdges = zeros(sum(~data.isBoundaryTriangle),4);
internalEdges(:,[1 3]) = data.triangles2tets(~data.isBoundaryTriangle,:);

t1 = tets(internalEdges(:,1),:);
t2 = tets(internalEdges(:,3),:);
[sq1, sq2] = sharedVertsVectorized(t1, t2);
internalEdges(:,2) = tetTriToIdx_vectorized(sq1);
internalEdges(:,4) = tetTriToIdx_vectorized(sq2);

% can make this more efficient...
% for ij = 1:size(internalEdges,1)
%     i = internalEdges(ij,1);
%     j = internalEdges(ij,3);
%         
%     t1 = tets(i, :);
%     t2 = tets(j, :); 
%     [sq1, sq2] = sharedVerts(t1, t2);
%     
%     if size(sq1,2) == 3
%         % allFaceEdgesSq = [allFaceEdgesSq; [sq1 sq2] ];
%         internalEdges(ij,[2 4]) = [tetTriToIdx(sq1) tetTriToIdx(sq2)];
% 
%         v1 = verts(t1(sq1(1)),:);
%         v2 = verts(t1(sq1(2)),:);
%         v3 = verts(t1(sq1(3)),:);
% 
%         sharedNormals = [sharedNormals; cross( v1-v2, v1-v3 ) ];
%     end 
% end

%%%
% Build up boundary data
%%%
% boundaryTris = ones(ntets, 4);
% for e = 1:size(internalEdges,1)
%     edge = internalEdges(e,:);
%     boundaryTris(edge(1), edge(2)) = 0;
%     boundaryTris(edge(3), edge(4)) = 0;
% end
boundaryTris2 = ones(ntets, 4);
boundaryTris2(sub2ind(size(boundaryTris2), internalEdges(:,1), internalEdges(:,2))) = 0;
boundaryTris2(sub2ind(size(boundaryTris2), internalEdges(:,3), internalEdges(:,4))) = 0;

% restore permutation to match previous implementation
[~, perm] = sortrows(internalEdges(:,[1 3]));
internalEdges = internalEdges(perm,:);

end