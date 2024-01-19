% a snapshot of paul's vectorized tet mesh loading code as of 1/7/2021
function data = processXTtetra_paul(T,X,lite,force,anglethresh)

    if ~exist('anglethresh','var') || anglethresh<0
        anglethresh = 25;
    end

    if(nargin == 2)
        lite = true;
        force = false;
    elseif(nargin == 3)
        force = false;
    end
    
    % lite means don't load all properties. force means don't set to lite
    % even if size is far beyond recommended. will take a long time to run.
    if ((size(T,1) > 330000 || size(X,1) > 61000) && ~force)
        'SETTING TO LITE. TOO MANY VERTICES/TETS.'
        lite = 1;
    end

    data.tetrahedra = T;
    data.vertices = X;
    data.numVertices = size(X,1);
    data.numTetrahedra = size(T,1);

    % 123 214 134 324
    tri = [T(:,1) T(:,2) T(:,3) ; T(:,2) T(:,1) T(:,4) ;
           T(:,1) T(:,3) T(:,4) ; T(:,3) T(:,2) T(:,4) ];
%     opposingV = [T(:,4); T(:,3) ;
%        T(:,2); T(:,1)];

    [~,IA,IC] = unique(sort(tri,2),'rows');
    data.triangles = tri(IA,:);
    
    triind = [(1:data.numTetrahedra)';(1:data.numTetrahedra)';(1:data.numTetrahedra)';(1:data.numTetrahedra)'];
%     data.triangles2tets = triind(IA);
    
    data.numTriangles = size(data.triangles,1);
    data.tetsToTriangles = reshape(IC,data.numTetrahedra,4);

%     triangles2opposingV = opposingV(IA);
%     data.tetsToTrianglesToOpposingV = triangles2opposingV(data.tetsToTriangles);
    
    counts = accumarray(IC,ones(size(IC)),[data.numTriangles 1]);
    data.isBoundaryTriangle = counts == 1;
    bt = data.triangles(find(data.isBoundaryTriangle),:);
    data.isBoundaryVertex = zeros(1,data.numVertices)==1; data.isBoundaryVertex(unique(sort(bt(:))))=true;
    data.boundaryTriangles = find(data.isBoundaryTriangle);

    tt = data.triangles;
    data.triangleNormals = -cross(X(tt(:,2),:)-X(tt(:,1),:),X(tt(:,3),:)-X(tt(:,1),:));
    data.triangleAreas = sqrt(sum(data.triangleNormals.^2,2))/2;
    data.triangleNormals = bsxfun(@rdivide,data.triangleNormals,2*data.triangleAreas);

    n = data.numTetrahedra;

    data.tetBarycenters = .25*(X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:)+X(T(:,4),:));
    data.triangleBarycenters = (X(data.triangles(:,1),:)+X(data.triangles(:,2),:)+X(data.triangles(:,3),:))/3;

    E1 = X(T(:,2),:) - X(T(:,1),:);
    E2 = X(T(:,3),:) - X(T(:,1),:);
    E3 = X(T(:,4),:) - X(T(:,1),:);
    data.tetVolumes = abs(dot(cross(E1,E2),E3,2)/6);

    %%
    % E is pairs of vertices for each triangle
    E = [data.triangles(:,1) data.triangles(:,2) ; data.triangles(:,2) data.triangles(:,3) ; data.triangles(:,3) data.triangles(:,1)];
    [~,IA,IC] = unique(sort(E,2),'rows');
    data.edges = E(IA,:);
    flipped = find(data.edges(:,1)>data.edges(:,2)); data.edges(flipped,:)=[data.edges(flipped,2) data.edges(flipped,1)];
    data.numEdges = size(data.edges,1);
    data.trianglesToEdges = reshape(IC,data.numTriangles,3);

    % get tetsToEdges
    E = [data.tetrahedra(:,1) data.tetrahedra(:,2) ; ...
        data.tetrahedra(:,1) data.tetrahedra(:,3) ; ...
        data.tetrahedra(:,1) data.tetrahedra(:,4) ; ...
        data.tetrahedra(:,2) data.tetrahedra(:,3) ; ...
        data.tetrahedra(:,2) data.tetrahedra(:,4) ; ...
        data.tetrahedra(:,3) data.tetrahedra(:,4)];

    [~,IA,IC] = unique(sort(E,2),'rows');
    TetEdges = E(IA,:);
    data.tetsToEdges = reshape(IC,data.numTetrahedra,6);

    % check that tets to edges is right. can make check random for speed. Or
    % prove but im not sure...
    tetsToEdgesWithDup = [reshape(data.trianglesToEdges(data.tetsToTriangles,1),size(data.tetsToTriangles)) ...
        reshape(data.trianglesToEdges(data.tetsToTriangles,2),size(data.tetsToTriangles)) ...
        reshape(data.trianglesToEdges(data.tetsToTriangles,3),size(data.tetsToTriangles))];
    
    %for j = 1:data.numTetrahedra
    %    assert(sum(sort(unique(tetsToEdgesWithDup(j,:))) == sort(data.tetsToEdges(j,:)))==6);
    %end

    data.isBoundaryEdge = zeros(data.numEdges,1);
    for i=1:3
        data.isBoundaryEdge(data.trianglesToEdges(data.boundaryTriangles,i)) = 1;
    end
    data.boundaryEdges = find(data.isBoundaryEdge);
    adj = sparse(repmat((1:n)',4,1),data.tetsToTriangles(:),ones(4*n,1)); % tets x tris
    adj_trilabeled = sparse(repmat((1:n)',4,1),data.tetsToTriangles(:),data.tetsToTriangles(:)); % tets x tris
    adj_tetlabeled = sparse(repmat((1:n)',4,1),data.tetsToTriangles(:),repmat([1:data.numTetrahedra]',4,1)); % tets x tris
    
    % compute tri x tri 2 tet
    trixtri2tet = (adj'*adj_tetlabeled);
    trixtri2tet(1:data.numTriangles+1:end)=0;
    data.trixtri2tet = trixtri2tet;
    %{
    for i=1:100000
        rtet = randsample(data.numTetrahedra,1);
        rtris = randsample(data.tetsToTriangles(rtet,:),2);
        assert(trixtri2tet(rtris(1),rtris(2))==rtet);
    end
    %}
    
    % compute triangles to tets. one tet will be 0 for boundary triangles
    bind = sum(adj);
    [bind1,perm]=sort(bind);
    adjp = adj(:,perm);
    boundarypart = adjp(:,1:(find(bind1==2,1)-1));
    intpart = adjp(:,find(bind1==2,1):end);
    [first_boundaryTri2Tet, ~]=find(boundarypart); 
    first_boundaryTri2Tet_padded = [first_boundaryTri2Tet zeros(numel(first_boundaryTri2Tet),1)];
    [ii,~] = find(intpart); first_intTri2Tet = reshape(ii,2,[])';
    first_tri2tet = [first_boundaryTri2Tet_padded; first_intTri2Tet];
    triangles2tets = first_tri2tet*0; 
    %josh comment 
    triangles2tets(perm,:) = first_tri2tet;
    data.triangles2tets = triangles2tets; % [tet1 0;...] for boundary tris. [tet1 tet2] for interior tris.
    
    tetxtet_trilabeled = adj_trilabeled*adj';
    tetxtet_trilabeled(1:data.numTetrahedra+1:end)=0;
    [ii,jj,kk]=find(tetxtet_trilabeled);
    data.tetXtri2tet = sparse(ii,kk,jj,data.numTetrahedra,data.numTriangles);
    data.tetXtri2tet(:,data.isBoundaryTriangle)=0;
    
    % verify data.tetXtri2tet is right.
    %{
    for i=1:100000
        rtet = randsample(data.numTetrahedra,1);
        rtri = randsample(data.tetsToTriangles(rtet,:),1);
        otet = data.tetXtri2tet(rtet,rtri);
        if data.isBoundaryTriangle(rtri)
            assert(data.tetXtri2tet(rtet,rtri)==0);
        else
            assert(data.tetXtri2tet(rtet,rtri)~=0);
            assert(rtri==intersect(data.tetsToTriangles(otet,:), data.tetsToTriangles(rtet,:)));
        end
    end
    %}
    
    % get trianglesToTets
    if(~lite)
        numberedAdj = adj.*[1:size(adj,1)]';
        collapsedAdj = numberedAdj(find(numberedAdj~=0));
        trianglesToTets = mat2cell(full(collapsedAdj), full(sum(adj)));
        data.trianglesToTets = trianglesToTets;
        
        %% compute nonboundaryTriToTets
        nonBoundaryTrianglesToTets = cell2mat(data.trianglesToTets(find(~data.isBoundaryTriangle)));
        nonBoundaryTrianglesToTets = reshape(nonBoundaryTrianglesToTets,2,numel(nonBoundaryTrianglesToTets)/2)';
        data.nonBoundaryTrianglesToTets =nonBoundaryTrianglesToTets ;
        
        dualGraphAdj = sparse(nonBoundaryTrianglesToTets(:,1),nonBoundaryTrianglesToTets(:,2),ones(size(nonBoundaryTrianglesToTets(:,1),1),1),data.numTetrahedra,data.numTetrahedra);
        dualGraphAdj = dualGraphAdj + dualGraphAdj' ~= 0;
        data.dualGraphAdjacency = dualGraphAdj;
        data.dualGraph = graph(data.dualGraphAdjacency);
    end
    
    % verts to tets
    tetsToVertsIndicator = sparse(repmat((1:data.numTetrahedra)',4,1),data.tetrahedra(:),ones(4*data.numTetrahedra,1),data.numTetrahedra, data.numVertices);
    data.tetsToVertsIndicator = tetsToVertsIndicator;
    
    % get edgesToTets
    tetsToEdgesIndicator = sparse(repmat((1:data.numTetrahedra)',6,1),data.tetsToEdges(:),ones(6*data.numTetrahedra,1),data.numTetrahedra, data.numEdges);
    data.tetsToEdgesIndicator = tetsToEdgesIndicator;
    if(~lite)
        
        numberedTetsToEdgesIndicator = tetsToEdgesIndicator.*[1:size(tetsToEdgesIndicator,1)]';
        collapsednumberedTetsToEdgesIndicator = numberedTetsToEdgesIndicator(find(numberedTetsToEdgesIndicator~=0));

        edgesToTets = mat2cell(full(collapsednumberedTetsToEdgesIndicator), full(sum(tetsToEdgesIndicator)));
        data.edgesToTets = edgesToTets;
    end
    
    % get edgesToTrianglesUnoriented
    trianglesToEdgesIndicator = sparse(repmat((1:data.numTriangles)',3,1),data.trianglesToEdges(:),repmat((1:data.numTriangles)',3,1));
    data.edgesToTrianglesIndicator = trianglesToEdgesIndicator';
    if(~lite)
        
        collapsedTrianglesToEdgesIndicator = trianglesToEdgesIndicator(find(trianglesToEdgesIndicator~=0));
    
        edgesToTriangles = mat2cell(full(collapsedTrianglesToEdgesIndicator), full(sum(trianglesToEdgesIndicator~=0)));
        data.edgesToTrianglesUnoriented = cellfun(@transpose,edgesToTriangles,'un',0);
    end
    
    %nonboundary triangle indices
    triangleSubindices = find(~data.isBoundaryTriangle);
    adj = adj(:,triangleSubindices);
    [Itet,Jtri,~] = find(adj);
    Jtri = triangleSubindices(Jtri);
    triTet = sortrows([Jtri Itet]);

    assert(isequal(triTet(1:2:end,1),triTet(2:2:end,1)));

    % get boundaryTets
    data.isBoundaryTet = sum(data.isBoundaryTriangle(data.tetsToTriangles),2)~=0;

    %tet1, edge, tet2
    entries = [];
    for j=1:3
        e = data.trianglesToEdges(triTet(1:2:end,1),j);
        possibleEntries = [];
        possibleEntries = [possibleEntries; triTet(1:2:end,2) e triTet(2:2:end,2) ];
        possibleEntries = [possibleEntries; triTet(2:2:end,2) e triTet(1:2:end,2) ];
        %idx = ~data.isBoundaryEdge(e);
        %possibleEntries = possibleEntries([idx;idx],:);
        entries = [entries;possibleEntries];
    end
    mtx = sortrows(entries);
    
    % mtx = sortrows([plusTet commonEdge minusTet]);
    nextTet1 = sparse(mtx(1:2:end,1),mtx(1:2:end,2),mtx(1:2:end,3),data.numTetrahedra,data.numEdges);
    nextTet2 = sparse(mtx(2:2:end,1),mtx(2:2:end,2),mtx(2:2:end,3),data.numTetrahedra,data.numEdges);

    % tets cycles
    edgeCycles = cell(data.numEdges,1);
    ne = data.numEdges;
    isBoundaryEdge = data.isBoundaryEdge;
    if(~lite)
%         tic
%         for i=1:data.numEdges
%             if mod(i,1000) == 1
%                 fprintf('Edge %d of %d...\n',i,ne);
%             end
% 
%             if isBoundaryEdge(i)
%                continue
%             end
% 
%             rtri = data.edgesToTrianglesUnoriented{i}(1);
%             adjTets = data.edgesToTets{i};
%             tets = data.trianglesToTets{rtri};
%             [~,ind] = ismember(tets,adjTets);
%             subG = subgraph(data.dualGraph, adjTets);
% %             subA = data.dualGraphAdjacency(adjTets,adjTets);
% %             subG = graph(subA); subG = subG.rmedge(ind(1),ind(2));
%             subG = subG.rmedge(ind(1),ind(2));
%             cycle = shortestpath(subG,ind(1),ind(2));
%             cycle = adjTets(cycle);
%             
%             assert(length(cycle) > 2);
%             edgeCycles{i} = cycle;
%         end
%         toc
%         tic
        for i=1:data.numEdges
            if mod(i,1000) == 1
                fprintf('Edge %d of %d...\n',i,ne);
            end

            if isBoundaryEdge(i)
               continue
            end

            cycle = find(nextTet1(:,i),1);
            cycle = [cycle nextTet1(cycle,i)];
            while (cycle(end) ~= cycle(1))
                if nextTet1(cycle(end),i) == cycle(end-1)
                    cycle = [cycle nextTet2(cycle(end),i)];
                else
                    cycle = [cycle nextTet1(cycle(end),i)];
                end
            end
            assert(length(cycle) > 2);
            edgeCycles{i} = cycle;
        end
%         toc
        

        for i=1:data.numEdges
            if mod(i,1000) == 1
                fprintf('Edge %d of %d...\n',i,ne);
            end

            if ~isBoundaryEdge(i)
               continue
            end

            unorderedTets = data.edgesToTets{i};
            if(numel(unorderedTets)==1)
                cycle = unorderedTets;
                edgeCycles{i} = cycle;
                continue;
            end
            % since tets can have multiple faces on boundary, generally, a boundary edge could have
            % more than 2 'boundary' tets.
            boundaryTets = unorderedTets(find(data.isBoundaryTet(unorderedTets)));

            % get the boundary tets that start and end the cycle for edge i.
            linedUpEdgeMatch = reshape(data.trianglesToEdges(data.tetsToTriangles(boundaryTets,:)',:)',12,numel(boundaryTets))' == i;
            TriangleWithBoundaryEdge = [sum(linedUpEdgeMatch(:,1:3),2) sum(linedUpEdgeMatch(:,4:6),2) sum(linedUpEdgeMatch(:,7:9),2) sum(linedUpEdgeMatch(:,10:12),2)] > 0;
            BoundaryTriangles = data.isBoundaryTriangle(data.tetsToTriangles(boundaryTets,:));
            TrueBoundaryTets = boundaryTets(find(sum(TriangleWithBoundaryEdge & BoundaryTriangles, 2)>0));
            assert(numel(TrueBoundaryTets)==2); % these boundary tets are the start and end of the cycle for edge i

            cycle = TrueBoundaryTets(1);

            while numel(cycle) ~= numel(unorderedTets)
                currentTet = cycle(end);
                triangs = data.tetsToTriangles(currentTet,:);
                inttri = triangs(find(~data.isBoundaryTriangle(triangs)));
                exttri = triangs(find(data.isBoundaryTriangle(triangs)));
                neighborTets = unique([reshape([data.trianglesToTets{inttri,:}],1,2*numel(inttri)) data.trianglesToTets{exttri}]);
                neighborTets(find(neighborTets == currentTet)) = [];
                viableNextTets = intersect(unorderedTets,neighborTets); % one or two tets
                assert(numel(viableNextTets)<=2);
                if(sum(find(cycle==viableNextTets(1)))>0)
                    cycle = [cycle viableNextTets(2)];
                else
                    cycle = [cycle viableNextTets(1)];
                end
                edgeCycles{i} = cycle;
            end
            assert(data.isBoundaryTet(cycle(end)));

        end

        data.edgeCycles = edgeCycles;

        isOrientedEdge = zeros(data.numEdges, 1);
        isOrientedTriangle = zeros(data.numTriangles, 1);
        startTriangle = find(~data.isBoundaryTriangle,1);
        startEdge = data.trianglesToEdges(startTriangle,1);
        isOrientedEdge(startEdge) = true;
        trianglesToTraverse = data.edgesToTrianglesUnoriented{startEdge};
        while numel(trianglesToTraverse)~=0
            currentTriangle = trianglesToTraverse(1);

            trianglesToTraverse = trianglesToTraverse(2:end);
            if(isOrientedTriangle(currentTriangle))
                continue;
            end

            verts = data.triangles(currentTriangle, :);
            edges = data.trianglesToEdges(currentTriangle, :);
            edgeverts = data.edges(edges,:);

            %if data.edgecycles{k} has only one tet, orientation is not well
            %defined.
            badOrientation = [numel(data.edgeCycles{edges(1)}) numel(data.edgeCycles{edges(2)}) numel(data.edgeCycles{edges(3)})]==1;

            orientedEdgeInd = find(isOrientedEdge(edges)' & ~badOrientation,1);
            assert(numel(orientedEdgeInd)==1);
            remainingEdges = edges(find([1:3]~=orientedEdgeInd));
            tets = data.trianglesToTets{currentTriangle};

            edgeOrientations = [findOrderOfAInB(edgeverts(1,:), verts) findOrderOfAInB(edgeverts(2,:), verts) findOrderOfAInB(edgeverts(3,:), verts)];
            forwardEdgeOrientation = edgeOrientations(orientedEdgeInd);

            tetOrientations = [findOrderOfAInB(tets, data.edgeCycles{edges(1)}) findOrderOfAInB(tets, data.edgeCycles{edges(2)}) findOrderOfAInB(tets, data.edgeCycles{edges(3)})];

            forwardTetOrientation = tetOrientations(orientedEdgeInd);
            ind1 = find(edges==remainingEdges(1));
            ind2 = find(edges==remainingEdges(2));

            if(~isOrientedEdge(edges(ind1)) && ...
                (   xor(tetOrientations(ind1), edgeOrientations(ind1)) ~= xor(forwardTetOrientation, forwardEdgeOrientation) ...
                    || numel(data.edgeCycles{edges(ind1)}) == 1 ))
                data.edgeCycles{edges(ind1)} = fliplr(data.edgeCycles{edges(ind1)});
            end

            if(~isOrientedEdge(edges(ind2)) && ...
                (xor(tetOrientations(ind2), edgeOrientations(ind2)) ~= xor(forwardTetOrientation, forwardEdgeOrientation) ...
                || numel(data.edgeCycles{edges(ind2)}) == 1))
                data.edgeCycles{edges(ind2)} = fliplr(data.edgeCycles{edges(ind2)});
            end

            % double check orientation is right now.
            tetOrientations = [findOrderOfAInB(tets, data.edgeCycles{edges(1)}) findOrderOfAInB(tets, data.edgeCycles{edges(2)}) findOrderOfAInB(tets, data.edgeCycles{edges(3)})];
            assert(badOrientation(ind1) || ~(xor(tetOrientations(ind1), edgeOrientations(ind1)) ~= xor(forwardTetOrientation, forwardEdgeOrientation)) || numel(data.edgeCycles{edges(ind1)}) == 1);
            assert(badOrientation(ind2) || ~(xor(tetOrientations(ind2), edgeOrientations(ind2)) ~= xor(forwardTetOrientation, forwardEdgeOrientation)) || numel(data.edgeCycles{edges(ind2)}) == 1);

            isOrientedEdge(edges(ind2))=true;
            isOrientedEdge(edges(ind1))=true;

            isOrientedTriangle(currentTriangle)=true;

            % badorientation edges don't need to have their triangles added to
            % the traversal list. They will be added when one of their non-bad edges is
            % oriented. If they are never added, that means all 3 of their
            % edges were badorienation, meaning had only 1 adjacent tet. This
            % corresponds to only one scenario: the tet is a single alone tet.
            % we're not interested in such a simple case.
            moreToTraverse = unique([data.edgesToTrianglesUnoriented{edges(find(~badOrientation))}]);
            moreToTraverse = moreToTraverse(~isOrientedTriangle(moreToTraverse));
            trianglesToTraverse = [trianglesToTraverse moreToTraverse];
        end
    
        %% verify that cycles are correctly oriented for interior edges
        for triIter = [1:data.numTriangles]
            edges = data.trianglesToEdges(triIter, :);
            if data.isBoundaryTriangle(triIter) || sum(data.isBoundaryEdge(edges))>0
                continue
            end

            edges = data.trianglesToEdges(triIter, :);
            edgeverts = data.edges(edges,:);
            verts = data.triangles(triIter, :);

            orient1 = [findOrderOfAInB(edgeverts(1,:), verts) findOrderOfAInB(edgeverts(2,:), verts) findOrderOfAInB(edgeverts(3,:), verts)];

            tets = data.trianglesToTets{triIter};

            tets1 = data.edgeCycles{edges(1)};
            tets2 = data.edgeCycles{edges(2)};
            tets3 = data.edgeCycles{edges(3)};

            orient2 = [findOrderOfAInB(tets, tets1) findOrderOfAInB(tets, tets2) findOrderOfAInB(tets, tets3)];

            assert(sum(orient1 == orient2)==3 | sum(orient1 == orient2)==0);
        end
    
        
        
        for i = 1:numel(data.edgeCycles)
            isbe = data.isBoundaryEdge(i);
            tetCycle = data.edgeCycles{i};
            neighboringTris = data.edgesToTrianglesUnoriented{i};
            ntris = numel(tetCycle)+1;
            tets2tris = data.tetsToTriangles(tetCycle,:);
            
            if(isbe)
                if(numel(tetCycle)==1)
                    % trying to carry orientation from tets. 1 tet alone
                    % has no orientation. gonna be lazy and hope users
                    % subdivide.
                    display('trouble! cannot have tet with 2 boundary faces!');
                end
                
                startTris = intersect(tets2tris(1,:),neighboringTris);
                startTri = startTris(find(data.isBoundaryTriangle(startTris)));
                triCycle = [startTri ARemoveB(startTris,startTri)];
                
                for j = 2:size(tets2tris,1)
                    tris = intersect(tets2tris(j,:),neighboringTris);
                    nextTri = ARemoveB(tris, triCycle(j));
                    triCycle = [triCycle nextTri];
                end
            else
                startTris = intersect(tets2tris(1,:),neighboringTris);
                nextTri = intersect(tets2tris(1,:),tets2tris(2,:));
                triCycle = [ARemoveB(startTris,nextTri) nextTri];
                for j = 2:(size(tets2tris,1)-1)
                    tris = intersect(tets2tris(j,:),neighboringTris);
                    nextTri = ARemoveB(tris, triCycle(j));
                    triCycle = [triCycle nextTri];
                end
            end
            data.edgeTriCycles{i}=triCycle;
        end
        
        % verify edge TriCycles
%         for i = 1:numel(data.edgeTriCycles)
%             triCycle = data.edgeTriCycles{i};
%             vs = data.vertices(data.edges(i,:),:); 
%             %plot3(vs(:,1),vs(:,2),vs(:,3),'r'); 
%             hold on; axis equal;
%             triloop = data.triangleBarycenters(triCycle,:);
%             %plot3(triloop(:,1),triloop(:,2),triloop(:,3),'r');
%             triloop = .2*triloop + .9*sum(vs)/2;
%             triverts = data.vertices(data.triangles(triCycle,:),:);
%             scatter3(triverts(:,1),triverts(:,2),triverts(:,3));
%             plot3(triloop(:,1),triloop(:,2),triloop(:,3),'b');
%             %hold off;
%         end
        
    end
    
    %% separate nonboundary and boundary edges
    data.NonBoundaryEdges = data.edges(find(~data.isBoundaryEdge),:);
    data.BoundaryEdges = data.edges(find(data.isBoundaryEdge),:);
    
    
    %% compute primal spanning tree of the volume. spans vertices
    % data.PrimalVolumeVertexSpanningTree = PrimalVolumeVertexSpanningTree(data.edges);
    
    % compute dual spanning tree of the volume. spans tets
    % data.DualVolumeVertexSpanningTree = DualVolumeVertexSpanningTree(data);
    
    % compute primal spanning tree of volume WITHOUT BOUNDARY EDGES
    %Inds=find(~data.isBoundaryEdge); Inds = Inds(PrimalVolumeVertexSpanningTree(data.NonBoundaryEdges));
    %data.BoundaryLessPrimalSpanningTreeRelToEdges = Inds;
    
    assert(all(data.edges(:,1)<=data.edges(:,2)));
    
    % compute incidence matrix
    data.primalIncidenceMatrix = sparse([1:data.numEdges 1:data.numEdges]', [data.edges(:,1);data.edges(:,2)], [ones(data.numEdges,1) -ones(data.numEdges,1)]    ,data.numEdges,data.numVertices);
    
    % compute primal vert/edge laplacian
    data.primalOneLaplacian = data.primalIncidenceMatrix'*data.primalIncidenceMatrix;
    
    %% get creased boundary edges
    BTris = data.triangles(data.isBoundaryTriangle==1,:);
    BTriNormals = data.triangleNormals(data.isBoundaryTriangle==1,:);
    BE2BTri = data.edgesToTrianglesIndicator(data.isBoundaryEdge==1,data.isBoundaryTriangle==1);
    [ii jj] = find(BE2BTri); [~, perm ] = sort(ii); jj=jj(perm);
    BE2BTriPair = reshape(jj,2,[])';
    isCreasedBE = dot(BTriNormals(BE2BTriPair(:,1),:), BTriNormals(BE2BTriPair(:,2),:),2) < cos(anglethresh*pi/180);
%     BE = data.edges(data.isBoundaryEdge==1,:);
    
    data.isCreasedBoundaryEdge = isCreasedBE;
%     data.boundaryEdges = BE;

    BTri2Verts = sparse(data.triangles(data.isBoundaryTriangle==1,:),repmat(1:sum(data.isBoundaryTriangle),1,3),ones(3*sum(data.isBoundaryTriangle),1),data.numVertices,sum(data.isBoundaryTriangle));
    BtriNormals = data.triangleNormals(data.isBoundaryTriangle==1,:);
    [jj,ii]=find(BTri2Verts(data.isBoundaryVertex,:));
    [~,ic]=unique(jj);
    data.vertNorms = BtriNormals(ii(ic),:);
    
    areaWeightedBTri2BVerts = BTri2Verts(data.isBoundaryVertex,:).*data.triangleAreas(data.isBoundaryTriangle)';
    data.areaWeightedVertNorms = areaWeightedBTri2BVerts * BtriNormals;
    data.areaWeightedVertNorms = data.areaWeightedVertNorms./vecnorm(data.areaWeightedVertNorms,2,2);
    
    %{
    figure; axis equal; hold all; rotate3d on;
    patch('faces',data.triangles(data.isBoundaryTriangle,:),'vertices',X,'faceColor','green','edgecolor','none')
    quiver3(X(data.isBoundaryVertex,1),X(data.isBoundaryVertex,2),X(data.isBoundaryVertex,3),data.areaWeightedVertNorms(:,1),data.areaWeightedVertNorms(:,2),data.areaWeightedVertNorms(:,3),'r','linewidth',1)
    %}
    
    data.edgeLengths = vecnorm(data.vertices(data.edges(:,1),:)-data.vertices(data.edges(:,2),:),2,2);
    
    %% compute sparse linear matrix that computes transition from linear vertex based function to constant tet based function.
    p1234 = permute(reshape([data.vertices(data.tetrahedra(:,1),:)'; ones(1,data.numTetrahedra);...
    data.vertices(data.tetrahedra(:,2),:)'; ones(1,data.numTetrahedra);...
    data.vertices(data.tetrahedra(:,3),:)'; ones(1,data.numTetrahedra);...
    data.vertices(data.tetrahedra(:,4),:)'; ones(1,data.numTetrahedra)],4,4,[]),[2 1 3]);
    linearVertsToConstantTets = multinv(p1234);
    linearVertsToConstantTetsGradOnly = linearVertsToConstantTets(1:3,:,:);
    II = permute(reshape(repmat([1:3*data.numTetrahedra]',1,4),3,[],4),[1 3 2]);
    JJ = repelem(data.tetrahedra,1,3)';
    linearVertsToConstantTetsOp = sparse(II(:),JJ(:),linearVertsToConstantTetsGradOnly(:),3*data.numTetrahedra, data.numVertices);
    data.linearVertsToConstantTetsOp = linearVertsToConstantTetsOp;
    
    %{
    %% validate linearVertsToConstantTetsOp operator is correct
    randf = randn(data.numVertices,1)*100;
    randdf = reshape(data.linearVertsToConstantTetsOp*randf,3,[])';
    for iii=1:data.numTetrahedra
        randtet = iii;
        tvi = data.tetrahedra(randtet,:);
        tv = data.vertices(tvi,:);
        fv = randf(tvi);
        aff_f = [tv ones(4,1)]\fv;
        gradf = aff_f(1:3);
        assert(norm(randdf(randtet,:)'-gradf)<.00001)
    end
    %}
    data.vertexWeights = (data.tetVolumes'*data.tetsToVertsIndicator/4)';
    
    data.tetsToTrianglesIndicator = sparse(data.tetsToTriangles(:),repmat(1:data.numTetrahedra,1,4),ones(4*data.numTetrahedra,1),data.numTriangles,data.numTetrahedra)';
    
    % hessian pattern to be given to fmincon if you want to save memory and
    % have a reasonable objective function.
    data.PrimalHessPattern1 = data.primalOneLaplacian~=0;
%    data.PrimalHessPattern2 = (data.PrimalHessPattern1*data.PrimalHessPattern1)~=0;
    data.DualHessPattern1 = (data.tetsToTrianglesIndicator*data.tetsToTrianglesIndicator')~=0;
%    data.DualHessPattern2 = (data.DualHessPattern1*data.DualHessPattern1)~=0;
    
    data.tet2vertShuffler{1} = sparse(data.tetrahedra(:,1),1:data.numTetrahedra,ones(data.numTetrahedra,1),data.numVertices,data.numTetrahedra);
    data.tet2vertShuffler{2} = sparse(data.tetrahedra(:,2),1:data.numTetrahedra,ones(data.numTetrahedra,1),data.numVertices,data.numTetrahedra);
    data.tet2vertShuffler{3} = sparse(data.tetrahedra(:,3),1:data.numTetrahedra,ones(data.numTetrahedra,1),data.numVertices,data.numTetrahedra);
    data.tet2vertShuffler{4} = sparse(data.tetrahedra(:,4),1:data.numTetrahedra,ones(data.numTetrahedra,1),data.numVertices,data.numTetrahedra);
    
    
    % dual laplacian
    idx = all(data.triangles2tets~=0,2);
    reltris = data.triangles2tets(idx,:);
    D = sparse(repmat(1:size(reltris,1),1,2)', reltris(:), [ones(size(reltris,1),1);-ones(size(reltris,1),1)]);
    reltriAreas = data.triangleAreas(idx);
    dualEdgeLen = vecnorm(data.tetBarycenters(reltris(:,1),:) - data.tetBarycenters(reltris(:,2),:),2,2);
    Mass_d1_p2 = diag(sparse(reltriAreas./dualEdgeLen));
    Mass_p3_d0 = diag(sparse(1./ data.tetVolumes));
    data.DualLap = Mass_p3_d0*D'*Mass_d1_p2*D;
    %data.DualLap = D'*D;
    % fiedler test; [a,b] = eigs(data.DualLap',2,'sm'); fiedler = a(:,2); scatter3(data.tetBarycenters(:,1),data.tetBarycenters(:,2),data.tetBarycenters(:,3),10,fiedler,'filled')
    
end









