function meshData = ProcessMesh(verts, tets, bdryAngleCutoff)

if nargin < 3
    bdryAngleCutoff = 0;
end

meshData.verts = verts;
meshData.tets = tets;

meshData.nt = size(meshData.tets, 1);
meshData.nv = size(meshData.verts, 1);
meshData.tetra = triangulation(meshData.tets, meshData.verts);

%% Exterior differential structure

[sortedTets, perm] = sort(meshData.tets, 2);
tetSign = levicivita(perm, 2);

[meshData.faces, ~, tet2face] = unique(reshape(sortedTets(:, nchoosek(1:4, 3)), [], 3), 'rows');
meshData.tet2face = reshape(tet2face, [], 4);
meshData.nf = size(meshData.faces, 1);

[meshData.edges, ~, face2edge] = unique(reshape(meshData.faces(:, [1 2; 2 3; 1 3]), [], 2), 'rows');
meshData.face2edge = reshape(face2edge, [], 3);
meshData.ne = size(meshData.edges, 1);

meshData.faceOrientation = [1 1 -1];
meshData.d0 = sparse(repmat((1:meshData.ne).', 1, 2), meshData.edges, repmat([-1 1], meshData.ne, 1), meshData.ne, meshData.nv);
meshData.d1 = sparse(repmat((1:meshData.nf).', 1, 3), meshData.face2edge, repmat(meshData.faceOrientation, meshData.nf, 1), meshData.nf, meshData.ne);
meshData.d2 = sparse(repmat((1:meshData.nt).', 1, 4), meshData.tet2face, tetSign .* [1 -1 1 -1], meshData.nt, meshData.nf);

meshData.edgeCenters = squeeze(sum(reshape(meshData.verts(meshData.edges.', :), 2, [], 3), 1) / 2);
meshData.faceCenters = squeeze((1/3) * sum(reshape(meshData.verts(meshData.faces.', :), 3, [], 3), 1));

%% Boundary

bdryFaces = freeBoundary(meshData.tetra);
bdryIdx = unique(bdryFaces(:));
warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
meshData.bdry = triangulation(bdryFaces, meshData.verts);
warning('on', 'MATLAB:triangulation:PtsNotInTriWarnId');
bdryNormals = vertexNormal(meshData.bdry, bdryIdx);
nb = length(bdryIdx);

% % % % Remove sharp creases (dihedral angle < bdryAngleCutoff)
% % % bdryVtxStars = vertexAttachments(meshData.bdry, bdryIdx);
% % % bdryVtxStarSizes = cellfun(@length, bdryVtxStars);
% % % bdryVtxStarIdx = repelem((1:nb)', bdryVtxStarSizes);
% % % bdryVtxTriIdx = cell2mat(bdryVtxStars')';
% % % bdryVtxStarNormals = repelem(bdryNormals, bdryVtxStarSizes, 1);
% % % bdryVtxTriNormals = faceNormal(meshData.bdry, bdryVtxTriIdx);
% % % bdryVtxCosAngles = dot(bdryVtxStarNormals, bdryVtxTriNormals, 2);
% % % bdryVtxMinCosAngle = accumarray(bdryVtxStarIdx, bdryVtxCosAngles, [], @min);
% % % bdryVtxSmooth = bdryVtxMinCosAngle >= cos(0.5 * (pi - bdryAngleCutoff));
% % % meshData.bdryIdx = bdryIdx(bdryVtxSmooth);
% % % meshData.bdryNormals = bdryNormals(bdryVtxSmooth, :);
% % % meshData.bdryEdgeIdx = find(sum(abs(meshData.d0(:, meshData.bdryIdx)), 2) == 2);
% % % 
% % % meshData.intIdx = setdiff((1:meshData.nv)', meshData.bdryIdx);

%% Laplacian

[meshData.L, meshData.M, meshData.star1] = GeometricPrimalLM(meshData);
meshData.MLump = spdiags(sum(meshData.M, 2), 0, meshData.nv, meshData.nv);

if false && gpuDeviceCount > 0
    % Inverse iteration to find second-smallest eigenvalue L*v = lambda*M*v
    v = randn(meshData.nv, 1, 'gpuArray');
    Lg = gpuArray(meshData.L);
    Mg = gpuArray(meshData.M);
    lambda = 0;
    for i = 1:100
        [v, ~] = pcg(Lg, v, 1e-6, 1000, [], [], v);
        v = v - mean(v);
        v = v / sqrt(v' * Mg * v);
        rayleigh = (v' * Lg * v);
        if abs(rayleigh - lambda) < 1e-5
            lambda = rayleigh;
            [v, ~] = pcg(@(w) Lg * w - lambda .* w, v, 1e-6, 1000, [], [], v);
            v = v / sqrt(v' * Mg * v);
            lambda = (v' * Lg * v);
            break;
        end
        lambda = rayleigh;
    end
    meshData.lambda1L = lambda;
    clear Lg Mg;
else
%     lambda = eigs(meshData.L, meshData.M, 2, 'smallestabs', 'IsSymmetricDefinite', true);
%     meshData.lambda1L = lambda(2);
end

end
