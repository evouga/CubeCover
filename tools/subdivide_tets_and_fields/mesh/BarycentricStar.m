function [star1, star2] = BarycentricStar(meshData)

tets = meshData.tets;
nt = meshData.nt;
faces = meshData.faces;
nf = meshData.nf;
verts = meshData.verts;
nv = meshData.nv;

%% Compute Star 2

intTri = sum(abs(meshData.d2), 1) == 2;
bdryTri = ~intTri;
tetCenters = squeeze(sum(reshape(verts(tets, :), nt, 4, 3), 2)) / 4;
triCenters = squeeze(sum(reshape(verts(faces, :), nf, 3, 3), 2)) / 3;
dualLengths = zeros(nf, 1);
dualLengths(intTri) = vecnorm(meshData.d2(:, intTri)' * tetCenters, 2, 2);
dualLengths(bdryTri) = vecnorm(abs(meshData.d2(:, bdryTri))' * tetCenters - triCenters(bdryTri, :), 2, 2);

a = verts(faces(:, 1), :);
b = verts(faces(:, 2), :);
c = verts(faces(:, 3), :);
areas = 0.5 * vecnorm(cross(a - b, c - b), 2, 2);

star2 = spdiags(dualLengths ./ areas, 0, nf, nf);

%% Compute Star 1

tetToEdge = nchoosek(1:4, 2);
tetEdges = reshape(tets(:, tetToEdge), nt, 6, 2);
tetEdgesOpp = flip(tetEdges, 2);

v0 = verts(tetEdges(:, :, 1), :);
v1 = verts(tetEdges(:, :, 2), :);
w0 = verts(tetEdgesOpp(:, :, 1), :);
w1 = verts(tetEdgesOpp(:, :, 2), :);

v = 0.5 * (v0 + v1);
w = 0.5 * (w0 + w1);
c = 0.5 * (v + w);
tetEdgeDualAreas = 0.5 * (vecnorm(cross(0.5 * (w0 - v), c - v), 2, 2) + vecnorm(cross(c - v, 0.5 * (w1 - v)), 2, 2));
tetEdgeLengths = vecnorm(v1 - v0, 2, 2);

L = sparse(tetEdges(:, :, 1), tetEdges(:, :, 2), tetEdgeDualAreas ./ tetEdgeLengths, nv, nv);
star1 = spdiags(nonzeros(L'), 0, meshData.ne, meshData.ne);

end