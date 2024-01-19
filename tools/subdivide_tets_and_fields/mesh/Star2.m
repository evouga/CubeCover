function star2 = Star2(meshData)

intTri = sum(abs(meshData.d2), 1) == 2;
bdryTri = ~intTri;
tetCenters = circumcenter(meshData.tetra, (1:meshData.nt)');
triCenters = circumcenter(triangulation(meshData.faces, meshData.verts), (1:meshData.nf)');
dualLengths = zeros(meshData.nf, 1);
dualLengths(intTri) = vecnorm(meshData.d2(:, intTri)' * tetCenters, 2, 2);
dualLengths(bdryTri) = vecnorm(abs(meshData.d2(:, bdryTri))' * tetCenters - triCenters(bdryTri, :), 2, 2);

a = meshData.verts(meshData.faces(:, 1), :);
b = meshData.verts(meshData.faces(:, 2), :);
c = meshData.verts(meshData.faces(:, 3), :);
areas = 0.5 * vecnorm(cross(a - b, c - b), 2, 2);

star2 = spdiags(dualLengths ./ areas, 0, meshData.nf, meshData.nf);

end