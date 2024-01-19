function exportMESH(filename, verts, tets)

vertTable = [array2table(["MeshVersionFormatted", "1", "", "";
                 "Dimension", "3", "", "";
                 "Vertices", size(verts,1), "", ""]);
             array2table([verts zeros(size(verts, 1), 1)])];
writetable(vertTable, filename, 'FileType', 'text', 'Delimiter', ' ', 'WriteVariableNames', false);

nt = size(tets, 1);
tetTable = [array2table(["Tetrahedra", nt, "", "", ""]);
            array2table([tets zeros(nt, 1)])];

tempFilename = filename + ".temp";
% tempFilename = "output_frames_dir/blahblah.tmp";
writetable(tetTable, tempFilename, 'FileType', 'text', 'Delimiter', ' ', 'WriteVariableNames', false);

% writetable(tetTable, , 'FileType', 'text', 'Delimiter', ' ', 'WriteVariableNames', false);
cat = 'cat ' + tempFilename + ' >> ' + filename;
system(cat);

fid = fopen(filename, 'at');
fprintf(fid, 'End');
fclose(fid);

system('rm ' + tempFilename);
end