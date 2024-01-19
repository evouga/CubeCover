function [] = exportFRA(f, nTets, frametype, filename)

if nargin == 0
    f = randi(10,18,1);
    nTets = 2;
    frametype = "--gl3";
    filename = "output_frames_dir/blah.fra";
end

% nvecs = size(f,1) / nTets;
nrows = size(f,1)/3;

out = reshape(f, [ 3, nrows ]);

frametypeId = 0;
if frametype == "--octa"
    frametypeId = 4;
end
    
    
framesPerTet = 3;
% 	0: 3 independent (globally-combed) vector fields
% 	1: R^(3 x 3)/O field
% 	2: fradeco field
% 	3: odeco field
% 	4: unit octahedral field
%   frametype = frametype;

[filepath,name,ext] = fileparts(filename); 
[~,~,~] = mkdir(filepath);

fileID = fopen(filename,'w');
fprintf(fileID,'FRA 1\n');
fprintf(fileID, '%i %i %i\n', nTets, framesPerTet, frametypeId);
%fprintf(fileID,'%6s %12s\n','x','exp(x)');
fprintf(fileID,'%0.10f %0.10f %0.10f\n',out);
fclose(fileID);

end