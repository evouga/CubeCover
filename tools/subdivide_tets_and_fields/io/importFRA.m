function frames_R9T = importFRA(frafile)

if nargin == 0
    frafile = "output_frames_dir/tetrahedron_100.fra";
end

[filepath, filename, ext] = fileparts(frafile);

if ~strcmpi(ext, '.fra')
    warning('wrong extension')
end

    fid = fopen(frafile,'r');
    
    version = textscan(fid, 'FRA %d', 'MultipleDelimsAsOne', true);
    metadata = textscan(fid, '%d %d %d', 1, 'MultipleDelimsAsOne', true); %, 'Whitespace', ' \t\n', 'MultipleDelimsAsOne', true);
    frames = textscan(fid, '%f %f %f', 'CollectOutput', true,'MultipleDelimsAsOne', true);
    
    fclose(fid);
    
    frames = double(frames{:});
    frames_R9T = frames';
    frames_R9T = frames_R9T(:);
    

    


end