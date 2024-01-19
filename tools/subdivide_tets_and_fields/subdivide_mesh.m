

function [] = subdivide_mesh(mesh_path, fra_path, subd_mesh_path)

    fra_char = char(fra_path);
    fra_extension = fra_char(end-3:end);
    
    assert(sum(fra_extension == '.fra')==4);
    
    scale = 1;


    [subd_mesh, sub_fra, mint_fra] = subdivideTetFace( mesh_path, fra_path, subd_mesh_path, false,  scale);

end