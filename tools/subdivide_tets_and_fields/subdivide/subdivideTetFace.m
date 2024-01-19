

function [sub_mesh, sub_fra, mint_fra] = subdivideTetFace(in_mesh, in_field, out_mesh, visualize, scale)

    if ~exist('out_mesh', 'var')
%         in_mesh = 'meshes/octahedron_2.mesh';
%         in_field = "output_frames_dir/subdiv_2.fra";
 in_mesh = 'meshes/octahedron_8.mesh';
        in_field = "output_frames_dir/subdiv.fra";

% in_mesh = 'output_frames_dir/sphere17_5.mesh';
% in_field = 'output_frames_dir/sphere17_balign_d32_0.6.fra';

        
        out_mesh = "output_frames_dir/subdiv_out.mesh";
%         out_field = "output_frames_dir/subdiv_out.fra";
        visualize = true;
    end
    
    if ~exist('visualize','var')
        visualize = 0;
    end
    
    if ~exist('scale', 'var')
        scale = 1;
    end

    %% Missing: 
    %% correct discrete transfer field across subdivision rule
    %% 
    %%%%
    
    
    
    
    
    mesh = ImportMesh(in_mesh);

    tets = mesh.tets;
    ntets = size(tets,1);

    verts = mesh.verts;
    nverts = size(verts,1);

    tC =  ( verts( tets(:,1), :) + verts( tets(:,2), :) + ...
                     verts( tets(:,3), :) + verts( tets(:,4), :) ) / 4;
%                  
%     verts = [ verts; tC ];
    
    [internalEdges, boundaryTris] = genDualMesh(tets, verts);
    

    
    
    
    
    %% make trimesh
    
    % insert extra vertex
    tet_tris(:,:,1) = tets(:, [1 2 3]);
    tet_tris(:,:,2) = tets(:, [1 2 4]);
    tet_tris(:,:,3) = tets(:, [1 3 4]);
    tet_tris(:,:,4) = tets(:, [2 3 4]);
    
    tet_tris = permute(tet_tris, [1 3 2]);
    tet_tris = reshape(tet_tris, 4*ntets, 3);
    
    % Need to establish order in which extra verts will be stored. 
    
    % establish ids for the added vertex per triangle barycenter.
    
    triToDualVert = zeros( ntets, 4 ) - nverts; % sub to make look like 0.

    btIdxs = find(boundaryTris);
    triToDualVert(btIdxs) = 1:size(btIdxs);
    bts = tet_tris(btIdxs, :);
    
    btC =  ( verts( bts(:,1), :) + verts( bts(:,2), :) + ...
             verts( bts(:,3), :) ) / 3;
         
%              triToDualVert2 = zeros( ntets, 4 ); 
    intIdxs = sub2ind( size(triToDualVert), internalEdges(:,1), internalEdges(:,2) );
    triToDualVert(intIdxs) = (1:size(intIdxs)) + size(btIdxs, 1);
    
    intIdxs2 = sub2ind( size(triToDualVert), internalEdges(:,3), internalEdges(:,4) );
    triToDualVert(intIdxs2) = (1:size(intIdxs)) + size(btIdxs, 1);
    
    its = tet_tris(intIdxs, :);
    itC =  ( verts( its(:,1), :) + verts( its(:,2), :) + ...
             verts( its(:,3), :) ) / 3;
    
    verts = [ verts; btC; itC ];
    triToDualVert = triToDualVert + nverts;

    
    
    
if visualize                   
    figure; hold all; axis image vis3d; rotate3d on;
    xlabel('XXXXX')
    ylabel('YYYYY')
    zlabel('ZZZZZ')
%     scatter3(tC(:,1),tC(:,2),tC(:,3),'filled');
    scatter3(verts(:,1), verts(:,2), verts(:,3), 'filled');
    
    scatter3(btC(:,1),btC(:,2),btC(:,3),'filled');   
    scatter3(itC(:,1),itC(:,2),itC(:,3),'filled');   
end 
    
%     [rows, cols] = ind2sub([ntets 4], find(boundaryTris))
%     bts = tet_tris(rows, :, cols)


%% Subdivide tet ids.  This should be modified for selective subd 
%% 
%%

    tets_proc = repelem(tets, 1, 1, 11);
%     tets_proc = flip(tets_proc, 3);
%     tets_proc = permute(tets_proc, [1 3 2]);
    

    
    for i = 1:ntets 
       
        p = tets(i,:);
        d = triToDualVert(i, :);
        
        tets_proc(i,3,1) = d(2);
        tets_proc(i,4,1) = d(1);
        
        tets_proc(i,2,2) = d(3);
         tets_proc(i,4,2) = d(1);
        
        tets_proc(i,1,3) = d(4);
        tets_proc(i,4,3) = d(1);
        
        tets_proc(i,2,4) = d(3);
        tets_proc(i,3,4) = d(2);
        
        tets_proc(i,1,5) = d(4);
        tets_proc(i,3,5) = d(2);
        
        tets_proc(i,1,6) = d(4);
        tets_proc(i,2,6) = d(3);
        
        tets_proc(i,:,7) = d;
        
        tets_proc(i,:,8) = d;
        tets_proc(i,:,9) = d;
        tets_proc(i,:,10) = d;
        tets_proc(i,:,11) = d;
        
        tets_proc(i,1,8) = p(4);
        tets_proc(i,2,9) = p(3);
        tets_proc(i,3,10) = p(2);
        tets_proc(i,4,11) = p(1);
        
    end
    
    % map frame face to subdivided tet id
    f2sd = zeros(4,3);
    f2sd(1,:) = [ 1 2 3 ];
    f2sd(2,:) = [ 1 4 5 ];
    f2sd(3,:) = [ 2 4 6 ];
    f2sd(4,:) = [ 3 5 6 ];
    
    
    tets_subd = permute(tets_proc, [1 3 2]);
    tets_subd = reshape(tets_subd, 11*ntets, 4);
    stC = ( verts( tets_subd(:,1), :) + verts( tets_subd(:,2), :) + ...
    verts( tets_subd(:,3), :) + verts( tets_subd(:,4), :) ) / 4;
    
    
if visualize
%     tets_proc = tets_proc(:,:,7);
    
    viz_tets_subd_face = tets_proc(:,:,1:6);
    viz_tets_subd_center = tets_proc(:,:,7);
    viz_tets_subd_joint = tets_proc(:,:,8:11);

    scatter3(stC(:,1),stC(:,2),stC(:,3),'filled');
    
    for i = 1:6
       show = i;
        patch('Faces',[tets_proc(:,1:3, show); 
                   tets_proc(:,2:4, show); 
                   tets_proc(:,[1 3 4], show); 
                   tets_proc(:,[1 2 4], show)],'Vertices',verts,...
        'EdgeColor','red','FaceColor','none','LineWidth',1);
        
    end

    show = 7;
    patch('Faces',[tets_proc(:,1:3, show); 
               tets_proc(:,2:4, show); 
               tets_proc(:,[1 3 4], show); 
               tets_proc(:,[1 2 4], show)],'Vertices',verts,...
    'EdgeColor','green','FaceColor','none','LineWidth',7);
    
    for i = 8:11
       show = i;
        patch('Faces',[tets_proc(:,1:3, show); 
                   tets_proc(:,2:4, show); 
                   tets_proc(:,[1 3 4], show); 
                   tets_proc(:,[1 2 4], show)],'Vertices',verts,...
        'EdgeColor','blue','FaceColor','none','LineWidth',3);
        
    end
    
        figure; hold all; axis image vis3d; rotate3d on;
    xlabel('XXXXX')
    ylabel('YYYYY')
    zlabel('ZZZZZ')
        show = 7;
    patch('Faces',[tets_proc(:,1:3, show); 
               tets_proc(:,2:4, show); 
               tets_proc(:,[1 3 4], show); 
               tets_proc(:,[1 2 4], show)],'Vertices',verts,...
    'EdgeColor','green','FaceColor','none','LineWidth',7);
        
end


%%
%%%%
%%
%% Missing: Subd border elements.
%%
%%
%%%%


    

%%
%%%%
%%
%% Subdivide the field.  
%% 
%% For now lets do uniform resampling + reproject







    field = importFRA(in_field);
    % lf = Lfield_R9T_R22T(field);
    % lf = reshape(lf, 22, []);
     
    field_subdiv = repelem(field, 1, 1, 11);

    
    %{
    %%Easy speedup: subdivide fields on subdiv.
    
% the edge adjacent tets sample thier neighbors
% the interior tets average these
% and the center tet is simply a constant map of the upper level tet.
    
% partially implemented below:
    
    neighborCounts = ones(ntets, 11);    
    
    for i = 1:size(internalEdges,1)
        e = internalEdges(i, :);
        
        sdtets = f2sd(e(2), :);
        lf_out(:, e(1), sdtets) = lf_out(:, e(1), sdtets) + repelem(lf(:, e(3)), 1, 1, 3);
        neighborCounts( e(1), sdtets ) = neighborCounts( e(1), sdtets ) + 1;
        
        sdtets = f2sd(e(4), :);
        lf_out(:, e(3), sdtets) = lf_out(:, e(3), sdtets) + repelem(lf(:, e(1)), 1, 1, 3);
        neighborCounts( e(3), sdtets ) = neighborCounts( e(3), sdtets ) + 1;
        
        
    end
    
    %}
    
    ntp = size(tets_subd, 1);
    
        vis_edges = zeros(ntp,1);
    

%     ave_frames = pi_gl3_fmin_freeboundary(lf_out, ...
%                                           field_subdiv ); 
                                
    temp = 1:(ntp);
    edges_subd = temp( vis_edges > 0 )';
    edges_subd = [ edges_subd vis_edges(vis_edges > 0 )];
    
    
    orient = dot(cross(verts(tets_subd(:,2), :)-verts(tets_subd(:,1), :),  ...
                       verts(tets_subd(:,3), :)-verts(tets_subd(:,1), :)), ...
                       verts(tets_subd(:,4), :)-verts(tets_subd(:,1), :), 2);
%     orient = -orient;
    
    to_flip = tets_subd( orient < 0, : );
    tets_subd( orient < 0, 3 ) = to_flip(:, 4);
    tets_subd( orient < 0, 4 ) = to_flip(:, 3);
   
%%
%% Save field
%%
    
    
    [filepath,name,ext] = fileparts(out_mesh); 
    output_name = filepath + "/" + name;
    sub_mesh = output_name + ".mesh";
    sub_fra = output_name + ".fra";
    
    mint_fra = output_name + "_int";
    

    exportMESH(sub_mesh, verts, tets_subd);
    exportFRA(field_subdiv(:) * scale, ntp, "--gl3", sub_fra);

%      exportFRA(ave_frames, ntp, "--gl3", sub_int_fra);
    
     
     %%
     %% Call mint, why not.  
     %%
     
     
% override_params.opt_steps = 40;
% override_params.post_steps = 20;
% override_params.global_solve_type = '--pgmres';
% 
% override_params.proj_solve_type = '--jennrich';
% override_params.opt_frame_type = '--fit 500'; % 
% override_params.delta_metric = '--M_dual';  
% % TODO: Add stress align term.
% override_params.lambda_delta = 1;
% override_params.lambda_fit_eye = .0001;  
% override_params.lambda_fit_rescale = 10;  
% 
% 
% % override_params.anneal_exp_arg = 4;
% override_params.opt_lambdas = [ 32 .00001 0 ];
% 
% % The smooth limit. 
% % override_params.lambda_delta = 0;
% override_params.output_slug = mint_fra;
% 
% mint_fit(sub_mesh, sub_fra, override_params);
%    sub_fra =   output_name;


    
%%
%% Show field
%%
     tets_proc = tets_subd;
if visualize  
                tC =  ( verts( tets_proc(:,1), :) + verts( tets_proc(:,2), :) + ...
                     verts( tets_proc(:,3), :) + verts( tets_proc(:,4), :) ) / 4;
                 
    visualize_frameField(field_subdiv(:), tC, false, 1);
%     visualize_frameField(ave_frames, tC, false, 1);
    
%     figure; hold all; axis image vis3d; rotate3d on;
%     xlabel('XXXXX')
%     ylabel('YYYYY')
%     zlabel('ZZZZZ')
%     visualize_frameField(importFRA(sub_int_fra), tC, false, 1);
%     scatter3(tC(:,1),tC(:,2),tC(:,3),'filled');
    scatter3(verts(:,1), verts(:,2), verts(:,3), 'filled');
    
    

    scatter3(tC(:,1),tC(:,2),tC(:,3),'filled');   
          
          

patch('Faces',[tets_proc(:,1:3); 
               tets_proc(:,2:4); 
               tets_proc(:,[1 3 4]); 
               tets_proc(:,[1 2 4])],'Vertices',verts,...
    'EdgeColor','red','FaceColor','none','LineWidth',2);
          
          patch('Faces',[tets(:,1:3); 
               tets(:,2:4); 
               tets(:,[1 3 4]); 
               tets(:,[1 2 4])],'Vertices',verts,...
    'EdgeColor','green','FaceColor','none','LineWidth',2);
    
    
    frames = importFRA(override_params.output_slug + ".fra");
%     frames = ave_frames; %field_subdiv;
    frames = frames(:);
    tetCenters = [tC; tC];
    scale = 1;
    
    nframes = size(frames,1)/9;
    frames_R3x3 = reshape(frames,[3,size(frames,1)/3])';
    for i = 1:nframes
        frames_R3x3((3*(i-1)+1):3*i, :) = frames_R3x3((3*(i-1)+1):3*i, :)';
    end
    
%     sFrame = 1;
%     nData  = size(frames,1)/3 ;
%     dim = 3;
%     M_tr = cell2mat(mat2cell(frames_R3x3, dim, nData*ones(sFrame,1)).').'

    frames_R3x3 = frames_R3x3 * scale;
    vecs1 = frames_R3x3(:,1);
    vecs2 = frames_R3x3(:,2);
    vecs3 = frames_R3x3(:,3);
    vecs1 = reshape(vecs1, [3, nframes])';
    vecs2 = reshape(vecs2, [3, nframes])';
    vecs3 = reshape(vecs3, [3, nframes])';
    
    patch('Vertices',tC,'Faces',edges_subd(:,[1 2 1]));
    
% %     figure; hold all; axis image vis3d; rotate3d on;
% %     xlabel('XXXXX')
% %     ylabel('YYYYY')
% %     zlabel('ZZZZZ')
%     scatter3(tC(:,1),tC(:,2),tC(:,3));
%         quiver3(tetCenters(:,1),tetCenters(:,2),tetCenters(:,3),... 
%             [vecs1(:,1); -vecs1(:,1)], ...
%             [vecs1(:,2); -vecs1(:,2)], ...
%             [vecs1(:,3); -vecs1(:,3)] ,'r');
%         quiver3(tetCenters(:,1),tetCenters(:,2),tetCenters(:,3),... 
%             [vecs2(:,1); -vecs2(:,1)], ...
%             [vecs2(:,2); -vecs2(:,2)], ...
%             [vecs2(:,3); -vecs2(:,3)],'g');
%         quiver3(tetCenters(:,1),tetCenters(:,2),tetCenters(:,3),... 
%             [vecs3(:,1); -vecs3(:,1)],...
%             [vecs3(:,2); -vecs3(:,2)],...
%             [vecs3(:,3); -vecs3(:,3)],'b');
        
end


    
end