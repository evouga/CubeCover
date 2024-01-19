function idx = tetTriToIdx_vectorized(sq)
    sq = sort(sq,2);
    [isFound, idx] = ismember(sq,[1 2 3; 1 2 4; 1 3 4; 2 3 4],'rows');
    if ~(all(isFound))
        error('unhandled case. inputs should represent intersection vertex indices between two adjacent tets. there should always be three and they should be valued 1 to 4.');
    end
end