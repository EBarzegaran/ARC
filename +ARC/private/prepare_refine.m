%%  first level refinement
function D2 = prepare_refine(D1,g,ind1,val1,V1,V2)
%INPUTS:
    % D1 is V1*1 matrix, contains data to be presented on the headplot
    % g is structure, includes vertices and faces before refinement
    % V1 is number of vertices before refinement
    % V2 is number of vertices after refinement
    
    F = g.faces;

    D2 = zeros(size(V2,1));D2(1:V1) = D1;
    for i=0:length(F)-1,
        vert1=F(i+1,1); 
        vert2=F(i+1,2); 
        vert3=F(i+1,3);

        index=ind1{vert1}; vals=val1{vert1};
        verta= vals(index==vert2);
        index=ind1{vert2}; vals=val1{vert2};
        vertb= vals(index==vert3);
        index=ind1{vert3}; vals=val1{vert3};
        vertc= vals(index==vert1);

        D2(verta) = (D2(vert1)+D2(vert2))/2;
        D2(vertb) = (D2(vert2)+D2(vert3))/2;
        D2(vertc) = (D2(vert1)+D2(vert3))/2;
    end

end