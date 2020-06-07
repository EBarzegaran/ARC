function [ind g4] = selec_hemisphere(g3,Hem)
% function to keep one hemisphere vertices and faces
    % INPUT: 
        % g3: all the vertices and faces in both hemisphere
        % Hem: hemisphere of interest: 'left' or 'right'
    % OUTPUT
        
        
    if strcmp(Hem,'left'),
        VOIk = find(g3.vertices(:,1)<0);%left : indices to keep
        VOIr = find(g3.vertices(:,1)>=0);%right : indices to remove
    else
        VOIk = find(g3.vertices(:,1)>=0);%right : indices to keep
        VOIr = find(g3.vertices(:,1)<0);%left : indices to remove
    end
    VOIr = sort(VOIr);
    VOIk = sort(VOIk);

  g4 = g3;
  g4.vertices(VOIr,:) = [];
  TI = zeros(1,size(g3.vertices,1));% for changing vertice number in the faces
  TI(VOIk)=1:numel(VOIk);
  %% remove extra faces
  Ind = 1;
  while ~isempty(Ind)
      Ind = [];%indices in faces to be keep
      for i = 1:3
            [~,I] = intersect(g4.faces(:,i),VOIr);% find the vertices in one hemisphere
            Ind = [Ind; I];
      end
      Ind = unique(Ind);
      g4.faces(Ind,:)=[];
  end
  g4.faces(:,1) = TI(g4.faces(:,1));g4.faces(:,2) = TI(g4.faces(:,2));g4.faces(:,3) = TI(g4.faces(:,3));
  
  ind = VOIk;
end