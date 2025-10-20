function output=unpackcell(c)
% This function unpacks cells and converts it to a matrix. 
% The cells are summarized into a matrix sized by the dimension of the
% original cell + 1 dimension. For instance, 10 cells of (2,2,2) matrixes
% would be (2,2,2,10)!
% =========================================================================

%Check if the input is cell or not!!
if ~iscell(c),error('Make sure to feed in cell to unpack!');end

tic
c=reshape(c,[],1);
numcell=size(c,1);
wcellsize=size(c{1});
output=zeros(horzcat(wcellsize,numcell));

for i=1:numcell,
    %This part is very dumb
    if length(size(output))==5,output(:,:,:,:,i)=c{i};end;
    if length(size(output))==4,output(:,:,:,i)=c{i};end;
    if length(size(output))==3,
        if numcell~=1
            output(:,:,i)=c{i};
        else
            output(:,:,:,i)=c{i};
        end
    end
end
toc



dim=ndims(c{1})+1;
output=cat(dim,c{:});

end