function labelsC=uniformLabel(word,vector,pos)
% uniformLabel allows you to create a string output by combining a certain
% word and vector in the range of limited size of character.
% Therefore, it creates a set of lables and make sure they all have the same length.
% (ex: A001,A019...and so on)//Optional to make word as prefix or suffix. 
% word: 1={a character to add on} // 2={a character to fill} // 3={limit length}
% vector: vector of numerical or character values
% pos= {where to put a word to(i.e., 'last', 'first')}
% EX):uniformLabel({'A','0',3},subjectID,'first')
% 


% Converting vector into a cell array and make sure they are characters
vecCell=num2cell(vector);
if ~ischar(vector),vecCell=cellfun(@(x) num2str(x),vecCell,'UniformOutput',false);end;

for i=1:length(vecCell)
    filler=word{2};limit=word{3};
    numnum=numplace(str2double(vecCell{i}));
    %keyboard
    fillerC=repmat(filler,1,limit-numnum);
    
    switch pos
        case {'first'},labelsC{i}=strcat(word{1},fillerC,vecCell{i});
        case {'last'},labelsC{i}=strcat(vecCell{i},word{1},fillerC);
    end
end

end