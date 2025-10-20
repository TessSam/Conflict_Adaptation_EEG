function tf_savebins(filename,EEG,datareq,varargin)
disp('~~~~~Saving data as .bin file!!!~~~~~');

% Default setting of triming is same as epoch duration
[tepochdur]=process_options(varargin,'tepochdur',[min(EEG.BINEPOCH.times),max(EEG.BINEPOCH.times)]);
ti=dsearchn(EEG.BINEPOCH.times',tepochdur');
EEG.BINEPOCH.times=EEG.BINEPOCH.times(ti(1):ti(2));%update time!
datafields={'eegraw','eegpower','eegphase'};
datareq=ismember(datafields,datareq);

for f = 1:length(datafields)
    switch datafields{f}
        case {'eegraw'},faddress='BINEPOCH';tDim=2;
        otherwise,faddress='BINEPOCH.WAVELET';tDim=3;
    end
    
    %Acess data and erace original data!
    data=eval(strcat('EEG.',faddress,'.(datafields{f})'));
    data=idxDR(data,[tDim],{ti(1):ti(2)});%trim data
    EEG.BINEPOCH.header.(datafields{f}) = size(data);%keep it above WAVELET
    eval(strcat('EEG.',faddress,'.',datafields{f},'=[];'));
    
    %saving
    if datareq(f)
        data = reshape(data,size(data,1),[]);
        dfile = strcat(filename, '_', datafields{f}, '.bin');
        binfile = fopen(dfile,'w');
        fwrite(binfile,data,'single');
        fclose(binfile);
        fprintf('Saved file %s\n',dfile);
    end
end

% Save EEG_set file to store setting!
matfile = strcat(filename, '.mat');
save(matfile,'EEG','-v7.3');
fprintf('Saved file %s\n',matfile);

end