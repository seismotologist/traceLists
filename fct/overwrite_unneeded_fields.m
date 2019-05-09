function overwrite_unneeded_fields(trList,fieldList)

% Because traceList.m is an object, fields cannot be removed. But you can
% still empty the space by replacing the field with a single integer.

verbose = false;

if verbose; fprintf(1,'Full size\t\t\t')
            trList.printObjectSize; 
end

%[~] = trList.printObjectSize;
ntr = numel(trList.eq.m);

% Separate regular fields and those from scalFeature
isScal          = logical(cellfun(@(x) sum(x), regexp(fieldList,'scalFeature')));
rmScalFieldList = fieldList( isScal);
rmRegFieldList  = fieldList(~isScal);

% Remove regular fields
nregfield = numel(rmRegFieldList);
for ifd = 1:nregfield
    
    thisField = rmRegFieldList{ifd};
    try
        eval(['trList.',thisField,' = -1*ones(1,1,''uint8'');'])
        if verbose; fprintf(1,['Removed field ',thisField,'\t\t'])
                    trList.printObjectSize; 
        end
    catch
    end
end

% Remove fields from scalFeature
%fprintf(1,sprintf('Removing field(s) %s\n',strjoin(rmScalFieldList)))
ptIdxList       = cell2mat(regexp(rmScalFieldList,'\.'));
rmScalFieldList = cellfun(@(x) x(ptIdxList+1:end), rmScalFieldList,'uniformOutput',0);
 
%scalTmp = trList.scalFeature;
nscalfield = numel(rmScalFieldList);
for ifd = 1:nscalfield 
    thisField = rmScalFieldList{ifd};
    try
        for itr = 1:ntr
            trList.scalFeature{itr} = rmfield(trList.scalFeature{itr},thisField);
        end
        if verbose; trList.printObjectSize; 
                    fprintf(1,['Removed field scalFeature.',thisField,'\t\t'])
        end
        %scalMod = cellfun(@(x) rmfield(x,rmScalFieldList{ifd}), scalTmp);
    catch
    end
end