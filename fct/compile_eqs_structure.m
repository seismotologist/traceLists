function eqs = compile_eqs_structure(trList,skipList)
% Uses trList.eq.name to find records from same earthquake, not like
% previous version of script, which used wform-filenames

if nargin<2; skipList = []; end

o.verbose = 0;

%fprintf(1,'Using hardcoded vp=6km/s')
vp = 6;

dm        = 1e-3;
thousands = linspace(1e3,1e6,1e3);

% A. Assign eqIdx to all entries of trList ................................
ntr          = numel(trList.eq.m);
trList.eq.idx = zeros(ntr,1,'int32');
id           = 0;

for itr = 1:ntr
    
    if trList.eq.idx(itr)==0
        fullName        = trList.fullName{itr};
        eqName          = trList.eq.name{itr};
        ids_thisEvent   = find(strcmp(trList.eq.name,eqName));
        %[ids_thisEvent2] = find_all_event_traces(fullName,trList);  % --> old version; very slow, and made some errors with NGA eqs
        if ~isempty(skipList)
            sids_thisEvent   = find(strcmp(skipList.eq.name,eqName));
            %[sids_thisEvent] = find_all_event_traces_NEW(fullName,skipList);
        end
            
        % Sometimes find_all_event_traces.m misclassifies traces. Check
        % if all associated traces have same magnitudes and if not,
        % assign them to different events. Because sometimes, m-values are
        % different only due to limited precision, magnitude differences
        % need to be compared to m+/-dm
        if isempty(skipList); mVect = trList.eq.m(ids_thisEvent);
        else                  mVect = [trList.eq.m(ids_thisEvent); skipList.eq.m(sids_thisEvent)];
        end
        [mValsRaw,istd] = unique(mVect,'first');
        
        if numel(mValsRaw)==1; mVals = mValsRaw;
        else
            
            % unique.m returns the unique values in sorted order. Undo the
            % sorting such that the order of the eqs-entries is consistent
            % with that of the traceList
            mValsRaw = mVect(sort(istd));
            % trList.eq.m(ids_thisEvent)
            
            mVals = [];
            while ~isempty(mValsRaw)
                thisM    = mValsRaw(1);                                    % First of remaning m-values
                mVals    = [mVals; thisM];                                 % Add it to list
                isSameM  = logical(mValsRaw>thisM-dm &mValsRaw<thisM+dm);  % Identify sufficiently similar entries
                mValsRaw = mValsRaw(~isSameM);                             % Take them out of mValsRaw 
            end 
        end
        
        for im = 1:numel(mVals)
            
            id = id+1;

            % TraceList
            itmp               = find(mVals(im) == trList.eq.m(ids_thisEvent));
            idx                = ids_thisEvent(itmp);
            trList.eq.idx(idx) = id;
            
            % SkipList
            if ~isempty(skipList); itmps                = find(mVals(im) == skipList.eq.m(sids_thisEvent));
                                   idxs                 = sids_thisEvent(itmps);
                                   skipList.eq.idx(idxs) = id;
                 if o.verbose; fprintf(1,sprintf('id = %i --- %i  records (+%i skipped) --- %s --- trList-index = %i/%i\n',id,numel(idx),numel(idxs),eqName,itr,ntr)); end    
            else if o.verbose; fprintf(1,sprintf('id = %i --- %i  records --- %s --- trList-index = %i/%i\n',id,numel(idx),eqName,itr,ntr)); end
            end
        end
    end
end


% B. Unique list of eqIds .................................................
% <eqs.eventId> = list of all event-idndices considered for inference
eqs.eventId = unique(trList.eq.idx);
neq         = numel(eqs.eventId);


% C. Save trace indices for all traces in each eqs.eventId ................
%    sorted wrt/ distance, along with the corresp. theretical p-arrival times
eqs.traceId = cell(neq,1);
eqs.tpx     = cell(neq,1);
eqs.lt      = cell(neq,1);
eqs.date    = cell(neq,1);
eqs.t0      = cell(neq,1);
eqs.region  = cell(neq,1);
eqs.name    = cell(neq,1);
%eqs.mType   = cell(neq,1);
eqs.lat     = zeros(neq,1);
eqs.lon     = zeros(neq,1);
eqs.z       = zeros(neq,1);

for ieq = 1:neq
    if ismember(ieq,thousands); fprintf(1,[num2str(ieq),'/',num2str(neq),'\n']);end
    
    trIdx            = find(eqs.eventId(ieq)==trList.eq.idx);      % traceList index of all traces from target event
    R                = trList.dist.hyp(trIdx);
    [Rstd,srtIdx]    = sort(R);
    tp               = Rstd/vp;
    
    eqs.traceId{ieq} = trIdx(srtIdx);
    eqs.tpx{ieq}     = tp;
    eqs.date{ieq}    = trList.eq.date{trIdx(1)};
    eqs.t0{ieq}      = trList.eq.t0    {trIdx(1)};
    eqs.name{ieq}    = trList.eq.name{trIdx(1)};
    %eqs.mType{ieq}   = trList.eq.mType {trIdx(1)};
    eqs.lat(ieq)     = trList.eq.lat (trIdx(1));
    eqs.lon(ieq)     = trList.eq.lon (trIdx(1));
    eqs.z(ieq)       = trList.eq.z   (trIdx(1));
    
    
    dsn    = trList.dataSetName{trIdx(1)};
    isJapa = strcmp(dsn,'kNet') |strcmp(dsn,'kikNet');
    isCali = strcmp(dsn,'scsn') |strcmp(dsn,'scsnPx');
    isNgaw = strcmp(dsn,'ngawest1');
    isWnch = strcmp(dsn,'wenchuan');
    
    if     isJapa; eqs.region{ieq} = 'japan';
    elseif isCali; eqs.region{ieq} = 'california';
    elseif isNgaw; eqs.region{ieq} = 'ngawest1';
    elseif isWnch; eqs.region{ieq} = 'wenchuan';
    end

    if ~isempty(skipList)
        skipIdx               = find(eqs.eventId(ieq)==skipList.eq.idx);
        Rskip                 = skipList.dist.hyp(skipIdx);
        [RskipStd,skipSrtIdx] = sort(Rskip);
        tpskip                = RskipStd/vp;
        
        eqs.skipId {ieq} = skipIdx(skipSrtIdx);
        eqs.Rskip{ieq}   = RskipStd;
        eqs.tpxSkip{ieq} = tpskip;
    end

end
eqs.eqIdx = trList.eq.idx;
eqs.m     = trList.eq.m(cellfun(@(x) x(1), eqs.traceId));