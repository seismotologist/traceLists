classdef traceList < handle
    %TRACELIST A list of wform meta-info; everything but the wform itself. 
    %   Version: i39/
    %   Last update: 190509
    %   mmeier@caltech.edu

    % Wishlist for i40
    % . more flexible station field, comp with structures in wformproc scripts
    % . make eq field comp with structures in wformproc scripts
    % . gm field
    % . files field
    % . t0 as datetime object
    % . flags field? e.g. for 'useme' flag
    % . file.name, file.fullname
    % . or: name.full, name.trace, name.xx
    % . one list entry for all three comps?
    % . more var fields

    properties
        
        fullName        % Waveform file name including absolute path
        eq              % Event information (hypocenter, origin time, etc.))
        dist            % Various recording distance measures
        station         % Recording station information
        dataSetName     % Name for identifying data sets of different origins
        px              % Pick times and indices
        fb              % Filter bank information
        pga             % Peak ground acceleration information
        pgv             % Peak ground velocity information
        pgd             % Peak ground displacement information
        mmi             % Modified Mercalli Intensity
        noise           % Pre-signal noise statistics
        comment         % e.g. skipReasons for skipped traces
        var             % Placeholder for unexpected variables 
        scalFeature     % Scalar waveform features
        vectFeature     % Vector waveform features
        prop            % Non-trace-specific info, e.g. fc, stored as first entry
    end
    
    methods
        
        % Invoke object   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
        function obj = traceList(ntr)
            
            obj.fullName    = cell(ntr,1);
            
            obj.eq.lat       = zeros(ntr,1,'single');
            obj.eq.lon       = zeros(ntr,1,'single');
            obj.eq.z         = zeros(ntr,1,'single');
            obj.eq.date      = cell(ntr,1);
            obj.eq.t0        = cell(ntr,1);
            obj.eq.ts        = zeros(ntr,1,'single');
            obj.eq.m         = zeros(ntr,1,'single');
            obj.eq.mType     = cell(ntr,1);
            obj.eq.type      = cell(ntr,1);
            obj.eq.mech      = cell(ntr,1);
            obj.eq.idx       = zeros(ntr,1,'single');
            obj.eq.name      = cell(ntr,1);
            
            obj.dist.epi = zeros(ntr,1,'single');
            obj.dist.hyp = zeros(ntr,1,'single');
            obj.dist.flt = zeros(ntr,1,'single');
            
            obj.station.name   = cell(ntr,1);
            obj.station.nw     = cell(ntr,1);
            obj.station.lat    = zeros(ntr,1,'single');
            obj.station.lon    = zeros(ntr,1,'single');
            obj.station.alt    = zeros(ntr,1,'single');
            obj.station.idx    = zeros(ntr,1,'single');
            obj.station.sr     = zeros(ntr,1,'single');
            obj.station.vs30   = zeros(ntr,1,'single');
            obj.station.bcode  = cell(ntr,1);
            obj.station.icode  = cell(ntr,1);
            obj.station.ocode  = cell(ntr,1);
            obj.station.filter = cell(ntr,1);
            
            obj.prop        = struct;
            obj.dataSetName = cell(ntr,1);
            obj.scalFeature = cell(ntr,1);
            obj.vectFeature = cell(ntr,1);
    
            obj.px.p.t     = zeros(ntr,1,'single');
            obj.px.p.idx   = zeros(ntr,1,'single');
            obj.px.p.hasPx = false(ntr,1);
            obj.px.s.t     = zeros(ntr,1,'single');
            obj.px.s.idx   = zeros(ntr,1,'single');
            obj.px.s.hasPx = false(ntr,1);
            obj.px.allPx   = cell(ntr,1);
            
            obj.noise.snr = zeros(ntr,1,'single');
            obj.noise.nb  = cell(ntr,1);
            obj.noise.dsp = cell(ntr,1);
            obj.noise.vel = cell(ntr,1);
            obj.noise.acc = cell(ntr,1);
            
            obj.comment     = cell(ntr,1);
            
            obj.var.v1      = cell(ntr,1);
            obj.var.v2      = cell(ntr,1);
            obj.var.v3      = cell(ntr,1);
            obj.var.v4      = cell(ntr,1);
            obj.var.v5      = cell(ntr,1);
            obj.var.v6      = cell(ntr,1);
            obj.var.v7      = cell(ntr,1);
            obj.var.v8      = cell(ntr,1);
            
            obj.fb.amax     = cell(ntr,1); % Maximum amplitudes in each frequency band
            obj.fb.cav      = cell(ntr,1); % Cumulative absolute velocities in each frequency band

            obj.pga.tsi   = cell(ntr,1);
            obj.pgv.tsi   = cell(ntr,1);
            obj.pgd.tsi   = cell(ntr,1);
            obj.pga.tszen = cell(ntr,1);
            obj.pgv.tszen = cell(ntr,1);
            obj.pgd.tszen = cell(ntr,1);
            
            obj.pgv.up = zeros(ntr,1,'single');
            obj.pgv.lo = zeros(ntr,1,'single');
            obj.pga.nb = cell(ntr,1);
            obj.pgv.nb = cell(ntr,1);
            obj.pgd.nb = cell(ntr,1);
            
            obj.mmi.tsi        = cell(ntr,1);
            obj.mmi.tszen      = cell(ntr,1);
            obj.mmi.pns.ampi   = zeros(ntr,1);
            obj.mmi.pns.ampzen = zeros(ntr,1);
            obj.mmi.pns.ampen  = zeros(ntr,1);
            %obj.mmi.iexi   = cell(ntr,1);
            %obj.mmi.iexzen = cell(ntr,1);
           
            obj.pga.p.ampi   = zeros(ntr,1,'single');
            obj.pga.p.idxi   = zeros(ntr,1,'single');
            obj.pga.p.ampzen = zeros(ntr,1,'single');
            obj.pga.p.idxzen = zeros(ntr,1,'single');
            obj.pga.p.ampen  = zeros(ntr,1,'single');
            obj.pga.p.idxen  = zeros(ntr,1,'single');
            
            obj.pgv.p.ampi   = zeros(ntr,1,'single');
            obj.pgv.p.idxi   = zeros(ntr,1,'single');
            obj.pgv.p.ampzen = zeros(ntr,1,'single');
            obj.pgv.p.idxzen = zeros(ntr,1,'single');
            obj.pgv.p.ampen  = zeros(ntr,1,'single');
            obj.pgv.p.idxen  = zeros(ntr,1,'single');
            
            obj.pgd.p.ampi   = zeros(ntr,1,'single');
            obj.pgd.p.idxi   = zeros(ntr,1,'single');
            obj.pgd.p.ampzen = zeros(ntr,1,'single');
            obj.pgd.p.idxzen = zeros(ntr,1,'single');
            obj.pgd.p.ampen  = zeros(ntr,1,'single');
            obj.pgd.p.idxen  = zeros(ntr,1,'single');
            
            obj.pga.s.ampi   = zeros(ntr,1,'single');
            obj.pga.s.idxi   = zeros(ntr,1,'single');
            obj.pga.s.ampzen = zeros(ntr,1,'single');
            obj.pga.s.idxzen = zeros(ntr,1,'single');
            obj.pga.s.ampen  = zeros(ntr,1,'single');
            obj.pga.s.idxen  = zeros(ntr,1,'single');
            
            obj.pgv.s.ampi   = zeros(ntr,1,'single');
            obj.pgv.s.idxi   = zeros(ntr,1,'single');
            obj.pgv.s.ampzen = zeros(ntr,1,'single');
            obj.pgv.s.idxzen = zeros(ntr,1,'single');
            obj.pgv.s.ampen  = zeros(ntr,1,'single');
            obj.pgv.s.idxen  = zeros(ntr,1,'single');
            
            obj.pgd.s.ampi   = zeros(ntr,1,'single');
            obj.pgd.s.idxi   = zeros(ntr,1,'single');
            obj.pgd.s.ampzen = zeros(ntr,1,'single');
            obj.pgd.s.idxzen = zeros(ntr,1,'single');
            obj.pgd.s.ampen  = zeros(ntr,1,'single');
            obj.pgd.s.idxen  = zeros(ntr,1,'single');
            
            obj.pga.pns.ampi   = zeros(ntr,1,'single');
            obj.pga.pns.idxi   = zeros(ntr,1,'single');
            obj.pga.pns.ampzen = zeros(ntr,1,'single');
            obj.pga.pns.idxzen = zeros(ntr,1,'single');
            obj.pga.pns.ampen  = zeros(ntr,1,'single');
            obj.pga.pns.idxen  = zeros(ntr,1,'single');
            
            obj.pgv.pns.ampi   = zeros(ntr,1,'single');
            obj.pgv.pns.idxi   = zeros(ntr,1,'single');
            obj.pgv.pns.ampzen = zeros(ntr,1,'single');
            obj.pgv.pns.idxzen = zeros(ntr,1,'single');
            obj.pgv.pns.ampen  = zeros(ntr,1,'single');
            obj.pgv.pns.idxen  = zeros(ntr,1,'single');
            
            obj.pgd.pns.ampi   = zeros(ntr,1,'single');
            obj.pgd.pns.idxi   = zeros(ntr,1,'single');
            obj.pgd.pns.ampzen = zeros(ntr,1,'single');
            obj.pgd.pns.idxzen = zeros(ntr,1,'single');
            obj.pgd.pns.ampen  = zeros(ntr,1,'single');
            obj.pgd.pns.idxen  = zeros(ntr,1,'single');
        end
        
        
        
        
        % -----------------------------------------------------------------
        % Select certain entries of GlobalList and output them as a new list
        % note: can also be used to obtain independent clone list to which
        % changes can be made without changing the original list
        function [subList] = selectSubList(obj,selectIndex) 
            if (islogical(selectIndex))
                error('Use indices instead of logical vectors to select entries from traceList.')
            end
            names   = fieldnames(obj);
            nf      = size(names,1);
            ntr     = numel(obj.(names{1}));
            ntr_sub = numel(selectIndex);
            subList = traceList(ntr_sub);      %generate new object
            for ifd = 1:nf

                fdName = names{ifd};
                
                % Determine field type
                isPropField   = strcmp(fdName,'prop');                  % Field with one structure entry per list, not per trace
                isEmptyField  = isempty(obj.(fdName));
                isSingleField = ~isstruct(obj.(fdName));                % Fields consisting of single matrix or cell, or zeros(1,1); 
                isMultiField  =  isstruct(obj.(fdName))&& ~isPropField; % Fields consisting of a structure of matrices and/or cells
                
                if     isPropField;   subList.(fdName) = obj.(fdName);
                elseif isEmptyField;  subList.(fdName) = [];
                elseif isSingleField; nentries = numel(obj.(fdName)(:,1));
                                      if nentries==ntr; subList.(fdName) = obj.(fdName)(selectIndex,:);
                                      else              subList.(fdName) = 0;
                                      end
                    
                elseif isMultiField
                    
                    % Subfields Level 1 ...................................
                    sfNames = fieldnames(obj.(fdName));
                    nsf     = numel(sfNames);
                    for isf = 1:nsf
                        
                        sfName   = sfNames{isf};
                        nentries = numel(obj.(fdName).(sfName));
                        
                        % Determine field type
                        isSingleField = ~isstruct(obj.(fdName).(sfName));
                        isMultiField  =  isstruct(obj.(fdName).(sfName));
                        isEmptyField  =  isempty( obj.(fdName).(sfName));
                        
                        if isEmptyField;  subList.(fdName).(sfName) = [];
                        
                        elseif isSingleField;
                            
                            if nentries==ntr; subList.(fdName).(sfName) = obj.(fdName).(sfName)(selectIndex,:);
                            else              subList.(fdName).(sfName) = 0;
                            end
                            
                        elseif isMultiField
                            
                            % Subfields Level 2 ...........................
                            % Only accept single fields at this level
                            ssfNames = fieldnames(obj.(fdName).(sfName));
                            nssf     = numel(ssfNames);
                            for issf = 1:nssf
                                ssfName  = ssfNames{issf};
                                nentries = numel(obj.(fdName).(sfName).(ssfName));
                                if nentries==ntr; subList.(fdName).(sfName).(ssfName) = obj.(fdName).(sfName).(ssfName)(selectIndex,:);
                                else              subList.(fdName).(sfName).(ssfName) = 0;
                                end
                            end
                        end
                    end % .......................................... end L1
                end
            end % ........................................... end main loop
        end % ................................................ end function
        
        
        
        
        % -----------------------------------------------------------------
        % Append lists to existing one   -  -  -  -  -  -  -  -  -  -  -  -  
        function appendList(obj,NewList)
            names = fieldnames(obj);
            nf    = size(names,1);
            ntr   = numel(NewList.(names{1}));
            
            for ifd = 1:nf
                
                fdName = names{ifd};
                
                % Determine field type
                isPropField   = strcmp(fdName,'prop');                  % Field with one structure entry per list, not per trace
                isEmptyField  = isempty(NewList.(fdName));
                isSingleField = ~isstruct(NewList.(fdName));                % Fields consisting of single matrix or cell, or zeros(1,1);
                isMultiField  =  isstruct(NewList.(fdName))&& ~isPropField; % Fields consisting of a structure of matrices and/or cells
                
                if isPropField &isempty(fieldnames(obj.prop)); 
                    obj.(fdName) = NewList.(fdName); 
                end %Overwrite prop-field with contents of new list?
                %if isEmptyField;  obj.(fdName) = obj.(fdName);
                if isSingleField; nentries = numel(NewList.(names{1}));
                    if nentries==ntr; obj.(fdName) = [obj.(fdName); NewList.(fdName)]; 
                    %else              obj.(fdName) = 0;
                    end
                
                elseif isMultiField
                    
                    % Subfields Level 1 ...................................
                    sfNames = fieldnames(obj.(fdName));
                    nsf     = numel(sfNames);
                    for isf = 1:nsf
                        
                        sfName   = sfNames{isf};
                        nentries = numel(NewList.(fdName).(sfName));
                        
                        % Determine field type
                        isSingleField = ~isstruct(obj.(fdName).(sfName));
                        isMultiField  =  isstruct(obj.(fdName).(sfName));
                        
                        if isSingleField;
                            if nentries==ntr; obj.(fdName).(sfName) = [obj.(fdName).(sfName); NewList.(fdName).(sfName)]; 
                            end
                            
                        elseif isMultiField
                            
                            % Subfields Level 2 ...........................
                            % Only accept single fields at this level
                            ssfNames = fieldnames(obj.(fdName).(sfName));
                            nssf     = numel(ssfNames);
                            for issf = 1:nssf
                                ssfName  = ssfNames{issf};
                                nentries = numel(NewList.(fdName).(sfName).(ssfName));
                                if nentries==ntr; obj.(fdName).(sfName).(ssfName) = [obj.(fdName).(sfName).(ssfName); NewList.(fdName).(sfName).(ssfName)]; 
                                end
                            end
                        end
                    end % .......................................... end L1
                end
            end % ........................................... end main loop
        end        

        
        
        % Sort lists  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
        function sortList(obj,sortIndex)
            names = fieldnames(obj);
            nf    = size(names,1);
            ntr   = numel(obj.(names{1}));
            
            for ifd = 1:nf
                
                fdName   = names{ifd};
                nentries = numel(obj.(fdName));
                
                % Determine field type
                isPropField   = strcmp(fdName,'prop');                  % Field with one structure entry per list, not per trace
                isSingleField = ~isstruct(obj.(fdName));                % Fields consisting of single matrix or cell
                isMultiField  =  isstruct(obj.(fdName))&& ~isPropField; % Fields consisting of a structure of matrices and/or cells
                
                if isSingleField
                    if nentries==ntr; obj.(fdName) = obj.(fdName)(sortIndex);
                    else              obj.(fdName) = 0;     % Overwritten fields
                    end
                    
                elseif isMultiField

                    % Subfields
                    sfNames = fieldnames(obj.(fdName));
                    nsf     = numel(sfNames);
                    for isf = 1:nsf
                        sfName   = sfNames{isf};
                        nentries = numel(obj.(fdName).(sfName));
                        
                        % Determine field type
                        isSingleField = ~isstruct(obj.(fdName).(sfName));
                        isMultiField  =  isstruct(obj.(fdName).(sfName));
                        
                        if isSingleField;
                            if nentries==ntr; obj.(fdName).(sfName) = obj.(fdName).(sfName)(sortIndex); 
                            end
                            
                        elseif isMultiField
                            
                            % Subfields Level 2 ...........................
                            % Only accept single fields at this level
                            ssfNames = fieldnames(obj.(fdName).(sfName));
                            nssf     = numel(ssfNames);
                            for issf = 1:nssf
                                ssfName  = ssfNames{issf};
                                nentries = numel(obj.(fdName).(sfName).(ssfName));
                                if nentries==ntr; obj.(fdName).(sfName).(ssfName) = obj.(fdName).(sfName).(ssfName)(sortIndex); 
                                end
                            end
                        end
                        
                            %                         %ORG
                            %                         fieldFullName = sprintf('%s.%s',fdName,sfName);
                            %                         evalCmd       = sprintf('obj.%s = obj.%s(sortIndex);',fieldFullName,fieldFullName);
                            %                         eval(evalCmd)
                    end
                elseif isPropField
                    % do nothing
                    1+1;
                else
                    error('Never should have arrived here. Check.')
                end

                %                 nentries = numel(obj.(names{ifd}));
                %                 ntr      = numel(sortIndex);
                %                 if ntr==nentries 
                %                     obj.(names{ifd}) = obj.(names{ifd})(sortIndex);
                %                 end
            end
        end
        
        
        
        
        % -----------------------------------------------------------------
        % Clear out empty entries
        function consolidate(obj)
            idx_lastEntry = find(cellfun(@(x) ~isempty(x), obj.fullName),1,'last');
            names         = fieldnames(obj);
            ntr           = numel(obj.(names{1}));
            if idx_lastEntry~=ntr;
                fprintf(1,'Need to update consolidate function. Turns out it is not obsolete...\n')
                pause
            end
            %             nf            = size(names,1);
            %             for ifd = 1:nf
            %                 
            %                 fdName = names{ifd};
            %                 nentries = numel(obj.(fdName)(:,1));
            %                 
            %                 % Determine field type
            %                 isPropField   = strcmp(fdName,'prop');                  % Field with one structure entry per list, not per trace
            %                 isSingleField = ~isstruct(obj.(fdName));                % Fields consisting of single matrix or cell
            %                 isMultiField  =  isstruct(obj.(fdName))&& ~isPropField; % Fields consisting of a structure of matrices and/or cells
            %                 
            %                 if     isPropField;   obj.(fdName) = obj.(fdName);
            %                     
            %                 elseif isSingleField;
            %                     if nentries==ntr; obj.(fdName) = obj.(fdName)(1:idx_lastEntry);
            %                     else              obj.(fdName) = 0;     % Overwritten fields
            %                     end
            %                     
            %                 elseif isMultiField
            %                     
            %                     % Subfields
            %                     sfNames = fieldnames(obj.(fdName));
            %                     nsf     = numel(sfNames);
            %                     for isf = 1:nsf
            %                         sfName       = sfNames{isf};
            %                         fieldFullName = sprintf('obj.%s.%s',fdName,sfName);
            %                         evalCmd       = sprintf('%s = %s(1:idx_lastEntry);',fieldFullName,fieldFullName);
            %                         eval(evalCmd)
            %                     end
            %                 else                           
            %                     fprintf(1,'Warning: untested object function for prop-field\n')
            %                 end
            %             end
        end
        
        
        %         % Remove skipped traces from traceList  -  -  -  -  -  -  -  -  -  
        %         % COULDNT THIS BE DONE WITH THE SKIPLIST FUNCTION?!
        %         function removeSkipped(obj,flg_skipTrace)
        %             
        %             nflgVals = numel(flg_skipTrace);
        %             names    = fieldnames(obj);
        %             nf       = size(names,1);
        %             ntr      = numel(obj.(names{1}));
        %             %ntr  = numel(obj.m);     
        %             if nflgVals~=ntr 
        %                 error('Use logicals/flags for specifing traces to skip, not indices\n')
        %             end
        %             
        %             dontSkip = (flg_skipTrace==0);
        %             
        %             
        %             for ifd = 1:nf
        %                 
        %                 fdName   = names{ifd};
        %                 nentries = numel(obj.(names{ifd}));
        %                 
        %                 % Determine field type
        %                 isPropField   = strcmp(fdName,'prop');                  % Field with one structure entry per list, not per trace
        %                 isSingleField = ~isstruct(obj.(fdName));                % Fields consisting of single matrix or cell
        %                 isMultiField  =  isstruct(obj.(fdName))&& ~isPropField; % Fields consisting of a structure of matrices and/or cells
        %                 
        %                 if isPropField
        %                     obj.(names{ifd}) = obj.(names{ifd});
        %                     
        %                 elseif isSingleField
        %                     if nentries==ntr; obj.(fdName) = obj.(fdName)(dontSkip);
        %                     else              obj.(fdName) = 0;     % Overwritten fields
        %                     end
        %                     
        %                 elseif isMultiField
        % 
        %                     % Subfields
        %                     sfNames = fieldnames(obj.(fdName));
        %                     nsf     = numel(sfNames);
        %                     for isf = 1:nsf
        %                         sfName       = sfNames{isf};
        %                         fieldFullName = sprintf('obj.%s.%s',fdName,sfName);
        %                         evalCmd       = sprintf('%s = %s(dontSkip);',fieldFullName,fieldFullName);
        %                         eval(evalCmd)
        %                     end
        %                 else
        %                     error('Never should have arrived here. Check.')
        %                 end
        %                 
        %                 %                 nentries = numel(obj.(names{ifd}));
        %                 %                 if nentries==ntr; obj.(names{ifd}) = obj.(names{ifd})(dontSkip);
        %                 %                 elseif            strcmp(names{ifd},'prop'); obj.(names{ifd}) = obj.(names{ifd});
        %                 %                 else              fprintf(1,sprintf('Note: Field %s only has %i entries; not removing any entries in this field\n',names{ifd},nentries))
        %                 %                 end
        %             end
        %         end

        
        
        
        % Add single or several entries     -  -  -  -  -  -  -  -  -  -  -  
        function addEntries(obj,sourceList,traceIdx)
            fprintf(1,'Warning: object function not updated to hanlde sub-fields and prop-field\n')
            idx_nextline = find(cellfun(@(x) isempty(x), obj.fullName),1,'first');
            nadds = numel(traceIdx);
            names = fieldnames(obj);
            nf    = size(names,1);
            for ifd = 1:nf
                if ~strcmp(names{ifd},'prop'); obj.(names{ifd})(idx_nextline:idx_nextline+nadds-1) = sourceList.(names{ifd})(traceIdx);
                else                           obj.(names{ifd}) = obj.(names{ifd});
                                               fprintf(1,'Warning: untested object function for prop-field\n')
                end
            end
        end
        
        % Add minimum info of single entry which is to be skipped  -  -  -
        function addSkipEntry(obj,fullName,reason,hypDist,m)
            fprintf(1,'Warning: object function not updated to hanlde sub-fields and prop-field\n')
            idx_nextline = find(cellfun(@(x) isempty(x), obj.fullName),1,'first');
            obj.fullName{idx_nextline}   = fullName;
            obj.comment{idx_nextline}    = reason;
            obj.hypDist(idx_nextline)    = hypDist;
            obj.m(idx_nextline)          = m;
        end
        
        % Return independent clone of the list
        % use obj.selectSublist for this purpose
        
        % Return a clone of the list, but containing only the "fullName" property
        function [nameList] = cloneNameList(obj)
            nameList          = traceNameList(size(obj.fullName,1));	%generate new object
            nameList.fullName = obj.fullName;
            nameList.ppxIdx   = obj.ppxIdx;
            nameList.tppx     = obj.tppx;
        end

        % Print most important info 
        function printSingleTraceSummary(obj,idx)
            fprintf(1,['\n\nFile name:\t\t\t',obj.fullName{idx},'\n'])
            fprintf(1,['EQ-coords (lat/lon/z):\t\t',num2str([obj.eq.lat(idx),obj.eq.lon(idx),obj.eq.z(idx)]),'km\n'])
            fprintf(1,['Magnitude:\t\t\t',num2str(obj.eq.m(idx)),'\t(of type ',obj.eq.mType{idx},')\n'])
            fprintf(1,['Hypocentral Distance:\t\t',num2str(obj.dist.hyp(idx),'%4.1f'),'km\n'])
            fprintf(1,['Epicentral  Distance:\t\t',num2str(obj.dist.epi(idx),'%4.1f'),'km\n'])
            fprintf(1,['Station-coords (lat/lon/alt):\t',num2str([obj.station.lat(idx),obj.station.lon(idx),obj.station.alt(idx)]),'m\n'])
            fprintf(1,['Station-specifications:\t\tband code ',obj.station.bcode{idx},'orientation code ',obj.station.ocode{idx},'\tinstrument code ',obj.station.icode{idx},'\n'])
            %fprintf(1,['Comments:\t\t\t',obj.comment{idx},'\n'])
            fprintf(1,['Comments:\t\t\tdeactivated. fix it. do it.\n'])
            fprintf(1,['-------------------------------------------------------------------------\n'])
        end

        % Print single line summary, e.g. for titles
        function printString = printSingleLineSummary(obj,idx)
            printString = [];
            if obj.eq.m(idx)~=0;       printString = sprintf('%sM%3.1f',printString,obj.eq.m(idx)); end
            if obj.dist.hyp(idx)~=0;   printString = sprintf('%s, R%5.1fkm',printString,obj.dist.hyp(idx)); end
            printString = sprintf('%s, channel %s%s%s',printString,obj.station.bcode{idx},obj.station.icode{idx},obj.station.ocode{idx});
            printString = sprintf('%s from %s',printString,obj.dataSetName{idx});
            printString = sprintf('%s, %s',printString,obj.eq.date{idx});
        end
        
        % Print single line summary, e.g. for filenames
        function printString = printSingleLineMiniSummary(obj,idx)
            printString = strrep(sprintf('M%3.1f_r%5.1fkm',obj.m(idx),obj.hypDist(idx)),'.','p');
            printString = strrep(printString,' ','');
        end
        
        %         % Update relative indices  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        %         function updateRelIdx(obj)
        %             fprintf(1,'Not yet written. do it.\n')
        %         end
        
        
        % Compute and print size of traceList  -  -  -  -  -  -  -  -  -  -
        function [totSize] = printObjectSize(obj)
            props   = properties(obj);
            totSize = 0;
            
            for ii=1:length(props)
                currentProperty = getfield(obj, char(props(ii)));
                s               = whos('currentProperty');
                totSize         = totSize + s.bytes;
            end
            
            if (totSize>1e9); fprintf(1, '%3.1f GB\n', totSize*1e-9)
            else              fprintf(1, '%3.1f MB\n', totSize*1e-6)
            end
        end
        
        
    end
end
