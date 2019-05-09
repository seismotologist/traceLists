% HOW TO USE TRACELIST OBJECTS
% 
% TraceLists are a way to organise seismic time series data in matlab.  
% Each record ('trace') has a list entry, with a large number of attributes 
% such as  event information (hypocenter, magnitude, etc.), station 
% information (station coordinates, SNCL code, etc.), peak ground motion
% statistics (PGA, PGV, etc.). The traceList object contains all relevant
% information for each record, and allows to easily manage the record data
% set. In particular, it is easy to extract subsets of the full data set,
% without having to loop over all attributes. This quick tutorial shows how 
% traceLists can be used to manage time series data sets efficiently.
%
% traceLists are handle class objects. To be able to load a traceList 
% the directory that contains the object definition file (traceList.m) 
% needs to be in the search path (use e.g. 
% "addpath(genpath('../path-to-traceList.m/'))"
%
% A mild warning: this is a piece of research code and contains some 
% ad-hoc, somewhat hacky solutions. But it works very well and I find it a 
% very useful tool for manipulating time series data; it greatly simplifies 
% scripts and processes. 
%
% mmeier@caltech.edu, last update 190509                              (\_/)


clear all
addpath(genpath('fct/'))

% Example traceList, populated with processed records from M>6 earthquakes
% from southern California
load data/trList.mat

% The traceList is called 'exampleList' & contains a number of main fields:
exampleList

% . fullName        % Waveform file name including absolute path
% . eq              % Event information (hypocenter, origin time, etc.))
% . dist            % Various recording distance measures
% . station         % Recording station information
% . dataSetName     % Name for identifying data sets of different origins
% . px              % Pick times and indices
% . fb              % Filter bank information
% . pga             % Peak ground acceleration information
% . pgv             % Peak ground velocity information
% . pgd             % Peak ground displacement information
% . mmi             % Modified Mercalli Intensity
% . noise           % Pre-signal noise statistics
% . comment         % e.g. skipReasons for skipped traces
% . var             % Placeholder for storing unforseen quantities
% . scalFeature     % Scalar waveform features
% . vectFeature     % Vector waveform features
% . prop            % Info that is not trace-specific, such as hyper params

% Each field is either a matrix or cell with ntr entries (where ntr is the 
% number of records in the list), or a nested structure containing such 
% matrices and cells. For example, the <eq> field contains various entries
% with event information:
exampleList.eq

% Plot epicenter map
clf; hold on; grid on; box on;
plot_Cali_map
plot(exampleList.eq.lon,exampleList.eq.lat,'pk','markerFaceColor','y','markerSize',15)

clf; histogram(exampleList.eq.m)
xlabel('Magnitude')
ylabel('No. of records in exampleList')

% The total number of records in list:
ntr = numel(exampleList.eq.m)



%% The property field
% The 'prop' field is the only field that does not contain one entry per 
% record. I use it to store data set information that is the same for all 
% records, for example names of files that were used in the waveform
% processing, or processing parameters
exampleList.prop
exampleList.prop.files      % File names
exampleList.prop.fc         % Filter bank corner frequencies



%% Select a subList
%  Let's extract separate traceLists for each sensor component (vertical, 
%  east, north):
fprintf(1,'\nSplitting exampleList into z-, e- & nList ... ')
idx_Z = find(strcmp(exampleList.station.ocode,'Z'));
idx_E = find(strcmp(exampleList.station.ocode,'E'));
idx_N = find(strcmp(exampleList.station.ocode,'N'));
zList = exampleList.selectSubList(idx_Z);
eList = exampleList.selectSubList(idx_E);
nList = exampleList.selectSubList(idx_N);
%clear exampleList
 
nz = numel(zList.eq.m)
ne = numel(eList.eq.m)
nn = numel(nList.eq.m)
fprintf(1,' done.\n')


% Extract vertical records with hypocentral distances <100km
idx    = find(zList.dist.hyp<100);
nsList = zList.selectSubList(idx);

% Extract vertical records with hypocentral distances <100km & magnitudes>6.5  
tmpList = zList.selectSubList(find(zList.dist.hyp<100 &zList.eq.m>6.5));

% You can also use the selectSubList function e.g. to sort the list in any
% order you prefer, e.g. with increasing hypocentral distances:
[val,idx]   = sort(zList.dist.hyp);
zsortedList = zList.selectSubList(idx);

clf; hold on;
plot(zList.dist.hyp,'k')                       % Original list, not sorted
plot(zsortedList.dist.hyp,'r','lineWidth',2)   % Sorted list
ylabel('Hypocentral Distance [km]')


% --> this selection function is one of the most useful features of the 
%     traceList format. You make a selection of records and you 
%     automatically retain all attributes for those records, without having 
%     to loop over various fields. I think that was why I started the 
%     format it in the first place :)




%% Various commands
%  Print main info on 51th trace
itr=51; % --> traceList index of randomly chosen record

% The file name of the waveform file including path is stored in the 
% fullName field
zList.fullName{itr}

% You can print a summary of any particular trace
zList.printSingleTraceSummary(itr);

% In the exampleList, the three compnents of each record were listed next
% to each other. After splitting the list, the three sublists have the same 
% order, i.e. you can use the same traceList index to find the three 
% components of the same record:
zList.printSingleTraceSummary(itr);
eList.printSingleTraceSummary(itr);
nList.printSingleTraceSummary(itr);

% The size of the object
zList.printObjectSize;





%% Find records of individual earthquakes
% Use the following function to find all records of any individual
% earthquake:
eqs = compile_eqs_structure(exampleList)

% The traceList-indices are stored in the field eqs.traceId, so you can
% make a new sublist with all records from the same earthquake:
ieq        = 6; % 1994, Mw6.7 Northridge California
eqname     = eqs.name{ieq};
tridx      = eqs.traceId{ieq};
nridgeList = exampleList.selectSubList(tridx);

clf; hold on; grid on; box on;
plot_Cali_map
plot(nridgeList.eq.lon     ,nridgeList.eq.lat     ,'pk','markerFaceColor','r','markerSize',15)
plot(nridgeList.station.lon,nridgeList.station.lat,'vk','markerFaceColor','y','markerSize',11)
set(gca,'xlim',[-122 -115],'ylim',[32 37])




%% Combine lists
%  e.g. make list of all vertical records of 1994 Mw6.7 Northridge (event 
% name 3144585) and 1992 Mw7.2 Landers (event name 3031111):
idx1            = find(cellfun(@(x) ~isempty(x), regexp(zList.fullName,'3144585')));
firstEventList  = zList.selectSubList(idx1);
idx2            = find(cellfun(@(x) ~isempty(x), regexp(zList.fullName,'3031111')));
secondEventList = zList.selectSubList(idx2);

bothEventsList = traceList(0);                  % Intitiate traceList
bothEventsList.appendList(firstEventList);      % Append first list
bothEventsList.appendList(secondEventList);     % Append second list




%% Clone list (IMPORTANT!)
%  Because traceLists are "handle classes", copying them with a command like
%  newzList = zList only creates a new reference to the same object, not a
%  new object itself. If you then make changes to newzList the same changes
%  will also be applied to zList itself. To create a new independent
%  instance of a traceList, "clone" the list using the selectSubList function:
nz                  = numel(zList.eq.m);
newIndependentzList = zList.selectSubList(1:nz);



%% Remove unneeded fields to save space
%  If the number of records is large, the traceLists can get very heavy. 
%  In this case you want to remove fields that you don't need to save 
%  space. This can be done using the following function, with a list of
%  fields you would like to discard:
overwrite_unneeded_fields(zList,{'var.v7';'fb.cav';'fb.amax'})

% Although traceLists are in no way optimised for efficient information
% storage, you can use them for manipulating up to millions of records at 
% once if you delete unnecessary fields.



%% Flexible fields
%  You can use the var.v1-8, as well as the scalFeature & vectFeature
%  fields to store all sorts of numeric data in a flexible way. I use these
%  fields to store things I hadn't thought of when I wrote the object
%  definition file. For example, I would use the 'acczen' subfield to store
%  quantities computed on the vector sum of all three acceleration 
%  components of a record. By convention I store all three component 
%  information with the  vertical component  
zList.scalFeature{33}.acczen.zhr

% You can then use cellfun to extract the information from all records in a
% list
zhr = cell2mat(cellfun(@(x) x.acczen.zhr', zList.scalFeature,'uniformOutput',0));



%% A final note
%  You can populate the traceList fields any way you like. The examples
%  here only show how I use them. If you're interested in the waveform 
%  processing scripts that I used to populate the exampleList, please 
%  send me an email. Also, I have an entire library of functions for
%  working with traceLists. If you send me a short description of what you
%  would like to do I can probably help out.


% If you've made it that far, olé! For comments and suggestions please send
% me an email, cheers.
%
%
%                             I liked ascii-tables better...
%             ,;::\::\        
%           ,'/' `/'`/       /
%       _\,: '.,-'.-':.     /
%      -./"'  :    :  :\/, /
%       ::.  ,:____;__; :-       
%       :"  ( .`-*'o*',);       
%        \.. ` `---'`' /
%         `:._..-   _.'
%         ,;  .     `.
%        /"'| |       \
%       ::. ) :        :
%       |" (   \       |
%       :.(_,  :       ;
%        \'`-'_/      /
%         `...   , _,'
%          |,|  : |
%          |`|  | |
%          |,|  | |
%      ,--.;`|  | '..--.
%     /;' "' ;  '..--. ))
%     \:.___(___   ) ))'                                              (\_/)