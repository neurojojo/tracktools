classdef tracksTableClass < handle

%
% Title: tracksTableClass
% Author: Jozsef Meszaros
%
% License: Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
% (CC BY-NC-SA 4.0)
%
% Creative Commons License
%
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
% To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
%
%--------------------------------------------------------------------------
% Usage:
%--------------------------------------------------------------------------
% output_table = tracksTableClass('file.mat');
%

    properties
       tracksTable; % [Public] TABLE with four columns (id, trackStart, x, y)
       metadata; % <b>[Public]</b> STRUCTURE with fields:
                 % obj_idx (double), 
                 % Directory (str), 
                 % fileStruct (loaded files), 
                 % Type (str), 
                 % Comments (str), 
                 % FilteredTracksRatio
                 % trackBoundaries, 2-D ARRAY containing x,y of the cell's boundary
                 % boundaryStats, TABLE with the Area and Perimeter as columns
    end
    
    properties(GetAccess='private')
       Ntracks; % [Private] Contains the number of tracks in the table before any filtering happens 
       dims=[0,256,0,256]; % [Private] This could change but for our data it seems constant
    end
    
    properties(GetAccess='private')
        defaultFiles = {'Tracking.mat','results.mat'};
        % These are properties for the table (currently 'id', 'trackStart', 'x', 'y')
        varnames = {'id','trackStart','x','y'};
        % These are the variable types for the table (currently double, double, cell, cell)
        vartypes = {'uint8','uint8','cell','cell'};
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % CONSTRUCTOR FUNCTION %
        %%%%%%%%%%%%%%%%%%%%%%%%
        function obj = tracksTableClass(trackingFile,varargin)
            
            % Check for errors
            assert(ischar(trackingFile),'Input must be a string that is the location of a Tracking.mat file');
            assert(exist(trackingFile)>0,'The file is not found at this location');
            assert(isstruct(varargin{1}),sprintf('Pass a structure as the second argument\n\nie. struct("height",256,"width",256)'));
            
            % Check if a structure with options is provided
            if nargin==2; options = varargin{1}; obj.dims = [0,options.height,0,options.width]; else options = struct(); end
            
            % Load the file, save metadata
            tracksFinal = load( trackingFile );
            trackingFile = split( trackingFile, '\' );
            obj.metadata.TLD = fullfile(trackingFile{1:end-1}); % Preserve the TLD for this folder, if provided
            fprintf('\nTracking file successfully loaded from within directory:\n\n%s\n',obj.metadata.TLD);
            
            % Create the table
            Ntracks = max( size(tracksFinal.tracksFinal) );
            %obj.tracksTable = table('Size',[Ntracks numel(obj.varnames)],'VariableTypes',obj.vartypes,'VariableNames',obj.varnames);
            %obj.tracksTable.id = [1:Ntracks]';
            obj.tracksTable = ...
                    obj.extractPosition(obj, tracksFinal.tracksFinal);
            
            obj.metadata.Type = 'tracks';
            obj.metadata.Comments='NA';
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % END CONSTRUCTOR FUNCTION %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Information functions %
        
        function analyzeTracks( obj )
            
            [x_end,y_end] = deal( rowfun( @(x) x{1}(end), obj.tracksTable(:,{'x'}), 'OutputFormat', 'uniform' ),...
                rowfun( @(y) y{1}(end), obj.tracksTable(:,{'y'}), 'OutputFormat', 'uniform' ) );
            [x_start,y_start] = deal( rowfun( @(x) x{1}(1), obj.tracksTable(:,{'x'}), 'OutputFormat', 'uniform' ),...
                rowfun( @(y) y{1}(1), obj.tracksTable(:,{'y'}), 'OutputFormat', 'uniform' ) );           
            
            D = pdist2( [x_start, y_start ], [x_end, y_end] ); % A distance matrix 
            D_lt_2 = find(D>0 & D<10);
            
            for i = 1:size(D_lt_2,1); [particle1(i),particle2(i)] = ind2sub( size(D), D_lt_2(i) ); end;
            
            particle_timeDiff = arrayfun( @(x) (obj.tracksTable( particle1(x),: ).trackLength + ... 
                obj.tracksTable( particle1(x),: ).trackStart-obj.tracksTable( particle2(x),: ).trackStart), [1:numel(particle1)] );
            particle_d1d2 = table( particle1', particle2',D(D_lt_2),...
                'VariableNames', {'Particle1','Particle2','timeDiff','distanceDiff'} );
            
            figure; imagesc( D );
        end
        
        function getTrackBoundaries(obj)
% Makes a mask and then retrieves the edges as obj.trackBoundaries

% (1) Places a 1 in EVERY x,y loc that any track touches at any time
% (2) Discards any 1's that are further than 1 pixel away from
% another 1
% (3) Calculates the location of any holes in where trajectories go
% (4) Obtains hole information

            fprintf('Currently on %i\n',obj.metadata.obj_idx);
            if gt(numel( obj.tracksTable ),0)
                
                N_tracks = size( obj.tracksTable, 1 );
                
                % Obtain mask over all tracks and image holes %
                unfilled_matrix = zeros(300,300);
                unfilled_matrix( sub2ind([300,300],max(1,fix(cell2mat(obj.tracksTable.x')')),max(1,fix(cell2mat(obj.tracksTable.y')'))) ) = 1;
                unfilled_matrix = arrayfun(@(x) lt( x, 1), bwdist( unfilled_matrix ));
                matrix = imfill(unfilled_matrix,'holes');
                image_holes =  matrix - unfilled_matrix;
                % Plot result 
                % figure; imagesc(matrix + image_holes);
                
                
                
                % Obtain boundaries of holes %
                % EVERY BOUNDARY CELL FROM EVERY HOLE ( N holes x 1 )
                all_hole_B = bwboundaries( image_holes );
                % EVERY PIXEL FROM EVERY HOLE ( k pixels x 2, within size of image )
                holeboundaries = cell2mat( all_hole_B );
                % EVERY PIXEL'S RESPECTIVE HOLE ( k pixels x 1, range [1:N holes] )
                hole_idx = cell2mat(arrayfun( @(x) x*ones(1,numel(all_hole_B{x})), [1:numel(all_hole_B)] , 'UniformOutput', false))';
                
                my_holes_pixeltable = table( 'Size',[ size(holeboundaries,1), 3 ],...
                    'VariableType',{'single','single','single'},...
                    'VariableNames',{'Pixel_idx','Hole','Disappearances'} );
                my_holes_pixeltable.Pixel_idx = arrayfun( @(i) sub2ind(size(matrix), holeboundaries(i,1),holeboundaries(i,2) ), [1:size(holeboundaries,1)] )';
                my_holes_pixeltable.Hole = cell2mat(arrayfun( @(x) x*ones( size(all_hole_B{x},1), 1 ), [1:numel(all_hole_B)] , 'UniformOutput', false)');
                % End of getting boundaries for holes table %
                
                % Track endings %
                % Get the x,y position of the track endings %
                [x_end,y_end] = deal( rowfun( @(x) x{1}(end), obj.tracksTable(:,{'x'}), 'OutputFormat', 'uniform' ),...
                    rowfun( @(y) y{1}(end), obj.tracksTable(:,{'y'}), 'OutputFormat', 'uniform' ) );
                [x_start,y_start] = deal( rowfun( @(x) x{1}(1), obj.tracksTable(:,{'x'}), 'OutputFormat', 'uniform' ),...
                    rowfun( @(y) y{1}(1), obj.tracksTable(:,{'y'}), 'OutputFormat', 'uniform' ) );
                % End Track endings %
                
                % Get boundary of the cell %
                [Boundary,~] = bwboundaries( matrix );
                [~,b] = max( cellfun(@(x) size(x,1), Boundary ) ); % Get only the largest boundary
                cellboundary = Boundary{b};
                obj.metadata.trackBoundaries = cellboundary;
                % End of getting boundary %
                
                % Distances between HOLES/EDGES and TRACKS %
                boundaries_to_seek = [cellboundary; holeboundaries]; % Test the track trajectories against the boundaries of the cell and holes
                Nboundary1 = size( cellboundary, 1 );
                Nboundary2 = size( holeboundaries, 1 );
                
                [ pdists_start, pdists_end ] = deal( pdist2( boundaries_to_seek, [x_start,y_start] ),...
                    pdist2( boundaries_to_seek, [x_end,y_end] ) );
                
                [ min_values_starts, min_value_locs_starts ] = min( pdists_start, [], 1 );
                [ min_values_ends, min_value_locs_ends ] = min( pdists_end, [], 1 );
                
                % The vector min_value_locs tells you into which of the
                % holes the track may have disappeared
                %
                % Condition them like so:
                d = 3;
                
                table_pdists = table( [1:N_tracks]', ...
                    min_values_starts', ...
                    min_value_locs_starts',...
                    min_values_ends',...
                    min_value_locs_ends',...
                    repmat( {'Central'}, N_tracks, 1 ),...
                    repmat( {'Central'}, N_tracks, 1 ),...
                    repmat( 0, N_tracks, 1 ),...
                'VariableNames', {'trackIdx',...
                            'min_values_starts',...
                            'min_value_locs_starts',...
                            'min_values_ends',...
                            'min_value_locs_ends',...
                            'start_loc','end_loc','keep_logical'} );
                
                %min_value_locs_starts'<=Nboundary1,...
                %min_value_locs_ends'<=Nboundary1,...
                %table_pdists( table_pdists.min_value_locs_starts'<=Nboundary1, : ).Start_Loc = 'Cell';
                
                % Custom for boundary types %
                crit1 = and( (min_values_starts<=d)' , (min_value_locs_starts'<=Nboundary1) );
                crit2 = and( (min_values_ends<=d)' , (min_value_locs_ends'<=Nboundary1) );
                crit3 = and( (min_values_ends>d)' , (min_values_starts>d)' );
                
                table_pdists( crit1 , :).start_loc = repmat({'Cell'},sum(crit1),1);
                table_pdists( crit2 , :).end_loc =  repmat({'Cell'},sum(crit2),1);
                table_pdists( crit3 , :).keep_logical =  repmat(1,sum(crit3),1);
                
                obj.tracksTable.NonEdgeTrack = table_pdists.keep_logical;
                
            end
        end
        
        function getBoundaryStats(obj)
% Retrieves the area and perimeter of a masked version of the track boundaries, stores the rest as a table in the boundaryStats property
            makemask = @(x) poly2mask( x(:,1), x(:,2), 256, 256 );
            area_perimeter = regionprops( makemask( obj.trackBoundaries ), 'Area', 'Perimeter');
            mytable = struct2table( area_perimeter ); 
            obj.metadata.boundaryStats = table(sum( mytable.Area ), sum( mytable.Perimeter ),'VariableNames',{'Area','Perimeter'});
        end
        
% Overloaded functions %
            
        function hits = find(obj,search_query)
% Function that takes an argument, which is an inequality based on some column of the table. Demo: findFoldersObj.tracks.obj_1.find('trackLength>1000')            
            equality_type = cell2mat(regexp(search_query,'(<=)|(>=)|(<>)|(==)|(<)|(>)','match'));
            value = str2double(regexp( search_query, '(\d+)$', 'match' ));
            parameter = regexp( search_query, '[a-zA-Z]+', 'match' ); parameter = parameter{1};
            
            switch equality_type
                case '<'
                    hits = find( lt(obj.tracksTable.(parameter), value) == 1 );
                case '>'
                    hits = find( gt(obj.tracksTable.(parameter), value) == 1 );
                case '<='
                    hits = find( le(obj.tracksTable.(parameter), value) == 1 );
                case '>='
                    hits = find( ge(obj.tracksTable.(parameter), value) == 1 );
                case '<>'
                    hits = find( ne(obj.tracksTable.(parameter), value) == 1 );
                case '=='
                    hits = find( eq(obj.tracksTable.(parameter), value) == 1 );
            end
            
        end
        
        function hist(obj,varargin)

            fig = figure('color','w');  
            parameter = 'trackLength';
            charargin = varargin( find(cellfun(@(x) isa(x,'char'), varargin )==1) );
            
            if any(strcmp(charargin,'NonEdgeTrack'))
                trackIdx = obj.find('NonEdgeTrack==1');
            end
            
            if any(strcmp(charargin,'EdgeTrack'))
                trackIdx = obj.find('NonEdgeTrack==0');
            end
            
            if any(strcmp(charargin,'trackLength'))
                parameter = 'trackLength';
            end
            
            if strcmp(parameter,'trackLength'); hist( obj.tracksTable(trackIdx,:).(parameter) ); end
            
        end
        
        function plot(obj,varargin)
% A function to plot 

           fig = figure('color','w');
           % Custom input parsing
           
           % (1) Check for any double arrays (which will be the indices of the tracks to be plotted)
           %
           % Any array present should be the indices of the tracks desired
           % to be plotted
           if any(cellfun(@(x) isa(x,'double'), varargin))
               trackIdx = varargin{ find(cellfun(@(x) isa(x,'double'), varargin)==1) };
           end
           
           % (2) Check for character arrays as inputs
           %
           % Possible arguments:
           % - bounds
           % - trackIdx
           % - color
           % - save
           % - NonEdgeTrack
           charargin = varargin( find(cellfun(@(x) isa(x,'char'), varargin )==1) );
           
           if any(cell2mat(strfind(charargin,'bounds')))
               linecolor = 'k';
               if any(cell2mat(strfind(charargin,'color')))
                  linecolor = 'w';
               end
               patch( obj.trackBoundaries(:,1), obj.trackBoundaries(:,2), linecolor, 'facealpha', 0.1 );
           end
           
           if ~exist('trackIdx'); % trackIdx hasn't been specified so give a random sets
               trackIdx = unique( obj.tracksTable.id(obj.tracksTable.trackLength>20) );
               if numel(trackIdx)>2000
                  trackIdx = trackIdx( randperm(numel(trackIdx),2000) );
               end
           end
           
           if any(strcmp(charargin,'NonEdgeTrack'))
                trackIdx = obj.find('NonEdgeTrack==1');
           end
           
           if any(strcmp(charargin,'EdgeTrack'))
                trackIdx = obj.find('NonEdgeTrack==0');
           end
           
           % Plot multiple tracks if they're present
           if gt(numel( trackIdx ),1)
               set(gca,'NextPlot','add');
               thisplot = arrayfun(@(thisIdx) plot( obj.tracksTable(thisIdx,:).x{1},...
                     obj.tracksTable(thisIdx,:).y{1} ), [1:numel(unique(trackIdx))] );
           else
               thisplot = plot( obj.tracksTable(trackIdx,:).x{1},...
                     obj.tracksTable(trackIdx,:).y{1} );
           end
           
           axis image
           
           mytitle = regexprep(regexprep(obj.metadata.Directory,'\',' '),'_',' ');
           
           set(gcf,'position',[120,350,1370,400]);
           title( sprintf('%s (%i)', mytitle, obj.metadata.obj_idx) );
           
           if any(cell2mat(strfind(charargin,'save')))
                saveas(gcf, sprintf('Trackboundaries_Figure_%i.svg',obj.metadata.obj_idx) );
                delete(fig)
           end
           
        end
        
    end

        methods(Static)
                  
            function tracksTable = extractPosition(obj, tracksFinal)
                % Function which dissects tracksFinal into x,y (called from within the constructor function)
                % Extract each x and y position from tracksFinal.tracksCoordAmpCG
                % [x0,y0,...(6 entries),x1,y1,...]
                Ntracks = max(size(tracksFinal));
                if isfield(tracksFinal,'ids'); tracksFinal = rmfield(tracksFinal,'ids'); end
                start_time = arrayfun(@(x) tracksFinal(x).seqOfEvents(1), [1:Ntracks] );
                xy = struct2cell( rmfield(tracksFinal,{'tracksFeatIndxCG','seqOfEvents'}) );
                % Auxiliary functions to pull x,y,I from the UTrack data structure
                % U-Track info: 
                % "For every frame, the matrix stores the particle’s x-coordinate, y-coordinate, z-coordinate
                % (0 if 2D),amplitude, x-coordinate standard deviation, ycoordinate standard deviation,
                % z-coordinate standard deviation (0 if 2D) and amplitude standard deviation."
                
                getX = @(x) [ x(1:8:size(x,2)) ]; 
                getY = @(x) [ x(2:8:size(x,2)) ];
                getI = @(x) [ x(4:8:size(x,2)) ];
                getI_sd = @(x) [ x(8:8:size(x,2)) ];
                
                [tmp_x,tmp_y,tmp_I,tmp_Isd] = deal( cellfun(getX, xy, 'UniformOutput', false),...
                                    cellfun(getY, xy, 'UniformOutput', false),...
                                    cellfun(getI, xy, 'UniformOutput', false),...
                                    cellfun(getI_sd, xy, 'UniformOutput', false) );
                [x,y,I,Isd] = deal( squeeze(tmp_x), squeeze(tmp_y), squeeze(tmp_I), squeeze(tmp_Isd) );
                if size(tmp_x,1)==1; x=tmp_x'; y = tmp_y'; I = tmp_I'; Isd = tmp_Isd'; end
                
                tracksTable_ = table( x,...
                                      y,...
                                      start_time',...
                                      I,...
                                      Isd );
                
                % Add column for length of each track
                tracksTable_.trackLength = arrayfun( @(idx) numel( tracksTable_(idx,:).x{1} ), ...
                        [1:size(tracksTable_,1)] )';
                tracksTable_.I_mean = cellfun(@(x) mean(x), tracksTable_.I );
            
            end
            % END OF extractPosition function %
            
            
            
            function fileStruct = findFiles(Directory, FilesToFind)
% This function locates desired files (called from within the constructor function)    
                for i = 1:numel(FilesToFind)
                    filetmp_ = dir( fullfile(Directory,'**',FilesToFind{i}) );
                    if ~isempty(filetmp_)
                        fprintf('Located %s\n', regexprep( fullfile(  filetmp_.folder, filetmp_.name ), regexprep(Directory,'\\','\\\'), '') ); 
                        
                        if size(filetmp_,1)>1
                            %fprintf('Located %i files\n', size(filetmp_,1) ); 
                            %for j = 1:size(filetmp_,1)
                            %   thisFile = fullfile( filetmp_(j).folder,filetmp_(j).name );
                            %   fprintf('%i) %s\n',j,thisFile); 
                            %end
                            %select_ = input('Which file to select? (Enter a number):');
                            select_ = size(filetmp_,1);
                        else
                            select_ = 1;
                        end
                        
                        fileStruct.( regexprep( FilesToFind{i}, '.mat', '') ) = 1;
                        fileStruct.( regexprep( strcat('address_',FilesToFind{i}), '.mat', '') ) = fullfile(  filetmp_(select_).folder, filetmp_(select_).name ); 
                    else
                        fileStruct.( regexprep( FilesToFind{i}, '.mat', '') ) = 0;
                    end
                end
            end
            % END OF findFiles function %
            
            
        end
        
    end