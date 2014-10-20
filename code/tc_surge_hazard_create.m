function hazard=tc_surge_hazard_create(hazard,hazard_set_file,suppress_plots)
% climada storm surge TS hazard event set
% NAME:
%   tc_surge_hazard_create
% PURPOSE:
%   create a storm surge (TS) hazard event, based on an existing
%   tropical cyclone (TC) hazard event set
%
%   two steps:
%   1) obtain (high res) bathymetry, calls etopo_get
%   2) convert all TC footprints into TS footprints
%   see CORE_CONVERSION in code below for the conversion formula
%
%   Note on bathymetry data:
%   ETOPO1.nc was originally named ETOPO1_Ice_g_gmt4.grd and is the full
%   high-res topography dataset (933.5 MB!) obtained directly from
%   http://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz
%
%   see tc_surge_TEST for a testbed for this code
%
% CALLING SEQUENCE:
%   hazard=tc_surge_hazard_create(hazard,hazard_set_file)
% EXAMPLE:
%   hazard=tc_surge_hazard_create
% INPUTS:
%   hazard: an already existing tropical cyclone (TC) hazard event set (a
%       TC hazard structure)
%       > prompted for if not given (for .mat file containing a TC hazard event set)
%   hazard_set_file: the name of the newly created storm surge (TS) hazard
%       event set
%       > promted for if not given
%   suppress_plots: =1, do not show any figures, =o: show plots (default)
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   hazard: a hazard event set, see core climada doc
%       also written to a .mat file (see hazard_set_file)
%       NOTE: for memory allocation reasons, the input hazard is used and
%       modified to create the output hazard
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20140421
% David N. Bresch, david.bresch@gmail.com, 20141017, module path relative
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('hazard','var'),hazard=[];end
if ~exist('hazard_set_file','var'),hazard_set_file=[];end
if ~exist('suppress_plots','var'),suppress_plots=0;end

% PARAMETERS
%
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];


% prompt for TC hazard event set if not given
if isempty(hazard) % local GUI
    TC_hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(TC_hazard_set_file, 'Select a TC hazard event set:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        TC_hazard_set_file=fullfile(pathname,filename);
    end
    load(TC_hazard_set_file);
end

% prompt for hazard_set_file if not given
if isempty(hazard_set_file) % local GUI
    hazard_set_file=[climada_global.data_dir filesep 'hazards' filesep 'TS_hazard.mat'];
    [filename, pathname] = uiputfile(hazard_set_file, 'Save new TS hazard event set as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard_set_file=fullfile(pathname,filename);
    end
end

% prep the region we need (rectangular region encompassing the hazard centroids)
centroids_rect=[min(hazard.lon) max(hazard.lon) min(hazard.lat) max(hazard.lat)];

% 1) create the bathymetry file
% -----------------------------

bb=1; % 1 degree of bathy outside centroids (to allow for smooth interp)
bathy_coords=[centroids_rect(1)-bb centroids_rect(2)+bb centroids_rect(3)-bb centroids_rect(4)+bb];

% cut the bathymtery data out of the global topography dataset

% file the bathymetry gets stored in ( might not be used later on, but
% saved to speed up re-creation of the TS hazard event set
[~,fN]=fileparts(hazard_set_file);
Bathymetry_file=[module_data_dir filesep fN '_bathy.mat'];

if ~exist(Bathymetry_file,'file')
    
    if ~exist('etopo_get','file')
        % safety to inform the user in case he misses the ETOPO module
        cprintf([1,0.5,0],'ERROR: no etopo_get function found. Please download from github and install the climada etopo module\nhttps://github.com/davidnbresch/climada_module_etopo\n');
        hazard=[];
        return
    end
    BATI=etopo_get(bathy_coords);
    
    fprintf('saving bathymetry as %s (you might later delete this file)\n',Bathymetry_file);
    save(Bathymetry_file,'BATI');
    
    if ~suppress_plots
        figure('Name','Bathymetry','Color',[1 1 1]);
        pcolor(BATI.x,BATI.y,BATI.h)
        hold on
        shading interp
        colorbar
        caxis([-200 20]) % a reasonable colorbar for elevation
        axis equal
        climada_plot_world_borders
        axis(bathy_coords);
    end
    
    % netcdf.close(bathy_data_file) % seems not be required
else
    fprintf('reading bathymetry from %s\n',Bathymetry_file);
    load(Bathymetry_file);
end

% 2) create the storm surge (TS) hazard event set
% -----------------------------------------------

% start from the 'mother' TC hazard event set

if isfield(hazard,'onLand'),hazard=rmfield(hazard,'onLand');end % re-create to be on the safe side
hinterp=interp2(BATI.x,BATI.y,BATI.h,hazard.lon,hazard.lat);
hazard.onLand=hazard.lon.*0+1; % allocate
hazard.onLand(hinterp<0)=0; % water points

hazard.peril_ID='TS'; % replace TC with TS
if isfield(hazard,'windfield_comment'),hazard=rmfield(hazard,'windfield_comment');end % remove
hazard.surgefield_comment=sprintf('created based on TC using proxy surge height and bathymetry %s',BATI.sourcefile);
hazard.comment=sprintf('TS hazard event set, generated %s',datestr(now));

% map windspeed onto surge height (the CORE_CONVERSION, see at botton of file, too)
% ===============================
arr_nonzero=find(hazard.intensity); % to avoid de-sparsify all elements
hazard.intensity(arr_nonzero)=0.1023*(max(hazard.intensity(arr_nonzero)-26.8224,0))+1.8288; % m/s converted to m surge height

% grid bathymetry to centroids (convert n x m BATI to vectors in call)
xx=double(reshape(BATI.x,numel(BATI.x),1));
yy=double(reshape(BATI.y,numel(BATI.y),1));
zz=double(reshape(BATI.h,numel(BATI.h),1));
hazard.bathy=griddata(xx,yy,zz,hazard.lon,hazard.lat);

hazard.elev=max(hazard.bathy,0); % only points above sea level

% subtract elevation above sea level from surge height
t0       = clock;
mod_step = 10; % first time estimate after 10 tracks, then every 100
n_events=size(hazard.intensity,1);
n_centroids=size(hazard.intensity,2);

% as the innermost loop is vectorized, it shall be the one repeated most:

if n_events<n_centroids % loop over events, since less events than centroids
    
    msgstr   = sprintf('processing %i events',n_events);
    fprintf('%s (updating waitbar with estimation of time remaining every 100th event)\n',msgstr);
    h        = waitbar(0,msgstr);
    for event_i=1:n_events
        arr_i=find(hazard.intensity(event_i,:)); % to avoid de-sparsify all elements
        hazard.intensity(event_i,arr_i)=max(hazard.intensity(event_i,arr_i)-hazard.elev(arr_i),0);
        
        if mod(event_i,mod_step)==0
            mod_step = 100;
            t_elapsed = etime(clock,t0)/event_i;
            n_remaining = n_events-event_i;
            t_projected_sec = t_elapsed*n_remaining;
            if t_projected_sec<60
                msgstr = sprintf('est. %3.0f sec left (%i/%i events)',t_projected_sec, event_i, n_events);
            else
                msgstr = sprintf('est. %3.1f min left (%i/%i events)',t_projected_sec/60, event_i, n_events);
            end
            waitbar(event_i/n_events,h,msgstr); % update waitbar
        end
        
    end % event_i
    
else % loop over centroids, since less centroids than events
    
    msgstr   = sprintf('processing %i centroids',n_centroids);
    fprintf('%s (updating waitbar with estimation of time remaining every 100th centroid)\n',msgstr);
    h        = waitbar(0,msgstr);
    for centroid_i=1:n_centroids
        arr_i=find(hazard.intensity(:,centroid_i)); % to avoid de-sparsify all elements
        hazard.intensity(arr_i,centroid_i)=max(hazard.intensity(arr_i,centroid_i)-hazard.elev(centroid_i),0);
        
        if mod(centroid_i,mod_step)==0
            mod_step = 100;
            t_elapsed = etime(clock,t0)/centroid_i;
            n_remaining = n_centroids-centroid_i;
            t_projected_sec = t_elapsed*n_remaining;
            if t_projected_sec<60
                msgstr = sprintf('est. %3.0f sec left (%i/%i centroids)',t_projected_sec, centroid_i, n_centroids);
            else
                msgstr = sprintf('est. %3.1f min left (%i/%i centroids)',t_projected_sec/60, centroid_i, n_centroids);
            end
            waitbar(centroid_i/n_centroids,h,msgstr); % update waitbar
        end
        
    end % event_i
    
end
close(h); % dispose waitbar

if isfield(hazard,'filename'),hazard.filename_source=hazard.filename;end
hazard.filename=hazard_set_file;
hazard.date=datestr(now);
hazard.matrix_density=nnz(hazard.intensity)/numel(hazard.intensity);
hazard.units='m'; % store the SI unit of the hazard intensity
if ~isfield(hazard,'orig_event_count') % fix a minor issue with some hazard sets
    if isfield(hazard,'orig_event_flag')
        fprintf('field hazard.orig_event_count inferred from hazard.orig_event_flag\n')
        hazard.orig_event_count=sum(hazard.orig_event_flag);
    else
        fprintf('WARNING: no field hazard.orig_event_flag\n')
    end
end

fprintf('saving TS surge hazard set as %s\n',hazard_set_file);
save(hazard_set_file,'hazard');

%%fprintf('TS: max(max(hazard.intensity))=%f\n',full(max(max(hazard.intensity)))); % a kind of easy check

if ~suppress_plots,climada_hazard_plot(hazard,0);end % show max surge over ALL events

return

% =====================================
% wind speed to surge height conversion (CORE_CONVERSION)
% =====================================

% copy/paste below into MATLAB command window to run

mph2ms=0.44704;
f2m=0.3048;

% the points read from the SLOSH graph
v0=60*mph2ms;
v1=140*mph2ms;
s0=6*f2m;
s1=18*f2m;

% the parameters for the linear function
a=(s1-s0)/(v1-v0)
s0
v0

figure('Name','windspeed to surge height conversion','Color',[1 1 1])
hold on
v=20:100;       plot(v,          a*(v-v0)          +s0,'-r','LineWidth' ,3);
vmph=60:20:140; plot(vmph*mph2ms,a*(vmph*mph2ms-v0)+s0,'.b','MarkerSize',10);
legend('conversion','SLOSH points')
xlabel('wind speed [m/s]')
ylabel('surge height [m]')
