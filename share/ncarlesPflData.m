classdef ncarlesPflData
% ncarlesPflData is a class for postprocessing the profile
%   data from NCAR LES
%
%   Qing Li, 160501
%            160510, support for nscl>1
%            160803, support for varying utau and stokes

%% properties
properties
    caseName        % case name
    dataDir         % file directory
    fPrefix         % file prefix
    tStart          % start time step
    tEnd            % end time step
    nnz             % vertical grid size
    nscl            % number of scalers in t
    ntime           % time dimension size
    time            % time dimension
    filename        % data filename
end
properties (Constant)
    % gravity (m/s^2)
    g = 9.81;
    % water density (kg/m^3)
    rho = 1000;
    % specific heat of water (J/Kg/K)
    cp = 4200;
    % temperature expansion coefficent (K^{-1})
    alpha = 1/5000;
    % von Karman's constant
    vk = 0.4;
    % parameter list
    pars = cellstr(['utau    ';'wtsfc   ';'fcor    ';'stokess ';...
                    'stokesw ']);
    % profile variable list
    vars = cellstr(['wtle    ';'wtsb    ';'englez  ';'engz    ';...
                    'engsbz  ';'t_dsle  ';'t_rprod ';'t_sprod ';...
                    't_stokes';'t_tran  ';'t_wq    ';'t_wp    ';...
                    'uwle    ';'uwsb    ';'vwle    ';'vwsb    ';...
                    'ups     ';'vps     ';'wps     ';'tps     ';...
                    'uxym    ';'vxym    ';'wxym    ';'txym    ';...
                    'wcube   ';'tcube   ';'wfour   ';'shrz    ';...
                    'uvle    ';'dudz    ';'dvdz    ';'t_diss  ';...
                    'uuwle   ';'uvwle   ';'uwwle   ';'vvwle   ';...
                    'vwwle   ';'udpdx   ';'vdpdx   ';'wdpdx   ';...
                    'udpdy   ';'vdpdy   ';'wdpdy   ';'udpdz   ';...
                    'vdpdz   ';'wdpdz   ';'ttau11  ';'ttau12  ';...
                    'ttau13  ';'ttau22  ';'ttau23  ';'ttau33  ';...
                    'dsle11  ';'dsle12  ';'dsle13  ';'dsle22  ';...
                    'dsle23  ';'dsle33  ';'utle    ';'vtle    ']);
    % number of variables that involve the scaler variable t
    nvscl = 5;
    % z flags: 1, z_w; 0, z_u
    f_z = [1;1;1;1;...
           0;1;1;1;...
           1;1;1;1;...
           1;1;1;1;...
           0;0;1;0;...
           0;0;1;0;...
           1;0;1;1;...
           0;1;1;1;...
           1;1;1;1;...
           1;0;0;0;...
           0;0;0;1;...
           1;1;0;0;...
           0;0;0;0;...
           0;0;0;0;...
           0;0;0;0];
end

%% methods
methods
    %% Constructor
    function obj = ncarlesPflData(varargin)
    % ncarlesPflData creates an object for the profile data from
    %   the NCAR LES, and initializes the property values
    %
    %   obj = ncarlesPflData('property',value,[...])

        % get field names
        fNames = fieldnames(obj);
        % count arguments
        nArgs = length(varargin);
        if round(nArgs/2)~=nArgs/2
            error('Needs propertyName/propertyValue pairs\n');
        end
        % fill in properties
        for pair = reshape(varargin,2,[])
            inpName = pair{1};
            if any(strcmp(inpName,fNames))
                obj.(inpName) = pair{2};
            else
                error(['%s is not a recognized '...
                      'property name, run fieldnames(obj) '...
                      'to get a list of property names\n'],inpName);
            end
        end
        if ~isempty(obj.caseName)
            inFile = obj.getFileName;
            obj.filename = inFile;
            ndim = obj.dSize;
            obj.nnz = ndim(1);
            obj.ntime = ndim(2);
            obj.nscl = ndim(3);
            tmp = obj.getVar('time');
            obj.time = tmp';
        end
    end % function ncarlesPflData

    %% getFileName
    function name = getFileName(obj)
    % getFileName returns the name of the file, determined from
    %   from the data directory, case name, file prefix and
    %   time stamps.

        % time stamp
        tStartStr = sprintf('%06d',obj.tStart);
        tEndStr = sprintf('%06d',obj.tEnd);
        %filename
        name = [obj.dataDir '/' obj.caseName '/',...
             obj.fPrefix '.' tStartStr '.' tEndStr '.nc'];
    end % function getFileName

    %% dSize
    function ndim = dSize(obj)
    % dSize returns the vertical, scaler index and time dimension sizes

        % get file name
        inFile = obj.filename;
        % file information
        finfo = ncinfo(inFile);
        % dimensions
        n = length(finfo.Dimensions);
        ndim = zeros([1,n-1]);
        for i=1:n
            switch finfo.Dimensions(i).Name
                case 'z_w'
                    ndim(1) = finfo.Dimensions(i).Length;
                case 'time'
                    ndim(2) = finfo.Dimensions(i).Length;
                case 'scl'
                    ndim(3) = finfo.Dimensions(i).Length;
            end
        end
        ndim(ndim==0) = NaN;
    end

    %% getVar
    function var = getVar(obj,fvar)
    % getVar read the profile data and return the timeseries of the
    %   profiles for variable fvar

        if ~ischar(fvar)
            error('Input argument should be a string.');
        end
%         if ~any(strcmp(fvar,obj.vars))
%             error(['Variable %s not found. Check f.vars for available ',...
%                 'variables.'],fvar);
%         end
        var = double(ncread(obj.filename,fvar));
    end % function getVar"

    %% getTimeIndex
    function ts = getTimeIndex(obj,varargin)
    % getTimeIndex returns the time index closest to the given time
    %
    %   getTimeIndex(time) in which time is measured in second since the
    %     start of simulation, or
    %   getTimeIndex(measure,time) in which measure could be
    %     'hour': hours since the start
    %     'second': seconds since the start

        % count arguments
        nArgs = length(varargin);
        if nArgs==1 && isnumeric(varargin{1})
            measure = 'second';
            t_in = varargin{1};
        elseif nArgs==2 && ischar(varargin{1}) && isnumeric(varargin{2})
            measure = varargin{1};
            t_in = varargin{2};
        else
            error('Invalid arguments\n');
        end
        % find the index
        if strcmpi(measure,'second')
            t_target = t_in;
        elseif strcmpi(measure,'hour')
            t_target = t_in.*3600;
        else
            error('Supported unit of time: second, hour. Got %s',measure);
        end
        [~,ts] = min(abs(obj.time-t_target));
    end % function getTimeIndex

    %% getParameters
    function dat = getParameters(obj)
    % getParameters returns the parameters in a structure

        % number of parameters
        nPar = numel(obj.pars);
        % loop over parameters
        for ip=1:nPar
            try
                dat.(obj.pars{ip}) = obj.getVar(obj.pars{ip});
            catch
                warning(['Parameter ' obj.pars{ip} ' is not found '...
                         'in file. Skipping...']);
            end
        end % par loop
        % Monin-Obukhov length scale
        dat.amonin = -dat.utau(end).^3 ...
                   ./(obj.vk.*obj.g.*obj.alpha.*dat.wtsfc);
        % Ekman-layer e-folding scale
        dat.lEkman = 0.25.*dat.utau(end)./dat.fcor;
%         % turbulent Langmuir number
%         dat.La_turb = sqrt(dat.utau(end)./dat.stokess);
%         % Stokes depth
%         if isinf(dat.La_turb)
%             dat.delta = NaN;
%         else
%             dat.delta = 0.5./dat.stokesw;
%         end
%         % Hoenikker number
%         dat.Ho = 2.*obj.alpha.*obj.g.*dat.wtsfc./dat.stokesw ...
%                ./dat.stokess./dat.utau(end).^2;
        % total time
        dat.tend = obj.time(obj.ntime);
    end % function getParameters

    %% getProfiles
    function dat = getProfiles(obj,varargin)
    % getProfiles returns a structure holding the averaged
    %   profiles of variables. The profile is first stretched at each
    %   time step based on the current mixed layer depth, then averaged.

        % count number of arguments
        nArgs = length(varargin);
        if nArgs ~= 2
            error('Require 2 arguments, %d is found.', nArgs);
        end
        % get dt
        dt = obj.time;
        dt(2:end-1) = 0.5.*(obj.time(3:end)-obj.time(1:end-2));
        dt(1) = dt(2);
        dt(end) = dt(end-1);
        % averaging window
        [its, ite] = obj.parsingAvgWindow(varargin{1:2});
        % number of profile variables
        nVar = numel(obj.vars);
        % read vertical grid
        dat.z_u = obj.getVar('z_u')';
        dat.z_w = obj.getVar('z_w')';
        tmp = obj.getVar('stokes');
        n = length(size(tmp));
        if n==1
            dat.stokes = tmp;
        elseif n==2
            dat.stokes = tmp(:,end)';
        end
        dat.stokes_dim = 'z_u';
        h = -obj.getBLDepth('max_Nsquared');
        mean_h = obj.wgtMean(h(its:ite),dt(its:ite),2);
        dat.hb = h;
        dat.mean_hb = mean_h;
        zz_in = zeros(1,obj.nnz+1);
        zn_u = dat.z_u./mean_h;
        zn_w = dat.z_w./mean_h;
        % loop over variables
        for iv=1:nVar
        try
            if obj.f_z(iv)==1
                dat.([obj.vars{iv} '_dim']) = 'z_w';
                zz_in(2:end) = dat.z_w;
                zz_out = zn_w;
            else
                dat.([obj.vars{iv} '_dim']) = 'z_u';
                zz_in(2:end) = dat.z_u;
                zz_out = zn_u;
            end
            tmp0 = obj.getVar(obj.vars{iv});
            ndim = size(tmp0);
            n = length(ndim);
            tmp2 = zeros(ndim);
            % do average
            if n==2
                % prevent NaNs for shrinking
                tmp = zeros(ndim(1)+1,ndim(2));
                tmp(2:end,:) = tmp0;
                tmp(1,:) = tmp0(1,:);
                % stretching in z
                for it = its:ite
                    tmp2(:,it) = interp1(zz_in./h(it),...
                        tmp(:,it),zz_out,'linear',NaN);
                end
                dat.(obj.vars{iv}) = ...
                    obj.wgtMean(tmp2(:,its:ite),dt(its:ite),2);
            elseif n==3
                % prevent NaNs for shrinking
                tmp = zeros(ndim(1)+1,ndim(2),ndim(3));
                tmp(2:end,:,:) = tmp0;
                tmp(1,:,:) = tmp0(1,:,:);
                % stretching in z
                for it = its:ite
                    tmp2(:,:,it) = interp1(zz_in./h(it),...
                        tmp(:,:,it),zz_out,'linear',NaN);
                end
                dat.(obj.vars{iv}) = ...
                    obj.wgtMean(tmp2(:,:,its:ite),dt(its:ite),3);
            else
                error(['Invalid dimension size %i of variable '...
                    '%s\n'],n,obj.vars{iv});
            end
        catch
            warning(['Variable ' obj.vars{iv} ' is not found '...
                     'in file. Set to zero...']);
            dat.(obj.vars{iv}) = 0;
            dat.([obj.vars{iv} '_dim']) = 'none';
        end
        end % var loop
    end % function getProfiles

    %% getStatEntrainment
    function dat = getStatEntrainment(obj,varargin)
    % getStatEntrainment returns the 25th, 50th, 75th percentile of the
    %       buoyancy flux at the base of the OSBL (boundary layer depth
    %       defined as the depth at which the buoyancy flux reaches the
    %       minimum), as well as the mean

        % count number of arguments
        nArgs = length(varargin);
        if nArgs ~= 2
            error('Require 2 arguments, %d is found.', nArgs);
        end
        % get dt
        dt = obj.time;
        dt(2:end-1) = 0.5.*(obj.time(3:end)-obj.time(1:end-2));
        dt(1) = dt(2);
        dt(end) = dt(end-1);
        % averaging window
        [its, ite] = obj.parsingAvgWindow(varargin{1:2});
        % z
        z_w = obj.getVar('z_w');
        % mixed layer depth
        h = obj.getBLDepth('max_Nsquared');
        mean_h = obj.wgtMean(h(its:ite),dt(its:ite),2);
        zn_w = z_w./mean_h;
        % buoyancy flux
        tmp = obj.getVar('wtle');
        wtle = squeeze(tmp(:,1,:));
        tmp = obj.getVar('wtsb');
        wtsb = squeeze(tmp(:,1,:));
        wbpfl = (wtle+wtsb).*obj.alpha.*obj.g;
        zz_in = z_w;
        zz_out = zn_w;
        tmp2 = wbpfl;
        for it = its:ite
            tmp2(:,it) = interp1(zz_in./h(it),...
                wbpfl(:,it),zz_out,'linear',NaN);
        end
        meandat = obj.wgtMean(tmp2(:,its:ite),dt(its:ite),2);
        [min_wb,inds] = min(meandat);
        wb = tmp2(inds,its:ite);
        dat.mean = min_wb;
        dat.p25 = prctile(wb,25);
        dat.median = median(wb);
        dat.p75 = prctile(wb,75);
    end % function getStatEntrainment

    %% wgtMean
    function datOut = wgtMean(~,datIn,wgt,idim)
    % wgtMean returns the time mean for uneven spaced data. Each sample
    %   is weighted by the time step

        % dimension size
        nsize = size(datIn);
        nd = numel(nsize);
        if idim > nd
            error('DIM exceeds the maximum dimension');
        end
        wsize = size(wgt);
        nw = numel(wgt);
        if nsize(idim) ~= nw
            error('The dimension of wgt and datIn do not match');
        elseif max(wsize) ~= nw
            error('wgt should be one dimensional');
        elseif wsize(2) == nw
            wgt = shiftdim(wgt,1);
        end
        nsft = idim-1;
        totWgt = sum(wgt(:));
        dat = shiftdim(datIn,nsft);
        tmp = sum(dat.*wgt)./totWgt;
        datOut = squeeze(tmp);
    end % function timeMean

    %% getTimeSeries
    function dat = getTimeSeries(obj,varargin)
    % getTimeSeries returns the time series of profile variables
    %   1) at a certain depth if one depth (m) is passed in
    %   2) averaged over a layer if the upper and lower bounds of
    %      the layer are passed in
    %   3) averaged over the whole boundary layer if no argument is
    %      passed in or a boundary layer definition is passed in

        % count arguments
        nArgs = length(varargin);
        if nArgs > 2
            error(['Maximum number of argument is 2, ' ...
                 '%d is found.'], nArgs);
        end
        % number of profile variables
        nVar = numel(obj.vars);
        % get input file name
        inFile = obj.getFileName;
        % read in depth
        zu = ncread(inFile, 'z_u');
        zw = ncread(inFile, 'z_w');
        % read in scl
        if ~isnan(obj.nscl)
            dat.scl = ncread(inFile, 'scl');
        end
        % read in time
        dat.time = ncread(inFile, 'time');
        if nArgs == 1 && isnumeric(varargin{1})
            d = varargin{1};
            if d < min(zw) || d > max(zu)
                error('Desired depth outside the domain');
            end
            % loop over variables
            for iv=1:nVar
                try
                    tmp   = ncread(inFile,obj.vars{iv});
                    % get timeseries, do interpolation
                    if obj.f_z(iv)==1
                        dat.(obj.vars{iv}) = ...
                            squeeze(interp1(zw,tmp,d,'linear'));
                    else
                        dat.(obj.vars{iv}) = ...
                            squeeze(interp1(zu,tmp,d,'linear'));
                    end
                catch
                    warning(['Variable ' obj.vars{iv} ...
                         ' is not found '...
                         'in file. Set to zero...']);
                    dat.(obj.vars{iv}) = zeros(1,obj.ntime);
                end
            end % var loop
            dat.depth = d;
        elseif nArgs == 2 && isnumeric(varargin{1}) && ...
                isnumeric(varargin{2})
            d1 = min(varargin{:});
            d2 = max(varargin{:});
            if d1 < min(zw) || d2 > 0
                error('Desired depth outside the domain');
            end
            for iv=1:nVar
                try
                    tmp   = ncread(inFile,obj.vars{iv});
                    % get timeseries, do interpolation
                    if obj.f_z(iv)==1
                        dat.(obj.vars{iv}) = obj.zAvg(zw,tmp,d1,d2,'w');
                    else
                        dat.(obj.vars{iv}) = obj.zAvg(zu,tmp,d1,d2,'u');
                    end
                catch
                    warning(['Variable ' obj.vars{iv} ...
                         ' is not found '...
                         'in file. Set to zero...']);
                    dat.(obj.vars{iv}) = zeros(1,obj.ntime);
                end
            end % var loop
            dat.depth1 = d1;
            dat.depth2 = d2;
        else
            % arguments are for getBLDepth, pass them on to getBLDepth
            h = obj.getBLDepth(varargin{:});
            % loop over variables
            dat.h = h;
            for iv=1:nVar
                try
                    tmp = ncread(inFile,obj.vars{iv});
                    % get timeseries, do interpolation
                    if obj.f_z(iv)==1
                        dat.(obj.vars{iv}) = obj.zAvg(zw,tmp,h,0,'w');
                    else
                        dat.(obj.vars{iv}) = obj.zAvg(zu,tmp,h,0,'u');
                    end
                catch
                    warning(['Variable ' obj.vars{iv} ...
                         ' is not found '...
                         'in file. Set to zero...']);
                    dat.(obj.vars{iv}) = zeros(1,obj.ntime);
                end
            end % var loop
        end % nArgs
    end % function getTimeSeries

    %% getBLDepth
    function h = getBLDepth(obj,varargin)
    % getBLDepth returns a time series of boundary layer depth
    %   calculated by method:
    %   't_threshold'  - the depth at which the mean temperature (txym)
    %                    differs by a critical value (default 0.1 K,
    %                    pass in a value to override it) from its
    %                    surface value
    %   'max_Nsquared' - the depth at which the gradient of the buoyancy
    %                    reaches its maximum
    %   'kpp'          - diagnose the KPP boundary layer (Large et al., 1994)
    %   'read'         - read from file

        % count arguments
        nArgs = length(varargin);
        if nArgs==0
            method = 'max_Nsquared';
            fprintf('getBLDepth: max_Nsquared, default\n');
        elseif nArgs==1
            if strcmpi(varargin{1},'t_threshold')
                method = 't_threshold';
                tcrit = 0.1;
                fprintf(['getBLDepth: t_threshold, ' ...
                    'default crit. t = %g K\n'],tcrit);
            elseif strcmpi(varargin{1},'max_Nsquared')
                method = 'max_Nsquared';
                fprintf('getBLDepth: max_Nsquared\n');
            elseif strcmpi(varargin{1},'max_Nsquared_smth')
                method = 'max_Nsquared_smth';
                fprintf('getBLDepth: max_Nsquared_smth\n');
            elseif strcmpi(varargin{1},'kpp')
                method = 'kpp';
                fprintf('getBLDepth: diagnose kpp BLDepth\n');
            elseif strcmpi(varargin{1},'read')
                method = 'read';
                fprintf('getBLDepth: read bld from file\n');
            else
                error('getBLDepth: Method %s not supported.\n',varargin{1});
            end
        elseif nArgs==2 % override the critical temperature
            if strcmpi(varargin{1},'t_threshold') && ...
                    isnumeric(varargin{2})
                method = 't_threshold';
                tcrit = varargin{2};
                fprintf('getBLDepth: t_threshold, crit. t = %g K\n',tcrit);
            else
                error(['Either method %s is not supported or ',...
                      'the critical temperature is not a number'],...
                      varargin{1});
            end
        else
            error(['No more than 2 arguments are allowed, '...
                '%d recieved'],nArgs)
        end
        % input file name
        nt = obj.ntime;
        nz = obj.nnz;
        h = zeros([1,nt]);
        switch method
        case 't_threshold'
            txym = obj.getVar('txym');
            ndim = size(txym);
            if length(ndim)==3
                txym = squeeze(txym(:,1,:));
            end
            zu = obj.getVar('z_u');
            for i = nt-1:nt
                h(i) = obj.getThresholdDepth(txym(:,i),zu,txym(1,i)-tcrit);
            end
        case 'max_Nsquared'
            txym = obj.getVar('txym');
            ndim = size(txym);
            if length(ndim)==3
                txym = squeeze(txym(:,1,:));
            end
            zu = obj.getVar('z_u');
            zw = obj.getVar('z_w');
            for i = 1:nt
                % To be robust, should be away from the upper and
                % lower boundary
                ibc0 = 3;
                ibc1 = 5;
                tmp = (txym(ibc0+1:end-ibc1,i)...
                    -txym(ibc0:end-ibc1-1,i)) ...
                    ./(zu(ibc0+1:end-ibc1)-zu(ibc0:end-ibc1-1));
                [~,ind] = max(abs(tmp));
                h(i) = zw(ibc0-1+ind);
                % debug
                % if h(i)<-50
                %     fprintf('h = %g\n',h(i));
                %     for k=1:numel(tmp)
                %         fprintf('z = %g, N2 = %g\n',zw(k),tmp(k));
                %     end
                % end
                % debug
            end
        case 'max_Nsquared_smth'
            txym = obj.getVar('txym');
            ndim = size(txym);
            if length(ndim)==3
                txym = squeeze(txym(:,1,:));
            end
            zu = obj.getVar('z_u');
            zw = obj.getVar('z_w');
            zwnew = obj.refineZ(zw,10);
            for i = 1:nt
                N2 = obj.NSquared(txym(:,i),zu);
                % To be robust, should be away from the upper and
                % lower boundary
                ibc0 = 3;
                ibc1 = 5;
                N2(1:ibc0) = N2(ibc0+1);
                N2(end-ibc1:end) = N2(end-ibc1-1);
                N2new = interp1(zw,N2,zwnew,'spline');
                [~,ind] = max(abs(N2new));
                h(i) = zwnew(ind);
                % debug
                % if h(i)<-50
                %     fprintf('h = %g\n',h(i));
                %     for k=1:numel(tmp)
                %         fprintf('z = %g, N2 = %g\n',zw(k),tmp(k));
                %     end
                % end
                % debug
            end
        case 'kpp'
            surf_layer_ext = obj.kpp_get_constant('surf_layer_ext');
            datpar = obj.getParameters;
            ustar = datpar.utau(end);
            bf = -obj.g.*obj.alpha.*datpar.wtsfc;
            uxym = obj.getVar('uxym');
            vxym = obj.getVar('vxym');
            bxym = obj.getVar('txym').*obj.g.*obj.alpha;
            ndim = size(bxym);
            if length(ndim)==3
                bxym = squeeze(bxym(:,1,:));
            end
            zu = obj.getVar('z_u');
            zw = obj.getVar('z_w');
            w_s = zeros(size(zu));
            for k=1:nz
                w_s(k) = obj.kpp_turbulent_scale(surf_layer_ext,...
                                -zu(k),bf,ustar,'s');
            end
            for i=1:nt
                h(i) = obj.kpp_boundary_layer_depth(bxym(:,i),...
                                                    uxym(:,i),...
                                                    vxym(:,i),...
                                                    zu,zw,w_s);
            end
        case 'read'
            tmp = obj.getVar('bld');
            h = tmp';
        end
    end % function getBLDepth

    %% getWe
    function w = getWe(obj,varargin)
    % getWe returns the timeseries of the entrainment velocity We=dh/dt

        h = obj.getVar('bld');
        t = obj.time;
        w = zeros(1,obj.ntime);
        w(2:end-1) = (h(3:end)-h(1:end-2))./(t(3:end)-t(1:end-2));
        w(1) = w(2);
        w(end) = w(end-1);
    end % function getWe

    %% parsingAvgWindow
    function [its, ite] = parsingAvgWindow(obj,varargin)
    % parsingAvgWindow returns the averaging window in time indices

        if ~(ischar(varargin{1}) && isnumeric(varargin{2}))
            error(['Please pass in the average window, in the format of '...
                'either {''hour''/''second'',[i_start, i_end]} '...
                'or {''inertial'',[i,j]} (averaging over the last jth to '...
                'ith inertial period(s)).']);
        end
        if strcmpi(varargin{1},'hour') || strcmpi(varargin{1},'second')
            measure = varargin{1};
            tavg_start = varargin{2}(1);
            tavg_end = varargin{2}(2);
            fprintf(['Averaging over perids (%s): %g to %g\n'],...
               measure,tavg_start,tavg_end);
            % average window, start and end indices
            its = obj.getTimeIndex(measure,tavg_start);
            ite = obj.getTimeIndex(measure,tavg_end);
        elseif strcmpi(varargin{1},'inertial')
            par = obj.getParameters;
            nip_start = varargin{2}(2);
            nip_end = varargin{2}(1);
            t_inertial_start = par.tend-2.*pi./par.fcor.*nip_start;
            t_inertial_end = par.tend-2.*pi./par.fcor.*nip_end;
            fprintf(['Averaging over the inertial '...
               'period(s): %g to %g s\n'],t_inertial_start,t_inertial_end);
            % average window, start and end indices
            its = obj.getTimeIndex('second',t_inertial_start);
            ite = obj.getTimeIndex('second',t_inertial_end);
        else
            error(['Please pass in the average window, in the format of '...
                'either {''hour''/''second'',[i_start, i_end]} '...
                'or {''inertial'',[i,j]} (averaging over the last jth to '...
                'ith inertial period(s)).']);
        end
    end % function parsingAvgWindow

    %% getAnisotropicBarycentricCoord
    function [c,z] = getAnisotropicBarycentricCoord(obj,varargin)
    % getAnisotropicBarycentricCoord returns the coordinates of the
    %   anisotropic barycentric map at each level
    %   arguments are for getProfiles()

        pf = obj.getProfiles(varargin{:});
        c = zeros(obj.nnz-1,3);
        for k=1:obj.nnz-1
            km1=k-1;
            u1u1 = pf.ups(k);
            u2u2 = pf.vps(k);
            u1u2 = pf.uvle(k);
            if km1==0
                u3u3 = 0.5.*pf.wps(k);
                u1u3 = 0.5.*pf.uwle(k);
                u2u3 = 0.5.*pf.vwle(k);
            else
                u3u3 = 0.5.*(pf.wps(k)+pf.wps(km1));
                u1u3 = 0.5.*(pf.uwle(k)+pf.uwle(km1));
                u2u3 = 0.5.*(pf.vwle(k)+pf.vwle(km1));
            end
            if any(isnan([u1u1, u2u2, u3u3, u1u2, u1u3, u2u3]))
                c(k,:) = NaN;
            else
                a = obj.anisotropy_tensor(u1u1,u2u2,u3u3,u1u2,u1u3,u2u3);
                c(k,:) = obj.barycentric_coord(a);
                obj.check_anisotropy_tensor(a);
            end
        end
        z = pf.z_u(1:obj.nnz-1);
    end % function getAnisotropicBarycentricCoord

    %% getAnisotropicTensor
    function [a,z] = getAnisotropicTensor(obj,varargin)
    % getAnisotropicTensor returns the anisotropic tensor
    %   at each level
    %   arguments are for getProfiles()

        nArgs = numel(varargin);
        if nArgs == 1
            pf = varargin{1};
        elseif nArgs == 2
            pf = obj.getProfiles(varargin{:});
        end
        a = zeros(obj.nnz-1,3,3);
        for k=1:obj.nnz-1
            km1=k-1;
            u1u1 = pf.ups(k);
            u2u2 = pf.vps(k);
            u1u2 = pf.uvle(k);
            if km1==0
                u3u3 = 0.5.*pf.wps(k);
                u1u3 = 0.5.*pf.uwle(k);
                u2u3 = 0.5.*pf.vwle(k);
            else
                u3u3 = 0.5.*(pf.wps(k)+pf.wps(km1));
                u1u3 = 0.5.*(pf.uwle(k)+pf.uwle(km1));
                u2u3 = 0.5.*(pf.vwle(k)+pf.vwle(km1));
            end
            if any(isnan([u1u1, u2u2, u3u3, u1u2, u1u3, u2u3]))
                a(k,:,:) = NaN;
            else
                a(k,:,:) = obj.anisotropy_tensor(u1u1,u2u2,u3u3,...
                    u1u2,u1u3,u2u3);
            end
        end
        z = pf.z_u(1:obj.nnz-1);
    end % getAnisotropicTensor

    %% check_anisotropy_tensor
    function r = check_anisotropy_tensor(obj,a)
    % check_anisotropy_tensor(a) checks the properties of I, II, and III
    %   following (21) of Banerjee et al., 2007

        eps = 1e-2;
        r = zeros(1,3);
        tmp = eig(a);
        lambda = sort(tmp,'descend');
        [I,II,III] = obj.anisotropy_invariant(a);
        r(1) = I;
        r(2) = (2.*(lambda(1).^2+lambda(1).*lambda(2)+lambda(2).^2)-II)./II;
        r(3) = (-3.*lambda(1).*lambda(2).*(lambda(1)+lambda(2))-III)./III;
        if abs(r(1)) > eps
            error('I = %g \n',r(1));
        elseif abs(r(2)) > eps
            error('(II(lambda)-II)/II = %g \n',r(2));
        elseif abs(r(3)) > eps
            error('(III(lambda)-III)/III = %g \n',r(3));
        end
    end % function check_anisotropy_tensor

    %% kpp_boundary_layer_depth
    function [h, Ri_bulk, R_shr] = kpp_boundary_layer_depth(obj, varargin)

    % kpp_boundary_layer_depth returns the diagnosed boundary layer depth
    %   based on the bulk Richardson number. See Large et al., 1994 for
    %   more detail.

        % count input arguments
        nArg = length(varargin);
        if nArg == 6
            l_langmuir = 'off';
        elseif nArg == 8
            l_langmuir = varargin{7};
            lasl = varargin{8};
        elseif nArg == 9
            l_langmuir = varargin{7};
            lasl = varargin{8};
            us0 = varargin{9};
        elseif nArg == 10
            l_langmuir = varargin{7};
            lasl = varargin{8};
            b0 = varargin{9};
            utau = varargin{10};
        else
            error(['Argument list:'...
                   ' bxym, uxym, vxym, zu, zw, w_s, '...
                   '[l_langmuir, lasl, [us0] / [b0, utau]]']);
        end
        bxym = varargin{1};
        uxym = varargin{2};
        vxym = varargin{3};
        zu = varargin{4};
        zw = varargin{5};
        w_s = varargin{6};

        % get constants
        Ri_crit = obj.kpp_get_constant('Ri_crit');
        surf_layer_ext = obj.kpp_get_constant('surf_layer_ext');
        beta_T = obj.kpp_get_constant('beta_T');
        c_s = obj.kpp_get_constant('c_s');
        kappa = obj.kpp_get_constant('kappa');
        factor = sqrt(-beta_T./c_s./surf_layer_ext)...
                ./Ri_crit./kappa.^2;
        % calculate stratification
        N2 = obj.NSquared(bxym,zu);
        N2(N2<0) = 0;
        N = sqrt(N2);
        Cv = ones(size(N)).*1.7;
        inds = N<0.002;
        Cv(inds) = 2.1-200.*N(inds);
        % calculate the surface layer averaged velocity and buoyancy
        br = obj.kpp_get_ref(zu,zw,bxym,surf_layer_ext);
        ur = obj.kpp_get_ref(zu,zw,uxym,surf_layer_ext);
        vr = obj.kpp_get_ref(zu,zw,vxym,surf_layer_ext);
        % calculate the unresolved shear squared Utsq
        switch l_langmuir
        case 'off'
            Utsq = -zu.*Cv.*N.*w_s.*factor;
        case 'LW16'
            % Li et al., 2016, (13)
            ef = obj.eFactorVR12(lasl,'LaSLProj',0);
            Utsq = -zu.*Cv.*N.*w_s.*factor.*ef;
        case 'LW16_En'
            % Li et al., 2016, (13)
            Us0sq = us0.^2;
            ef = obj.eFactorVR12(lasl,'LaSLProj',0);
            Utsq = -zu.*Cv.*N.*w_s.*factor.*ef+Us0sq;
        case 'RW16'
            % Reichl et al., 2016, (38)-(39)
            ef = 1+2.3*lasl.^(-0.5);
            Ri_crit = 0.235;
            factor = sqrt(-beta_T./c_s./surf_layer_ext)...
                ./Ri_crit./kappa.^2;
            Utsq = -zu.*Cv.*N.*w_s.*factor.*ef;
        case 'LF17'
            cb1 = 0.17;
            cb2 = 0.15;
            cb3 = 0.083;
            bfactor = cb2.*b0.*(-zu)+cb1.*utau.^3+cb3.*lasl.^(-2).*utau.^3;
            Utsq = -zu.*Cv.*N./Ri_crit.*sqrt(bfactor./w_s);
        case 'LF17ck'
            cb1 = 0.17;
            cb2 = 0.15;
            cb3 = 0.0;
            bfactor = cb2.*b0.*(-zu)+cb1.*utau.^3+cb3.*lasl.^(-2).*utau.^3;
            Utsq = -zu.*Cv.*N./Ri_crit.*sqrt(bfactor./w_s);
        otherwise
            error('Invalid Langmuir option %s',l_langmuir);
        end
        scaling = 1-0.5*surf_layer_ext;
        % the bulk Richardson number
        Ri_bulk = -zu.*scaling.*(br-bxym)./...
                ((ur-uxym).^2+(vr-vxym).^2+Utsq);
        R_ru = ((ur-uxym).^2+(vr-vxym).^2)./Utsq;
%         R_ru = (ur-uxym).^2+(vr-vxym).^2;
%        h = obj.getThresholdDepth(Ri_bulk, zu, Ri_crit);
        [h, R_shr] = obj.getThresholdDepthM(Ri_bulk, zu, Ri_crit, R_ru);
    end

    %% kpp_turbulent_scale
    function w_sm = kpp_turbulent_scale(obj, sigma, hb, bf, utau, sm)
    % kpp_turbulent_scale returns the turbulent velocity scale

        kappa = obj.kpp_get_constant('kappa');
        surf_ext = obj.kpp_get_constant('surf_layer_ext');
        switch sm
            case 's'
                con_sm = 'c_s';
            case 'm'
                con_sm = 'c_m';
            otherwise
                error('Invalid case, should be s or m');
        end
        c_sm = obj.kpp_get_constant(con_sm);
        w_sm = zeros(size(sigma));
        if utau ~= 0
            zeta = min(surf_ext, sigma).*hb.*bf.*kappa./utau.^3;
            w_sm = kappa.*utau./obj.kpp_phi(zeta,sm);
        else
            if bf < 0
                w_sm = kappa...
                    .*(-c_sm.*kappa.*min(surf_ext, sigma).*hb.*bf).^(1/3);
            end
        end

    end % function kpp_turbulent_scale

    %% kpp_phi
    function phi = kpp_phi(obj,zeta, sm)
    % kpp_phi returns the similarity function used in KPP

        zeta_s = obj.kpp_get_constant('zeta_s');
        zeta_m = obj.kpp_get_constant('zeta_m');
        a_s = obj.kpp_get_constant('a_s');
        a_m = obj.kpp_get_constant('a_m');
        c_s = obj.kpp_get_constant('c_s');
        c_m = obj.kpp_get_constant('c_m');
        switch sm
            case 's'
                if zeta >= 0
                    phi = 1 + 5.*zeta;
                elseif zeta >= zeta_s
                    phi = (1 - 16.*zeta).^(-0.5);
                else
                    phi = (a_s - c_s.*zeta).^(-1/3);
                end
            case 'm'
                if zeta >= 0
                    phi = 1 + 5.*zeta;
                elseif zeta >= zeta_m
                    phi = (1 - 16.*zeta).^(-0.25);
                else
                    phi = (a_m - c_m.*zeta).^(-1/3);
                end
            otherwise
                error('Invalid case, should be s or m');
        end
    end % function kpp_phi
    %%
end % methods

%% Static method
methods(Static)
    %% anisotropy_tensor
    function a = anisotropy_tensor(u1u1,u2u2,u3u3,u1u2,u1u3,u2u3)
    % anisotropy_tensor(u1u1,u2u2,u3u3,u1u2,u1u3,u2u3) returns a 3x3
    %   matrix a_{ij} measuring the anisotropy of the turbulence.
    %   a_{ij}=\frac{\overline{u_i u_j}}{2TKE}-\delta_{ij}/3,
    %   TKE=\frac{\overline{u_i u_i}}{2}
    %   (See Lumley, 1979)

        a = zeros(3,3);
        tke = 0.5.*(u1u1+u2u2+u3u3);
        a(1,1) = 0.5.*u1u1./tke-1./3;
        a(2,2) = 0.5.*u2u2./tke-1./3;
        a(3,3) = 0.5.*u3u3./tke-1./3;
        a(1,2) = 0.5.*u1u2./tke;
        a(1,3) = 0.5.*u1u3./tke;
        a(2,1) = 0.5.*u1u2./tke;
        a(2,3) = 0.5.*u2u3./tke;
        a(3,1) = 0.5.*u1u3./tke;
        a(3,2) = 0.5.*u2u3./tke;
    end % function anisotropy_tensor

    %% barycentric_coord
    function c = barycentric_coord(a)
    % barycentric_coord returns the coordinates of the barycentric map
    %   given the anisotropy tensor a_{ij}
    %   (See Banerjee et al., 2007)

        c = zeros(1,3);
        tmp = eig(a);
        lambda = sort(tmp,'descend');
        c(1) = lambda(1)-lambda(2);
        c(2) = 2.*(lambda(2)-lambda(3));
        c(3) = 3.*lambda(3)+1;
    end % function barycentric_coord

    %% anisotropy_invariant
    function [I,II,III] = anisotropy_invariant(a)
    % anisotropy_invariant(a) returns the three anisotropy invariants
    %   for the anisotropy tensor a (See Lumbey, 1979)

        I = trace(a);
        II = 0;
        III = 0;
        for i=1:3
            for j=1:3
                II = II + a(i,j).*a(j,i);
                for n=1:3
                    III = III + a(i,j).*a(i,n).*a(j,n);
                end
            end
        end
    end % function anisotropy_invariant

    %% kpp_get_constant
    function var = kpp_get_constant(vname)
    % kpp_get_constant returns the value of the kpp constant vname
    %   See Large et al., 1994 for more detail

        switch vname
            case 'kappa'
                var = 0.4;
            case 'surf_layer_ext'
                var = 0.1;
            case 'beta_T'
                var = -0.2;
            case 'Ri_crit'
                var = 0.3;
            case 'zeta_s'
                var = -1.0;
            case 'zeta_m'
                var = -0.2;
            case 'a_s'
                var = -sqrt(17)*7; % ~-28.86
            case 'a_m'
                var = 4.2^(-0.25)*1.8; % ~1.26;
            case 'c_s'
                var = sqrt(17)*24; % ~98.96
            case 'c_m'
                var = 4.2^(-0.25)*12; % ~8.38;
            otherwise
                error('Constant %s not found', vname);
        end
    end % function kpp_get_constant

    %% kpp_get_ref
    function vr = kpp_get_ref(zu, zw, var, surf_layer_ext)
    % kpp_get_ref returns an array with same size of zu, containing the
    %   surface layer averaged values, assuming a boundary layer of zu(i)

        zr = zu.*surf_layer_ext;
        dsize = size(zu);
        dz = zeros(dsize);
        dz(2:end) = zw(1:end-1)-zw(2:end);
        dz(1) = 0-zw(1);
        nz = numel(zu);
        vr = zeros(dsize);
        vr(1) = var(1);
        for i=2:nz
            inds = sum(zw > zr(i));
            if inds > 0
                tmp = sum(var(1:inds).*dz(1:inds))...
                    + var(inds+1)*(zw(inds)-zr(i));
                vr(i) = -tmp/zr(i);
            else
                vr(i) = var(1);
            end
        end
    end % function kpp_get_ref

    %% NSquared
    function N2 = NSquared(bxym,zu)
    % NSqured returns N2 from buoyancy

        nb = numel(bxym);
        nz = numel(zu);
        if nb~=nz
            error('Dimensions of bxym and zu do not match');
        end
        N2 = zeros(size(bxym));
        N2_tmp = (bxym(2:end)-bxym(1:end-1))./...
                (zu(2:end)-zu(1:end-1));
%         N2(2:end-1) = 0.5.*(N2_tmp(1:end-1)+N2_tmp(2:end));
        N2(2:end-1) = N2_tmp(2:end);
        N2(1) = N2(2);
        N2(end) = N2(end-1);
    end % function NSquared

    %% getThresholdDepth
    function dat = getThresholdDepth(val, zz, val0)
    % getThresholdDepth returns the shallowest depth at which val profile
    %   reaches the threshold value val0. Linear interpolation.

        if val0 > max(val) || val0 < min(val)
            error('Threshold outside the range of valid values');
        end

        diff = val - val0;
        [~, ind] = min(abs(diff));
        if diff(ind)*diff(ind-1) >= 0
            ind2 = ind + 1;
        else
            ind2 = ind - 1;
        end
        wgt2 = abs(diff(ind))./abs(val(ind)-val(ind2));
        wgt = 1 - wgt2;
        dat = wgt.*zz(ind)+wgt2.*zz(ind2);
    end % function getThresholdDepth

    %% getThresholdDepthM
    function [dat, dat2] = getThresholdDepthM(val, zz, val0, val2)
    % getThresholdDepthM returns the shallowest depth at which val profile
    %   reaches the threshold value val0. Linear interpolation.
    %   Also return the value of val2 at the same depth

        if val0 > max(val) || val0 < min(val)
            error('Threshold outside the range of valid values');
        end

        diff = val - val0;
        [~, ind] = min(abs(diff));
        if diff(ind)*diff(ind-1) >= 0
            ind2 = ind + 1;
        else
            ind2 = ind - 1;
        end
        wgt2 = abs(diff(ind))./abs(val(ind)-val(ind2));
        wgt = 1 - wgt2;
        dat = wgt.*zz(ind)+wgt2.*zz(ind2);
        dat2 = wgt.*val2(ind)+wgt2.*val2(ind2);
    end % function getThresholdDepthM

    %% getUsSL
    function dat = getUsSL(stokes, zz, h)
    % getUsSL returns the surface layer (defined as h/5) averaged Stokes
    %   drift given the Stokes drift profile. The discrete Stokes drift
    %   profile is first interpolated to a finer vertical grid using
    %   spline interpolation, then averaged over the surface layer.

        z0 = zz(1);
        z1 = zz(end);
        dz = (zz(2)-zz(1))/10;
        znew = z0:dz:z1;
        usnew = interp1(zz,stokes,znew,'spline');
        [~,ind] = min(abs(znew+h/5));
        if ind == 1
            ustot = abs(z0).*stokes(1);
        else
            ustot = abs(z0).*stokes(1)-trapz(znew(1:ind),usnew(1:ind));
        end
        dat = ustot./h.*5;
    end % function getUsSL

    %% getHAvgPfl
    function dat = getHAvgPfl(pfl, zz, h)
    % getHAvgPfl returns the boundary layer averaged data given
    %   the profile. The discrete profile is first interpolated
    %   to a finer vertical grid using spline interpolation.

        % remove NaNs
        tmp = pfl;
        pfl(isnan(tmp)) = [];
        zz(isnan(tmp)) = [];
        z0 = zz(1);
        z1 = zz(end);
        dz = (zz(2)-zz(1))/10;
        znew = z0:dz:z1;
        pflnew = interp1(zz,pfl,znew,'spline');
        [~,ind] = min(abs(znew+abs(h)));
        dat = (abs(z0).*pfl(1)-trapz(znew(1:ind),pflnew(1:ind)))./h;
    end % function getHAvgPfl

    %% getVInt
    function dat = getVInt(pfl, zz)
    % getVInt returns the vertical integrated profile data, given the
    %   profile and veritcal coordinate.

        % remove nan
        tmp = pfl;
        pfl(isnan(tmp)) = [];
        zz(isnan(tmp)) = [];
        z0 = zz(1);
        dat = abs(z0).*pfl(1)-trapz(zz,pfl);
    end

    %% getVIntH
    function dat = getVIntH(pfl, zz, h)
    % getVInt returns the vertical integrated profile data below a certian
    %   depth, given the profile, veritcal coordinate, and the depth.
    %   The discrete profile data is first interpolated
    %   to a finer vertical grid using spline interpolation.

        % remove nan
        tmp = pfl;
        pfl(isnan(tmp)) = [];
        zz(isnan(tmp)) = [];
        z0 = zz(1);
        z1 = zz(end);
        dz = (zz(2)-zz(1))/10;
        znew = z0:dz:z1;
        pflnew = interp1(zz,pfl,znew,'spline');
        if h == 0
            dat = abs(z0).*pfl(1)-trapz(znew,pflnew);
        else
            [~,inds] = min(abs(znew+h));
            dat = -trapz(znew(inds:end),pflnew(inds:end));
        end
    end

    %% refineZ
    function dat = refineZ(zz, nr)
    % refineZ returns the refined grid with dz_new = 1/nr * dz
        z0 = zz(1);
        z1 = zz(end);
        dz = (zz(2)-zz(1))/nr;
        dat = z0:dz:z1;
    end

    %% getStokesDepth
    function dat = getStokesDepth(stokes, zz)
    % getDs returns the Stokes depth given the discrete Stokes
    %   drift profile.

        z0 = zz(1);
        z1 = zz(end);
        dz = (zz(2)-zz(1))/10;
        znew = z0:dz:z1;
        usnew = interp1(zz,stokes,znew,'spline');
        dat = (abs(z0).*stokes(1)-trapz(znew,usnew))./stokes(1);
    end % function getStokesDepth

    %% getPosDepth
    function dat = getPosDepth(val, zz)
    % getPosDepth returns the depth at which val profile shift
    %   from negative to positive. Linear interpolation.

        nsize = numel(zz);
        l_neg1 = 0;
        for i=2:nsize
            l_neg0 = l_neg1;
            if val(i)<0
                l_neg1 = 1;
            else
                l_neg1 = 0;
            end
            if l_neg0 == 1 && l_neg1 == 0
                ind = i;
                ind2 = i+1;
                break;
            end
        end

        wgt2 = abs(val(ind))./abs(val(ind)-val(ind2));
        wgt = 1 - wgt2;
        dat = wgt.*zz(ind)+wgt2.*zz(ind2);
    end % function getPosDepth

    %% eFactorVR12
    function ef = eFactorVR12(la,casename,alpha)
    % eFactorVR12 returns the enhancement factor from scaling in
    %   Van Roekel et al., 2012
    %   ef = eFactorVR12(la, casename, alpha)
    %   Input:
    %           la: Langmuir number
    %           casename: 'LaTurb' or 'LaSLProj'
    %           alpha: angle between wind and LCs (optional if 'LaTurb')

        funcname = 'eFactorVR12';
        if (strcmp(casename,'LaTurb'))
            ef = sqrt(1+(3.1.*la).^(-2)+(5.4.*la).^(-4));
        elseif (strcmp(casename,'LaSLProj'))
            if (nargin==3)
                ef = cos(alpha).*sqrt(1+(1.5.*la).^(-2)+(5.4.*la).^(-4));
            else
                error([funcname ': not enough input arguments...']);
            end
        else
            error([funcname ': unknown case name...']);
        end
    end  % function eFactorVR12

    %% getAlphaL
    function alpha = getAlphaL(u, v, us, vs, zz, hx)
    % getAlphaL returns the angle in degrees between LC and wind
    %   diagnosed from the detph-averaged Lagrangian shear

        z0 = zz(1);
        z1 = zz(end);
        dz = (zz(2)-zz(1))/10;
        znew = z0:dz:z1;
        unew = interp1(zz, u, znew, 'linear');
        vnew = interp1(zz, v, znew, 'linear');
        usnew = interp1(zz, us, znew, 'linear');
        vsnew = interp1(zz, vs, znew, 'linear');
        ushr = (unew(1:end-1)-unew(2:end))./(znew(1:end-1)-znew(2:end));
        vshr = (vnew(1:end-1)-vnew(2:end))./(znew(1:end-1)-znew(2:end));
        usshr = (usnew(1:end-1)-usnew(2:end))./(znew(1:end-1)-znew(2:end));
        vsshr = (vsnew(1:end-1)-vsnew(2:end))./(znew(1:end-1)-znew(2:end));
        dl = 0.2.*hx;
        [~, inds] = min(abs(znew+dl));
        ratioshr = mean(vshr(1:inds)+vsshr(1:inds))...
                ./ mean(ushr(1:inds)+usshr(1:inds));
        alpha = atand(ratioshr);

    end % getAlphaL

    %% getAlphaLOW
    function alpha = getAlphaLOW(utau, us0, hx, theta, z1)
    % getAlphaLOW returns the angle in degrees between LC and wind
    %   estimated from the law of the wall

        kappa = 0.4;
        num = sind(theta);
        den = utau./us0./kappa.*log(abs(hx./z1))+cosd(theta);
        alpha = atand(num./den);
    end

    %% zAvg
    function val = zAvg(z,dat,z1,z2,gtype)
    % zAvg averages dat on grid z over a layer from z1 to z2
    %   It is required that z1<z2
    %   If numel(z1)>1, and numel(z1)==size(dat)()

        % find dz
        nz = length(z);
        nd = size(dat);
        nn = length(nd);
        % check dimensions
        if numel(z)~=nz
            error('z should be one-dimensional');
        elseif nz~=nd(1)
            error(['The first dimension of dat should be the same '...
                   'with z']);
        elseif nn>3
            error('dat should be no larger than 3-dimensional');
        end
        z = reshape(z,[1,nz]);
        dz0 = z(1:end-1)-z(2:end);
        if strcmp(gtype,'w')
            dz = zeros([1,nz]);
            dz(1) = 0-z(1);
            dz(2:end) = dz0;
        elseif strcmp(gtype,'u')
            dz1 = zeros([1,nz]);
            dz2 = zeros([1,nz]);
            dz1(1) = 0-z(1);
            dz1(2:end) = 0.5.*dz0(:);
            dz2(1:end-1)  = 0.5.*dz0;
            dz2(end) = 0.5.*dz0(end);
            dz = dz1 + dz2;
        end
        nz1 = numel(z1);
        if nz1==1
            [~,ind1] = min(abs(z-z1));
            [~,ind2] = min(abs(z-z2));
            id1 = min(ind1,ind2);
            id2 = max(ind1,ind2);
            if nn==2
                val = dz(1,id1:id2)*dat(id1:id2,:)./sum(dz(id1:id2));
            elseif nn==3
                val = zeros([nd(2),nd(3)]);
                for i=1:nd(2)
                    val(i,:) = dz(1,id1:id2)*squeeze(dat(id1:id2,i,:))...
                             ./sum(dz(id1:id2));
                end
            end
        elseif nn==2 && nz1==nd(2)
            val = zeros([1,nz1]);
            for i=1:nz1
                [~,ind1] = min(abs(z-z1(i)));
                [~,ind2] = min(abs(z-z2));
                id1 = min(ind1,ind2);
                id2 = max(ind1,ind2);
                val(i) = dz(1,id1:id2)*dat(id1:id2,i)./sum(dz(id1:id2));
            end
        elseif nn==3 && nz1==nd(3)
            val = zeros([nd(2),nz1]);
            for i=1:nz1
                [~,ind1] = min(abs(z-z1(i)));
                [~,ind2] = min(abs(z-z2));
                id1 = min(ind1,ind2);
                id2 = max(ind1,ind2);
                val(:,i) = dz(1,id1:id2)*dat(id1:id2,:,i)./sum(dz(id1:id2));
            end
        else
            error('The dimension of z1 should match the time dimension');
        end
    end % function zAvg
%%
end % methods

end % classdef


