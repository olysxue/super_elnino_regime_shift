%   stars.m
%
%   Sequential T-test for Analysis of Regime Shifts (change point detection in
%   mean level).  This is adapted from VBA code provided in Sergei Rodionov's
%   Microsoft Excel plugin.
%
%   Usage:
%       [R, opts] = stars(x, varargin);
%
%       R       -   rsi index values; non-zero indicate regime shift points
%       opts    -   structure with STARS parameters used in this analysis
%
%   The varargin term is used to pass in STARS arguments as keyword, value
%   pairs.  Valid keywords are:
%
%       'L'      -   length of analysis window (default = 10)
%       'p'      -   significance, e.g., 0.05 for 95% confidence (default = 0.1)
%       'h'      -   Huber weight parameter for outliers (default = 1)
%       'sigL'   -   L-pt variance estimate (default: computed from x)
%
%   Example:
%       > x = rand(100, 1);                 % baseline [0,1] noise
%       > x(18:35) = 3.5;                   % +step transient
%       > x(62:64) = -7.7;                  % brief -step transient
%       > [rsi, opts] = stars(x, 'L', 10, 'p', 0.05, 'h', 1);
%
%   References:
%       Rodionov, S. N. (2004), "A sequential algorithm for testing climate
%           regime shifts," Geophys. Res. Lett., 21, L09204.
%
%       Regime detection plugin v. 2.1 for Microsoft Excel, 
%           http://www.beringclimate.noaa.gov/regimes/
%
%   See Also:
%       starsMeans, starsDemo_002
%
%   Andrew B. Gardner
%   <andywocky@gmail.com>
%   Feb 21, 2006
%
%   updates:
%       defaultParms changed from cell to structure in response to Matlab v7 struct warnings
%       minor syntax changes
%   question: is the Huber weight parameter used?
%   Andrew R. Friedman <andfried@gmail.com>
%   04 August 2013
%
% -------------------------------------------------------------------------
%	NOTES
%
%	1.	Verified against the excel plugin for two, 110-pt sequences
%			(abg, 5/25/2006)
%	2.	There are significant differences between this code and the original
%			stars paper, e.g., (i) Huber weighting, (ii) the use of a running
%			mean, etc.  (abg, 5/26/2006)
%	3.	Here is a quick lookup-table mapping variable names between this
%			m-file and Sergei's VB code:
%
%					m-file			VB-code
%					------			-------
%					  N				rng.Count
%					  p				prob
%					  L				cutoff
%					 sigL			dblAveVar
%					 cp				ChangePoint
%					 lvl			dblRegimeMean
%					thresh			dblDiff
%					  x				rngX
%					  h				frmMain.txtHuberParam.Value
%
% -------------------------------------------------------------------------
%
%
function [R, opts] = stars(x, varargin)
    %
    %   Input parsing
    %
    defaultParms = struct('L', 10, 'p', 0.1, 'h', 1);
    sOpts = cell2struct(varargin(2:2:end), varargin(1:2:end-1), 2);
    opts = mergeStruct(defaultParms, sOpts);

    L = opts.L;
    p = opts.p;
    h = opts.h;

    %
	%	Initialization
	%
    %   variance can be passed-in, or we can compute
    if ~isfield(opts, 'sigL')
        opts.sigL = estimateSigma(L);           % sample L-pt variance
    end
    sigL = opts.sigL;
    thresh = getThreshold(L, p, sigL);          % critical threshold
    lvl = hWeightedAverage(x(1:L));             % initial mean level
    R = zeros(size(x));                         % rsi values
    cp = 1;                                     % current change-point index
    N = length(x);                              % number of samples

    %
	%	Main routine: iterate over the time series
	%
    for k = 2:N;
        R(k) = rsi(k);

        %   too few samples to confirm last regime change?
        if (abs(R(k)) > 0) && (k > (N-L+1))
            return
        end

        %   test for regime shifts and update current regime mean (unless we
        %   are within L-pts of most recent change-point)
        if R(k) == 0
            if k >= (cp + L)
                lvl = hWeightedAverage(x(cp:k));    % same regime, far enough from cp
            end
        else
            cp = k;                             % regime change
            lvl = hWeightedAverage(x(k:k+L-1)); % same regime, far enough from cp
        end
    end


    % -------------------------------------------------------------------------
    %   rsi(k)
    %
    %   Compute the rsi for a given sample index, regime mean, and critical
    %   threshold.
    % -------------------------------------------------------------------------
    function r = rsi(k);
        if x(k) > (lvl + thresh)
            r = cusumUp(k);
        elseif x(k) < (lvl - thresh)
            r = cusumDown(k);
        else
            r = 0;
        end
    end


    % -------------------------------------------------------------------------
    %   cusumUp()
    %
    %   Compute the L-pt cusum for positive level changes.  For a positive
	%	regime change to be accepted, we require the L-pt lookahead samples
	%	to produce a cusum sequence which does not go negative.
    % -------------------------------------------------------------------------
    function cs = cusumUp(k)
		% LL sets the look ahead length: L, or the number of points until the
		% end of the signal.
        LL = min(L, N-k+1);
		
		% dblXdev is the length-LL vector of normalized deviations of x outside
		% of the range lvl +/- thresh
        dblXdev = (x(k:k+LL-1) - (lvl+thresh)*ones(LL,1)) / sqrt(sigL) + eps;
		
		% these are Huber weight values, so large deviations are deemphasized
		dblXweight = min( [ones(size(dblXdev)), 1./abs(dblXdev)]' )';
		
		% the cusum is the integral of the weighted deviations; we normalize
		% here, too, by dividing by the sum of the weights
        cs = filter(1, [1 -1], dblXweight .* dblXdev ) / sum(dblXweight);
			
		% we check for cusum values below zero, which would indicate a failed
		% regime change; otherwise, we have a positive shift
        if length(find(cs < 0) > 0)
            cs = 0;
        else
            cs = cs(LL);
        end
    end


    % -------------------------------------------------------------------------
    %   cusumDown()
    %
    %   Compute the L-pt cusum for negative level changes.
    % -------------------------------------------------------------------------
    function cs = cusumDown(k);
		% LL sets the look ahead length: L, or the number of points until the
		% end of the signal.
        LL = min(L, N-k+1);
		
		% dblXdev is the length-LL vector of normalized deviations of x outside
		% of the range lvl +/- thresh
        dblXdev = (x(k:k+LL-1) - (lvl-thresh)*ones(LL,1)) / sqrt(sigL) + eps;
		
		% these are Huber weight values, so large deviations are deemphasized
		dblXweight = min( [ones(size(dblXdev)), 1./abs(dblXdev)]' )';
		
        % the cusum is the integral of the weighted deviations; we normalize
		% here, too, by dividing by the sum of the weights
        cs = filter(1, [1 -1], dblXweight .* dblXdev ) / sum(dblXweight);
		
		% we check for cusum values above zero, which would indicate a failed
		% regime change; otherwise, we have a negative shift
        if length(find(cs > 0) > 0)
            cs = 0;
        else
            cs = cs(LL);
        end
    end


    % -------------------------------------------------------------------------
    %   getThreshold()
    %
    %   Calculate the critical threshold of deviation that signals regime
    %   changes.  This does not change over the signal.
    % -------------------------------------------------------------------------
    function thresh = getThreshold(L, p, sigL)
        dof = 2*L-2;                        % number degrees freedom
        t = tinv(1 -p/2, dof);              % crit 2-sided t-value
        thresh = t*sqrt(2*sigL/L);          % crit deviation
    end


    % -------------------------------------------------------------------------
    %   estimateSigma(x,L)
    %
    %   Estimate the long-term, L-pt variance (assume homoskedastic signals).
    %   Improvements would be to implement a state machine that for each
    %   sample, looks at the L-pt window with probability p s.t. the mean
    %   number of looks gives the correct estimate.
    %
    %   Could be important to estimate var from past samples, only-- this is
    %   NOT THE CASE IN THIS FUNCTION.
    % -------------------------------------------------------------------------
    function sigL = estimateSigma(L)
        %   Estimate the long-term length-L variance.  If the signal is >> the
        %   length of the analysis window, sample to estimate variance.
        nx = length(x);
        if (nx/L) > 300;
            ix = round(rand(1,100)*(nx-2*L)+L);
        else
            ix = (L:nx);
        end

        s = 0;
        for i = 1:length(ix);
            xwin = x(ix(i)-L+1:ix(i));
            s = s + var(xwin);
        end
        sigL = s / length(ix);
    end


    % -------------------------------------------------------------------------
    %   hWeightedAverage(xwin)
    %
    %   Calculates the mean estimate for a given range using Huber's weights.
    % -------------------------------------------------------------------------
    function dblWeightedAve = hWeightedAverage(xwin)
        % simple estimate of the regime mean for the windowed clip
        dblEstAve = mean(xwin);

        % The procedure consists of two iterations.  Initially, the estimate of the
        % regime average is the simple unweighted arithmetic mean.  In the second
        % iteration it is weighed average from the first iteration.
        for k = 1:2;
            dblSumofWeights = 0;
            dblSumAve = 0;

            dblXdev = (xwin-dblEstAve)/sqrt(sigL)+eps;
            dblXweight = min( [ones(size(dblXdev)), 1./abs(dblXdev)]' )';

            dblSumofWeights = sum(dblXweight);
            dblSumAve = dblXweight' * dblXdev;

            dblSumAve = dblSumAve / dblSumofWeights;
            dblSumAve = dblSumAve * sqrt(sigL) + dblEstAve;
            dblEstAve = dblSumAve;
        end

        % our return value
        dblWeightedAve = dblEstAve;
    end
end
