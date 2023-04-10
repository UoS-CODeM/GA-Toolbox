function [mpgastate, mpgaoptions] = mpga(mpgaoptions, varargin)
% implements a multi-population genetic algorithm optimisation
% 
% Syntax
%
% [mpgastate, mpgaoptions] = mpga(mpgaoptions)
% [...] = mpga(..., 'Parameter', value) 
%
% Description
%
% implements the multi-population GA derived from the GA Toolbox produced
% by the Department of Automatic Control and Systems Engineering at the
% University of Sheffiedl. See the GA Toolbox manual for details of the
% algorithms employed. The mpga minimises an objective function.
%
% Input:
%
%  mpgaoptions - a structure containing some or all of the following
%   members:
%
%   OBJ_F : string or function handle to the objective function
%
%   XOV_F : string or function handle to the recombination function for 
%    individuals
%
%   MUT_F : string or function handle to the mutation function
%
%   SEL_F : string or function handle to the selection function
%
%   GGAP : (optional scalar) Generation gap, how many new individuals are 
%    created, default is 0.8
%
%   INSR : (optionsl scalar) Insertion rate, how many of the offspring are
%    inserted, default is 0.9
%
%   XOVR : (optional scalar) Crossover rate, default is 1
%
%   SP : (optional scalar) Selective Pressure, default is 2
%
%   MUTR : (optional scalar) Mutation rate, default value is calculated
%    based on the number of variables in the objective function
%
%   MIGR : (optional scalar) Migration rate between subpopulations,
%    default is 0.1
%
%   MAXGEN : (optional scalar) Max number of generations, default is
%    calculated from the number of variables
%
%   MIGGEN : (optional scalar) Number of generations between migration
%    (isolation time), default is calculated as ceil(MAXGEN / 10)
%        
%   TERMEXACT : (optional scalar) Value for termination if minimum
%    reached, default is 1e-8
%
%   TERMIMPROVEGENS : (optional scalar) Value for number of generation
%    where there has been an improvement less than TERMIMPROVETOL after
%    which the optimisation is terminated. Default is inf meaning the full
%    number of generations specified in MAXGEN will be performed.
%
%   TERMIMPROVETOL: (optional scalar) Value for termination if the
%    improvement in the best score is less than this value for the number
%    of generations specified in TERMIMPROVEGENS, default is 1e-8.
%
%
%   SUBPOP : (optional scalar) Number of subpopulations, default is
%    calculated based on the number of variables
%
%   NIND : (optional scalar) Number of individuals per subpopulations,
%    default is calculated from the number of variables
%
%   MUTR : (scalar) Mutation rate depending on NVAR
%
%   STEP : (optional scalar) The variables of the best individuals will be 
%    plotted for each generation sparated by an interval of size STEP.
%    Default is set to ceil(MAXGEN / 10)
%
%   DISPLAYMODE : (optional scalar) If 1, some graphs are displayed at
%    each interval in STEP, if 2, some text output is displayed at each
%    interval, but no graphs are plotted except the final results. If 3,
%    same as 2, but no final results are plotted. Otherwise, no output is
%    displayed
%
%   SAVEMODE : (optional scalar) If 0, no data is saved, otherwise, if no 
%    filename is supplied data is perioically saved in a default file name:
%    "OBJ_F" + _output.mat which can be used to resume the algorithm from
%    the current generation. If a filename is supplied, the data is saved
%    in that file.
%
%   FILENAME : optional name of a file in which to periodically save data 
%    about the progress of the  GA, if not supplied, and SAVEMODE is not
%    zero, a default file name will be used. The default filename will be
%    the name of the objective function with the string '_output.mat'
%    appended.
%
%   RESUMEFILE : name of file where data from a previous run can be found. 
%    If supplied, the GA will resume from this previous starting point
%
%   EVALUATE : true or false flag, the objective function will be
%    evaluated and returned with all output arguments if true. Defaults to
%    false.
%
%   EVALINDS : if EVALUATE is true, this is an optional vector of
%    indices in the latest Chomosome to evaluate. Set this value the the
%    string ':' to evaluate the entire chromosome.
%
%
% In addition to the mpgaoptions structure you may also supply optional
% additional input arguments to the objective function, selection function
% and recombination function. These optional arguments must be stored in a
% cell array, and are supplied to mpga as parmater-value pairs (i.e. a
% string and the argument, see example below). The posssible additional
% arguments are:
%
%  'ObjectiveArgs' - cell array of additional argumetns to be passed to the
%    objective function (in mpgaoptions.OBJ_F). It will then be evaluated
%    like so:
%
%    scores = feval (mpgaoptions.OBJ_F, objargs{:})
%
%    where objargs is the cell array of additional input arguements
%
%  'SelectionArgs' - cell array of additional argumetns to be passed to the
%    selection function (in mpgaoptions.SEL_F). It will then be evaluated
%    like so: 
%
%    feval(SEL_F, [], NSel, selargs{:})
%
%    where selargs is the cell array of additional input arguements
%    supplied. For more details on the selection function, see the manual
%    and the default selection functions.
%
%  'RecombinArgs' - cell array of additional argumetns to be passed to the
%    selection function (in mpgaoptions.SEL_F). It will then be evaluated
%    like so: 
%
%    feval(XOV_F, [], RecOpt, xovargs{:})
%
%    where xovargs is the cell array of additional input arguements
%    supplied for the recombination (crossover) function. For more details
%    on the recombination function, see the manual and the default
%    recombination functions.
%
% Example:
%
% The user has an objective function which depends on a some data he does
% not wish to load in and out of memory frequently, and some arguments
% dictating the behaviour of the objective function he would prefer not to
% hard-code into the objective function.
% 
% The user then creates an objective function with the calling syntax
% below:
% 
%      ObjVal = objtestfun(Chrom, rtn_type, someData, objmode)
% 
% To optimise this function with the GA the following script is then
% required (using the minimum number of optional arguments):
% 
%       mpgaoptions.OBJ_F = 'objtestfun'
%       mpgaoptions.SEL_F = 'sus';
%       mpgaoptions.XOV_F = 'recint';
%       mpgaoptions.MUT_F = 'mutbga';
% 
%       someData = load('myData.dat')
% 
%       mode = 0;
%       
%       objfunargs = {someData, mode};
% 
%       [mpgastate, mpgaoptions] = mpga(mpgaoptions, 'ObjectiveArgs', objfunargs)
%
% Additional function arguments can be supplied to the objective, selection
% and recombination functions using the parameters 'ObjectiveArgs',
% 'SelectionArgs' and 'RecombinArgs' respectively. The calling format of
% any functions must match those of existing toolbox functions, other than
% the ability to accept additional arguments.
%
% Output:
% 
% by defualt mpga outputs a structure, mpgastate containing information
% about the GA run, the following fields will be included:
%
%   mpgastate.ScoreTracking - (n X 3) matrix. The first column is the best
%       individual's objective function score at each generation. The
%       second column is the mean value of the Objective function scores
%       for each generation. The thrid column is the number of function
%       iteration performed for each generation.
%
%   mpgastate.BestChroms - (n x NVAR) matrix. The variable values of the best 
%       individual at each generation.
%
%   mpgastate.Chrom -  (n x NVAR) matrix. The current chromosomes, i.e. the 
%       variables of all members of the latest population.
%
%   mpgastate.ObjV - (n x 1) vector. The values of the objective function for
%       the current generation, i.e. the individuals in Chom.
%
%   mpgastate.Generation - scalar value of the last generation number
%       evaluated
%
%   GAoptions - structure containing the mpgaoptions used in the GA
%
% If the evaluate flag is true (see above), a single arguemet is returned
% which is a cell array of cell arrays. Each element of the outer cell
% array contains a cell array with the arguments returned by the objective
% function. The objective funciton must expect the rtn_type value 4 when
% this is being requested, and return the arguments in  cell array.
%
%

% Author: Richard Crozier, The University of Edinburgh 2009
% Modified: 09 Mar 2012
%

    % parse the optional arguments that can be passed to the various GA
    % functions
    mpgastate.ObjectiveArgs = {};
    mpgastate.SelectionArgs = {};
    mpgastate.RecombinArgs = {};
    
    mpgastate = parse_pv_pairs(mpgastate, varargin);
    
    % Validate and complete the mpgaoptions structure
    mpgaoptions = parsempgaoptions(mpgaoptions, mpgastate.ObjectiveArgs);

    if mpgaoptions.EVALUATE
        
        % load the stored mpga progress information
        if exist (mpgaoptions.RESUMEFILE, 'file')
            load(mpgaoptions.RESUMEFILE);
        else
            error('mpga resume file does not appear to exist.')
        end
        
        % check it was a valid GA state file
        if ~exist('mpgastate', 'var')
            error('mpgastate could not be loaded from resume file.')
        end
        
        % now evaluate the objective function
        if isempty (mpgaoptions.EVALINDS)
            objargs = [{mpgastate.BestChroms(mpgastate.Generation,:), 4}, mpgastate.ObjectiveArgs];
        else
            objargs = [{mpgastate.Chrom(options.EVALINDS,:), 4}, mpgastate.ObjectiveArgs];
        end
        mpgastate = feval(mpgaoptions.OBJ_F, objargs{:});
        mpgaoptions = [];
        
        return;
        
    elseif mpgaoptions.RESUME
        % If a file from a previous run is supplied, we resume the previous
        % run by loading the variables into memory, otherwise, an initial
        % population is generated and evaluated
        newMaxGen = mpgaoptions.MAXGEN;
        
        % load the stored mpga progress information
        load(mpgaoptions.RESUMEFILE);
        
        % check it was a valid GA state file
        if ~exist('mpgastate', 'var')
            error('mpgastate could not be loaded from resume file.')
        end
        
        if newMaxGen > mpgaoptions.MAXGEN
            expandvars = newMaxGen - mpgaoptions.MAXGEN;
            mpgaoptions.MAXGEN = newMaxGen;
            mpgastate.ScoreTracking = [ mpgastate.ScoreTracking; NaN * ones(expandvars, 3) ];
            mpgastate.BestChroms = [ mpgastate.BestChroms; ones(expandvars, size(mpgastate.Chrom, 2)) .* NaN];
        end
        
        % re-run mpgaoptions in case any default have been added since mpga
        % state was saved
        mpgaoptions = parsempgaoptions(mpgaoptions, mpgastate.ObjectiveArgs);
        
        fprintf(1, '\nResuming GA from generation %d\n', mpgastate.Generation);
        
        if nargin > 1
            warning('MPGA:ResumeArgs', ...
                'Extra arguments in mpga are being ignored as GA is being resumed from a previous run');
        end
        
        % construct a cell array containing the objective function
        % arguments
        objargs = [{mpgastate.Chrom, []}, mpgastate.ObjectiveArgs];

    else
        % initialise various variables and perform the initial evaluation
        % of the population
        
        % Get boundaries of objective function
        objargs = [{[], 1}, mpgastate.ObjectiveArgs];
        mpgastate.FieldDR = feval(mpgaoptions.OBJ_F, objargs{:});

        % Get value of minimum, defined in objective function
        objargs{2} = 3;
        mpgastate.GlobalMin = feval(mpgaoptions.OBJ_F, objargs{:});

        % Get title of objective function, defined in objective function
        objargs{2} = 2;
        mpgastate.FigTitle = [feval(mpgaoptions.OBJ_F, objargs{:}) '   (' int2str(mpgaoptions.SUBPOP) ':' int2str(mpgaoptions.MAXGEN) ') '];

        % Clear mpgastate.ScoreTracking and storing matrix
        % Preallocate Matrix for storing best results with all NaN values
        mpgastate.ScoreTracking = NaN * ones(mpgaoptions.MAXGEN, 3);
        % set the running count of function evaluations to zero
        mpgastate.ScoreTracking(:,3) = zeros(size(mpgastate.ScoreTracking,1),1);

        % Create real population
        mpgastate.Chrom = crtrp(mpgaoptions.SUBPOP * mpgaoptions.NIND, mpgastate.FieldDR);

         % Matrix for storing best individuals, preallocate with all NaNs
        mpgastate.BestChroms = ones(mpgaoptions.MAXGEN, size(mpgastate.Chrom, 2)) .* NaN;
        
        % reset count variables
        mpgastate.Generation = 0;
        mpgastate.TermOpt = 0;
        mpgastate.GensWithoutImprovement = 0;
    
        % Calculate objective function for entire population. Note that
        % mpga expects this output to be a column vector of the same
        % height as mpgastate.Chrom
        objargs = [{mpgastate.Chrom, []}, mpgastate.ObjectiveArgs];
        mpgastate.ObjV = feval(mpgaoptions.OBJ_F, objargs{:});

        % add first count of number of objective function evaluations
        mpgastate.ScoreTracking(1,3) = mpgastate.ScoreTracking(1,3) + mpgaoptions.NIND;
        
        % Store the best and mean objective values and the best individual
        [mpgastate.ScoreTracking(1,1), ix] = min(mpgastate.ObjV);
        mpgastate.ScoreTracking(1,2) = mean(mpgastate.ObjV);
        % Add the chromosome of the current best individual to the end of
        % the mpgastate.BestChroms matrix
        mpgastate.BestChroms(1,:) = mpgastate.Chrom(ix,:);
        % store the current best individual and their score in mpgastate.BestInd
        mpgastate.BestInd = [mpgastate.Chrom(ix,:), mpgastate.ScoreTracking(1,1)];
        
        if mpgaoptions.SAVEMODE == 1;
            try
                save(mpgaoptions.FILENAME, 'mpgaoptions', 'mpgastate');
            catch
                warning('MPGA:savestatefailed', ...
                    'Attempt to save mpga state file failed at gen %d', ...
                    mpgastate.Generation);
            end
        end

    end

    hfig_iter = [];
    
    % Iterate subpopulation till termination or MAXGEN
    while ((mpgastate.Generation < mpgaoptions.MAXGEN) && (mpgastate.TermOpt == 0))

        % Fitness assignment to whole population
        FitnV = ranking(mpgastate.ObjV, [2 0], mpgaoptions.SUBPOP);

        % Select individuals from population
        SelCh = select(mpgaoptions.SEL_F, mpgastate.Chrom, FitnV, mpgaoptions.GGAP, mpgaoptions.SUBPOP, mpgastate.SelectionArgs);

        % Recombine selected individuals
        SelCh = recombin(mpgaoptions.XOV_F, SelCh, mpgaoptions.XOVR, mpgaoptions.SUBPOP, mpgastate.RecombinArgs);

        % Mutate offspring
        SelCh = mutate(mpgaoptions.MUT_F, SelCh, mpgastate.FieldDR, [mpgaoptions.MUTR], mpgaoptions.SUBPOP);

        % Calculate objective function for offspring
        objargs{1} = SelCh;
        ObjVOff = feval(mpgaoptions.OBJ_F, objargs{:});
        % store the number of function evaluations
        mpgastate.ScoreTracking(mpgastate.Generation+1,3) = mpgastate.ScoreTracking(mpgastate.Generation+1,3) + size(SelCh,1);

        % Insert best offspring in population replacing worst parents
        [mpgastate.Chrom, mpgastate.ObjV] = reins(mpgastate.Chrom, SelCh, mpgaoptions.SUBPOP, [1 mpgaoptions.INSR], mpgastate.ObjV, ObjVOff);
        
        % Store the best and average objective values and the best individual
        [mpgastate.ScoreTracking(mpgastate.Generation+1,1),ix] = min(mpgastate.ObjV);
        if mpgastate.Generation > 0
            if (mpgastate.ScoreTracking(mpgastate.Generation+1,1) < mpgastate.ScoreTracking(mpgastate.Generation,1)) ...
                    && abs(mpgastate.ScoreTracking(mpgastate.Generation,1) - mpgastate.ScoreTracking(mpgastate.Generation+1,1)) > mpgaoptions.TERMIMPROVETOL ...
                mpgastate.GensWithoutImprovement = 0;
            else
                mpgastate.GensWithoutImprovement = mpgastate.GensWithoutImprovement + 1;
            end
        end
        mpgastate.ScoreTracking(mpgastate.Generation+1,2) = mean(mpgastate.ObjV);
        mpgastate.BestChroms(mpgastate.Generation+1,:) = mpgastate.Chrom(ix,:);
        mpgastate.BestInd = [mpgastate.Chrom(ix,:), mpgastate.ScoreTracking(mpgastate.Generation+1,1)];
        
        % Increment the generation number
        mpgastate.Generation = mpgastate.Generation + 1;
        
        % save the current state of the GA, if desired
        if mpgaoptions.SAVEMODE == 1
            try
                save(mpgaoptions.FILENAME, 'mpgaoptions', 'mpgastate');
            catch
                warning('MPGA:savestatefailed', ...
                    'Attempt to save mpga state file failed at gen %d', ...
                    mpgastate.Generation);
            end
        end
        
        % run an output function if supplied
        if ~isempty(mpgaoptions.OUTPUTFCN)
            feval(mpgaoptions.OUTPUTFCN, mpgastate, mpgaoptions, mpgaoptions.OUTPUTFCNARGS{:})
        end

        % process any outputs the user has requested during the GA
        % evalaution
        hfig_iter = iteroutput(mpgastate, mpgaoptions, hfig_iter);

        % Check, if best objective value near mpgastate.GlobalMin ->
        % termination criterion and compute difference between
        % mpgastate.GlobalMin and best objective value
        ActualMin = abs(min(mpgastate.ObjV) - mpgastate.GlobalMin);
        
        % if ActualMin smaller than TERMEXACT --> termination
        if ((ActualMin < (mpgaoptions.TERMEXACT * abs(mpgastate.GlobalMin))) || (ActualMin < mpgaoptions.TERMEXACT))
            mpgastate.TermOpt = 1;
        end
        
        if mpgastate.GensWithoutImprovement > mpgaoptions.TERMIMPROVEGENS
            mpgastate.TermOpt = 1;
        end

        % migrate individuals between subpopulations
        if ((mpgastate.TermOpt ~= 1) && (rem(mpgastate.Generation, mpgaoptions.MIGGEN) == 0))
            [mpgastate.Chrom, mpgastate.ObjV] = migrate(mpgastate.Chrom, mpgaoptions.SUBPOP, [mpgaoptions.MIGR, 1, 0], mpgastate.ObjV);
        end

    end

    % save the final state of the algorithm
    if mpgaoptions.SAVEMODE == 1
        try
            save(mpgaoptions.FILENAME, 'mpgaoptions', 'mpgastate');
        catch
            warning('MPGA:savefinalstatefailed', ...
                'Attempt to save final mpga state file failed');
        end
    end
    
    % Results
    
    % add number of objective function evaluations
    Results = cumsum(mpgastate.ScoreTracking(1:mpgastate.Generation,3));
    
    % number of function evaluation, best and mean results
    Results = [Results mpgastate.ScoreTracking(1:mpgastate.Generation,1) mpgastate.ScoreTracking(1:mpgastate.Generation,2)];

    if mpgaoptions.DISPLAYMODE == 1 || mpgaoptions.DISPLAYMODE == 2

        % Plot final results
        finalresplot(Results, mpgastate.BestChroms, mpgastate.FigTitle, mpgastate.Generation, mpgaoptions.STEP);

    end

    if mpgaoptions.DISPLAYMODE == 1 || mpgaoptions.DISPLAYMODE == 2 || mpgaoptions.DISPLAYMODE == 3

        % Display the best individual and score
        [minObjV, ix] = min(mpgastate.ObjV);

        fprintf(['\nTotal No. of Function Evaluations:  ', num2str(Results(mpgastate.Generation,1)),'\n\n']);

        fprintf(['mpgastate.ScoreTracking score:  ', num2str(minObjV, 8),'\n\nMost Successful Individual:\n\n']);

        for k = 1:length(mpgastate.Chrom(ix,:))

            disp(['Var. ', num2str(k), ':   ',num2str(mpgastate.Chrom(ix,k), 6)])

        end

        fprintf('\n\n');

    end
    
    

end
% End of function



function hfig = iteroutput(mpgastate, mpgaoptions, hfig)
% display GA ouput at current iteration depending on user preferences

    

    % Plot some results, rename title of figure for graphic output
    if (mpgaoptions.DISPLAYMODE == 1) ...
            && ( (rem(mpgastate.Generation, mpgaoptions.STEP) == 1) ...
                    || (rem(mpgastate.Generation, mpgaoptions.MAXGEN) == 0) ...
                    || (mpgastate.TermOpt == 1) ...
               )

        if isempty(hfig)
            hfig = figure;
        end

        set(hfig,'Name',[mpgastate.FigTitle ' in ' int2str(mpgastate.Generation)], 'Renderer', 'painters');

        resplot( mpgastate.Chrom(1:2:size(mpgastate.Chrom,1),:),... % Selects every second value in mpgastate.Chrom
                 mpgastate.BestChroms(max(1,mpgastate.Generation-39):size(mpgastate.BestChroms,1),:),... % Selects the last 40 generations' best individuals variable values
                 mpgastate.ObjV,... % The current generations individuals' scores
                 mpgastate.ScoreTracking(max(1,mpgastate.Generation-19):mpgastate.Generation,[1 2]),... % Selects the last 20 generations' best individuals scores and mean scores
                 mpgastate.Generation, ... % The current generation number
                 hfig ); 

    elseif (mpgaoptions.DISPLAYMODE == 2 || mpgaoptions.DISPLAYMODE == 3) ...
            && ( mpgaoptions.STEP == 1 ...
                  || (rem(mpgastate.Generation, mpgaoptions.STEP) == 1) ...
                  || (rem(mpgastate.Generation, mpgaoptions.MAXGEN) == 0) ...
                  || (mpgastate.TermOpt == 1) ...
               )

        [minObjV, ix] = min(mpgastate.ObjV);

        fprintf(['\nGen: ', int2str(mpgastate.Generation),' best score:  ', num2str(minObjV, 8),'\n\nMost Successful Individual:\n\n']);

        for k = 1:length(mpgastate.Chrom(ix,:))

            disp(['Var. ', num2str(k), ':   ',num2str(mpgastate.Chrom(ix,k), 6)])

        end

        fprintf('\n');

    end
        
end



function options = parsempgaoptions(options, objargs)
% checks and fills the options structure necessary for mpga.m with default
% values where missing
%
% Input:
%
%   options - a structure containing some of the following members:
%
%       OBJ_F - string or function handle to the objective function
%
%       XOV_F - string or function handle to the recombination function for 
%               individuals
%
%       MUT_F - string or function handle to the mutation function
%
%       SEL_F - string or function handle to the selection function
%
%       GGAP - (optional scalar) Generation gap, how many new individuals
%              are created, default is 0.8
%
%       INSR - (optionsl scalar) Insertion rate, how many of the offspring
%              are inserted, default is 0.9
%
%       XOVR - (optional scalar) Crossover rate, default is 1
%
%       SP - (optional scalar) Selective Pressure, default is 2
%
%       MUTR - (optional scalar) Mutation rate, default value is calculated
%              based on the number of variables in the objective function
%
%       MIGR - (optional scalar) Migration rate between subpopulations,
%              default is 0.1
%
%       MAXGEN - (optional scalar) Max number of generations, default is
%                calculated from the number of variables
%
%       MIGGEN - (optional scalar) Number of generations between migration
%                (isolation time), default is calculated as 
%                ceil(MAXGEN / 10)
%                
%       TERMEXACT - (optional scalar) Value for termination if minimum
%                   reached, default is 1e-8
%
%       SUBPOP - (optional scalar) Number of subpopulations, default is
%                calculated based on the number of variables
%
%       NIND - (optional scalar) Number of individuals per subpopulations,
%              default is calculated from the number of variables
%
%       MUTR - (scalar) Mutation rate depending on NVAR
%
%       STEP - (optional scalar) The variables of the best individuals will
%              be plotted for each generation sparated by an interval of
%              size STEP. Default is set to ceil(MAXGEN / 10)
%
%       DISPLAYMODE - (optional scalar) If 1, some graphs are displayed at
%                     each interval in STEP, if 2, some text output is
%                     displayed at each interval, but no graphs are
%                     plotted. Otherwise, no output is displayed
%
%       SAVEMODE - (optional scalar) If 0, no data is saved, otherwise, if no
%                  filename is supplied data is perioically saved in a
%                  default file name: "OBJ_F" + _output.mat which can be
%                  used to resume the algorithm from the current
%                  generation. If a filename is supplied, the data is saved
%                  in that file.
%
%       FILENAME - optional name of a file in which to periodically save data
%                  about the progress of the  GA, if not supplied, and SAVEMODE
%                  is not zero, a default file name will be used. The
%                  default filename will be the name of the objective
%                  function with the string '_output.mat' appended.
%
%       RESUMEFILE - name of file where data from a previous run can be found.
%                    If supplied, the GA will resume from this previous
%                    starting point
%
%       EVALUATE - true or false flag, the objective function will be
%                  evaluated and returned with all output arguments if true
%
% Output:
%
%   options - fully filled options structure using the default values where
%             they have been omitted by the user
%

    % handle case where user just supplied the resume file name
    if ischar (options) || (isstring (options) && isscalar (options))
        
        ga_save_file = options;

        if exist (ga_save_file, 'file') == 2
            
            options = struct ();
            
            options.EVALUATE = true;
            options.RESUMEFILE = ga_save_file;

        else
            error ('The supplied GA state file: %s does not appear to exist.', ga_save_file);
        end
    elseif ~isstruct (options)
        error ('mpgaoptions must be a structure or a character vector or scalar string containing the path of a GA state file.');
    end
        
    % convert all field names to upper case
    options = upperfnames(options);

    if ~isfield (options, 'EVALINDS')
        options.EVALINDS = [];
    end
    
    if ~isfield (options, 'EVALUATE')
        options.EVALUATE = false;
    end
    
    if options.EVALUATE
        return;
    end
    
    % get objective function
    if ~isfield(options, 'OBJ_F')
        error('Must supply an objective function');
    else
        if ischar(options.OBJ_F)
            options.OBJ_F = str2func(options.OBJ_F);
        end
    end

    % get crossover function
    if ~isfield(options, 'XOV_F')
        error('You must supply a crossover function')
    else
        if ischar(options.XOV_F)
            options.XOV_F = str2func(options.XOV_F);
        end
    end

    % get the mutation function
    if ~isfield(options, 'MUT_F')
        error('You must supply a mutation function')
    else
        if ischar(options.MUT_F)
            options.MUT_F = str2func(options.MUT_F);
        end
    end

    % Get the selection function
    if ~isfield(options, 'SEL_F')
        % default to stochastic universal selection
        options.SEL_F = @sus;
    else
        if ischar(options.SEL_F)
            options.SEL_F = str2func(options.SEL_F);
        end
    end

    % default generation gap is 0.8
    if ~isfield(options, 'GGAP')
        options.GGAP = 0.8;
    else
        if options.GGAP < 0 || options.GGAP > 1
            error('Invalid generation gap value, must be between 0 and 1')
        end
    end

    % default insertion rate is 0.9
    if ~isfield(options, 'INSR')
        options.INSR = 0.9;
    else
        if options.INSR <= 0 || options.INSR > 1
            error('Invalid insertion rate, must be greater than zero and less than 1')
        end
    end

    % default crossover rate if 1
    if ~isfield(options, 'XOVR')
        options.XOVR = 1;
    end

    if ~isfield(options, 'SP')
        options.SP = 2;
    end
    
    % Get number of variables from objective function
    tempArgs = [{[],1}, objargs];
    FieldDR = feval(options.OBJ_F, tempArgs{:});
    NVAR = size(FieldDR,2);

    if ~isfield(options, 'NIND')
        % Number of individuals per subpopulations
        options.NIND = 20 + 5 * floor(NVAR/50);
    end

    if ~isfield(options, 'MAXGEN')
        % Max number of generations
        options.MAXGEN = 300 * floor(sqrt(NVAR));
    end

    if ~isfield(options, 'SUBPOP')
        % Number of subpopulations
        options.SUBPOP = 2 * floor(sqrt(NVAR));
    else
        if options.SUBPOP < 1
            error('You must specify at least 1 subpopulation')
        end
    end

    % Default migration rate is 0.1
    if ~isfield(options, 'MIGR')
        options.MIGR = 0.1;
    else
        if options.MIGR <= 0 || options.MIGR > 1
            error('Invalid migration rate, must be greater than zero and less than 1')
        end
    end

    % Calculate mutation rate based on the number of variables
    if ~isfield(options, 'MUTR')
        % Mutation rate depending on NVAR
        options.MUTR = 1 / NVAR;
    end

    % number of generations between migration events
    if ~isfield(options, 'MIGGEN')
        options.MIGGEN = ceil(options.MAXGEN / 10);
    end

    % default terminal value is 1e-8
    if ~isfield(options, 'TERMEXACT')
        options.TERMEXACT = 1e-8;
    end
    
    % default terminal value is 1e-8
    if ~isfield(options, 'TERMIMPROVETOL')
        options.TERMIMPROVETOL = 1e-8;
    end
    
    if ~isfield(options, 'TERMIMPROVEGENS')
        options.TERMIMPROVEGENS = inf;
    end

    if ~isfield(options, 'DISPLAYMODE')
        options.DISPLAYMODE = 0;
        % Determines how frequently graphs of progress will be displayed
        if ~isfield(options, 'STEP')
            options.STEP = options.MAXGEN;
        end
    else
        if options.DISPLAYMODE ~= 0 && options.DISPLAYMODE ~= 1 && options.DISPLAYMODE ~= 2
            
            options.DISPLAYMODE = 0;
            % Determines how frequently graphs of progress will be displayed
            if ~isfield(options, 'STEP')
                options.STEP = options.MAXGEN;
            end
            warning('DISPLAYMODE was an invalid value, therefore this has been set to 0')
            
        elseif options.DISPLAYMODE == 1
            
            if ~isfield(options, 'STEP')
                options.STEP = ceil(options.MAXGEN / 10);
            end      
            
        elseif options.DISPLAYMODE == 2
            
            if ~isfield(options, 'STEP')
                options.STEP = 1;
            end
            
        end
    end

    if ~isfield(options, 'SAVEMODE')
        options.SAVEMODE = 0;
    else
        if options.SAVEMODE ~= 0 && options.SAVEMODE ~= 1
            options.SAVEMODE = 0;
            warning('SAVEMODE was an invalid value, therefore this has been set to 0')
        elseif options.SAVEMODE == 1
            if ~isfield(options, 'FILENAME')
                options.FILENAME = [func2str(options.OBJ_F), '_output.mat'];
            end
        end
    end

    if ~isfield(options, 'RESUMEFILE')
        options.RESUME = 0;
    else
        options.RESUME = 1;
    end
    
    if ~isfield(options, 'OUTPUTFCN')
        options.OUTPUTFCN = [];
    end
    
    if ~isfield(options, 'OUTPUTFCNARGS')
        options.OUTPUTFCNARGS = {};
    end

end
