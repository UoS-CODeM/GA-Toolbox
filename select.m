function SelCh = select2(SEL_F, Chrom, FitnV, GGAP, SUBPOP, SEL_F_ARGS)
% SELECT.M          (universal SELECTion)
%
% This function performs universal selection. The function handles
% multiple populations and calls the low level selection function
% for the actual selection process.
%
% Syntax:  SelCh = select(SEL_F, Chrom, FitnV, GGAP, SUBPOP)
%
% Input parameters:
%    SEL_F     - Name of the selection function
%    Chrom     - Matrix containing the individuals (parents) of the current
%                population. Each row corresponds to one individual.
%    FitnV     - Column vector containing the fitness values of the
%                individuals in the population.
%    GGAP      - (optional) Rate of individuals to be selected
%                if omitted 1.0 is assumed
%    SUBPOP    - (optional) Number of subpopulations
%                if omitted 1 subpopulation is assumed
%    SEL_F_ARGS - (optional) Additional arguments to be passed to the
%                 selection function in the form of a cell array
% Output parameters:
%    SelCh     - Matrix containing the selected individuals.
%
% Author:     Hartmut Pohlheim
% History:    10.03.94     file created
%             22.01.03     tested under MATLAB v6 by Alex Shenfield
%             08.03.10     modified by Richard Crozier to allow the user to
%                          specify additional input arguments to the selection 
%                          function

% Check parameter consistency
   if nargin < 3, error('Not enough input parameter'); end

   % Identify the population size (Nind)
   [NindCh,Nvar] = size(Chrom);
   [NindF,VarF] = size(FitnV);
   if NindCh ~= NindF, error('Chrom and FitnV disagree'); end
   if VarF ~= 1, error('FitnV must be a column vector'); end
  
   if nargin < 5, SUBPOP = 1; end
   if nargin > 4,
      if isempty(SUBPOP), SUBPOP = 1;
      elseif isnan(SUBPOP), SUBPOP = 1;
      elseif length(SUBPOP) ~= 1, error('SUBPOP must be a scalar'); end
   end

   if (NindCh/SUBPOP) ~= fix(NindCh/SUBPOP), error('Chrom and SUBPOP disagree'); end
   Nind = NindCh/SUBPOP;  % Compute number of individuals per subpopulation

   if nargin < 4, GGAP = 1; end
   if nargin > 3,
      if isempty(GGAP), GGAP = 1;
      elseif isnan(GGAP), GGAP = 1;
      elseif length(GGAP) ~= 1, error('GGAP must be a scalar');
      elseif (GGAP < 0), error('GGAP must be a scalar bigger than 0'); end
   end

% Compute number of new individuals (to select)
   NSel=max(floor(Nind*GGAP+.5),2);

   if nargin < 6, SEL_F_ARGS = {}; end
% Initialize the selection function arguments
   tempSelArgs = [{[], NSel}, SEL_F_ARGS];
   
% Select individuals from population
   SelCh = [];
   for irun = 1:SUBPOP,
      FitnVSub = FitnV((irun-1) * Nind+1:irun*Nind);
      tempSelArgs{1} = FitnVSub;
      ChrIx = feval(SEL_F, tempSelArgs{:}) + (irun-1)*Nind;
      %ChrIx = feval(SEL_F, FitnVSub, NSel) + (irun-1)*Nind;
      SelCh = [SelCh; Chrom(ChrIx,:)];
   end
 
end
% End of function