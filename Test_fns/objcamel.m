% OBJCAMEL.M      (OBJective function for Camel function)
%
% This function implements the HARVEST PROBLEM.
%
% Syntax:  ObjVal = objharv(Chrom,rtn_type)
%
% Input parameters:
%    Chrom     - Matrix containing the chromosomes of the current
%                population. Each row corresponds to one individual's
%                string representation.
%                if Chrom == [], then special values will be returned
%    rtn_type  - if Chrom == [] and
%                rtn_type == 1 (or []) return boundaries
%                rtn_type == 2 return title
%                rtn_type == 3 return value of global minimum
%
% Output parameters:
%    ObjVal    - Column vector containing the objective values of the
%                individuals in the current population.
%                if called with Chrom == [], then ObjVal contains
%                rtn_type == 1, matrix with the boundaries of the function
%                rtn_type == 2, text for the title of the graphic output
%                rtn_type == 3, value of global minimum
%                
%
% Author:     Hartmut Pohlheim
% History:    18.02.94     file created (copy of vallinq.m)
%             01.03.94     name changed in obj*
%             14.01.03     updated for MATLAB v6 by Alex Shenfield

function ObjVal = objcamel(Chrom,rtn_type);

% global gen;

% Dimension of objective function
   Dim = 2;
   
% Compute population parameters
   [Nind,Nvar] = size(Chrom);

% Check size of Chrom and do the appropriate thing
   % if Chrom is [], then define size of boundary-matrix and values
   if Nind == 0
      % return text of title for graphic output
      if rtn_type == 2
         ObjVal = ['CAMEL FUNCTION-' int2str(Dim)];
      % return value of global minimum
      elseif rtn_type == 3
         ObjVal = -2;
      % define size of boundary-matrix and values
      else   
         % lower and upper bound, identical for all n variables        
         ObjVal1 = [-2; 2];
         ObjVal = rep(ObjVal1,[1 Dim]);
      end
   % if Dim variables, compute values of function
   elseif Nvar == Dim
       
      ObjVal = CamelFun(Chrom(:,1),Chrom(:,2));
      
   % otherwise error, wrong format of Chrom
   else
      error('size of matrix Chrom is not correct for function evaluation');
   end   


% End of function