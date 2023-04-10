function [mpgaoptions, mpgastate] = expandpop(method, n, mpgaoptions, mpgastate)
% expands the population of a ga
%
% Syntax
%
% [mpgaoptions, mpgastate] = expandpop(method, n, mpgaoptions, mpgastate)
%
% Input
%
%  method - A string containing the method by which the population may be
%    expanded. This can be either 'subpops' or 'nind'. 
%
%    If 'subpops', entire subpopulations of new individuals are added to
%    the population, each containing the same number of individuals as are
%    present in the existing subpopulations (i.e. the value in
%    mpgaoptions.NIND). The number of the populations added is the value
%    supplied in 'n'.
%
%    If 'nind' new individuals are added to every subpopulation already
%    present. 'n' new individuals are added to every subpopulation.
% 
%  n - Either the number of subpopulations to be added, or the number of
%   individuals to be added to each subpopulation depending on the value of
%   'method'.
% 
%  mpgaoptions - The existing mpgaoptions structure for the ga.
% 
%  mpgastate - The existing mpgastate structure for the ga.
%
%

    if strcmpi(method, 'subpops')
        
        % Create new individuals
        newChrom = crtrp(n * mpgaoptions.NIND, mpgastate.FieldDR);
        
        % evaluate the new individuals
        objargs = [{newChrom, []}, mpgastate.ObjectiveArgs];
        newObjV = feval(mpgaoptions.OBJ_F, objargs{:});
        
        mpgastate.Chrom = [ mpgastate.Chrom; newChrom ];
        mpgastate.ObjV = [ mpgastate.ObjV; newObjV ];
        
        mpgaoptions.SUBPOP = mpgaoptions.SUBPOP + n;
        
    elseif strcmpi(method, 'nind')
        
        % Create new individuals
        newChrom = crtrp(mpgaoptions.SUBPOP * n, mpgastate.FieldDR);
        
        % evaluate the new individuals
        objargs = [{newChrom, []}, mpgastate.ObjectiveArgs];
        newObjV = feval(mpgaoptions.OBJ_F, objargs{:});
        
        newpopChrom = [];
        newpopObjV = [];
        
        % perform insertion for each subpopulation
        for ind = 0:(mpgaoptions.SUBPOP-1)
            % Calculate positions in old subpopulation, where offspring are
            % inserted
            newpopChrom = [ newpopChrom; 
                            mpgastate.Chrom(((ind*mpgaoptions.NIND)+1):((ind+1)*mpgaoptions.NIND));
                            newChrom(((ind*n)+1):((ind+1)*n));
                           ];
                       
            newpopObjV = [ newpopObjV;
                           mpgastate.ObjV(((ind*mpgaoptions.NIND)+1):((ind+1)*mpgaoptions.NIND));
                           newObjV(((ind*n)+1):((ind+1)*n));
                         ];
            
        end
        
        mpgastate.Chrom = newpopChrom;
        mpgastate.ObjV = newpopObjV;
        
        % note the new subpopulation sizes in the mpgaoptions structure
        mpgaoptions.NIND = mpgaoptions.NIND + n;
        
    else
        
        error('Unrecognised expansion method.');
        
    end

end