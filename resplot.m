function resplot(Chrom,IndAll,ObjV,Best,gen,hfig)
% RESPLOT.M       (RESult PLOTing)
%
% This function plots some results during computation.
%
% Syntax:  resplot(Chrom,IndAll,ObjV,Best,gen)
%
% Input parameters:
%    Chrom     - Matrix containing the chromosomes of the current
%                population. Each line corresponds to one individual.
%    IndAll    - Matrix containing the best individual (variables) of each
%                generation. Each line corresponds to one individual.
%    ObjV      - Vector containing objective values of the current
%                generation
%    Best      - Matrix containing the best and average Objective values of 
%                each generation, [best value per generation,average value 
%                per generation]
%    gen       - Scalar containing the number of the current generation
%
% Output parameter:
%    no output parameter
%
%  Author:    Hartmut Pohlheim
%  History:   27.11.93     file created
%             29.11.93     decision, if plot or not deleted
%                          yscale not log
%             15.12.93     MutMatrix as parameter and plot added
%             16.03.94     function cleaned, MutMatrix removed, IndAll added
%             22.01.03     tested under MATLAB v6 by Alex Shenfield

  set(0,'CurrentFigure',hfig);

   % plot of best and mean value of passed in generations
      subplot(2,2,1), plot(1:size(Best,1), Best(:,1)', '-r', 1:size(Best,1), Best(:,2)',  '--b');
      set(gca, 'xlim', [0 (size(Best,1) + 1)])
      title(sprintf(['Best and Mean Objective Values\n(last ', num2str(size(Best,1)), ' Generations)']));
      xlabel('Generation'), ylabel('Objective Values');

   % plot of best individuals variable values in passed in generations
   
      subplot(2,2,2), bar(IndAll, 'group');
      set(gca, 'xlim', [0 (size(IndAll,1) + 1)])
      set(gca, 'XTickLabel', (gen-size(IndAll,1)):5:(gen+1))
      title(sprintf(['Best Individuals Variable Values\n(last ', num2str(size(IndAll,1)), ' Generations)']));
      xlabel('Generation'), ylabel('Value of Variables');

   % plot of variables of selected individuals in current generation
      subplot(2,2,3), bar(Chrom, 'group'); % plot(Chrom');
      set(gca, 'xlim', [0 (size(Chrom,1) + 1)])
      title(['Variables of Selected Individuals Gen: ',num2str(gen)]);
      xlabel('Individual Index'), ylabel('Variable Values');
    
   % plot of all objective values in current generation
      subplot(2,2,4), plot(ObjV,'*b');
      set(gca, 'xlim', [0 (length(ObjV) + 1)])
      title(['All Objective Values Gen: ',num2str(gen)]);
      xlabel('Current Generation'), ylabel('Objective Values');

      drawnow;

end
% End of function