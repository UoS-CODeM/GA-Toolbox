function finalresplot(Results, IndAll, FigTitle, gen, step)
% Plot GA Results and show best individuals => optimum

    %scrsz = get(0,'ScreenSize');
    %figure('Name',['Results of ' FigTitle], 'Position',[1 scrsz(4) scrsz(3) scrsz(4)], 'Renderer', 'painters')
    
    figure('Name',['Results of ' FigTitle], 'Renderer', 'painters');

    % First plot best scores per generation
    subplot(2,2,1), plot(1:gen,Results(:,2),'-r');
    title('Best Scores Versus Generation');
    xlabel('Generation'), ylabel('Best Score');

    % Next plot mean scores per generation
    subplot(2,2,2), plot(1:gen,Results(:,3),'--b');
    title('Mean Scores Versus Generation');
    xlabel('Generation'), ylabel('Mean Score');

    % Plot best individual's variables at every third generation
    subplot(2,2,3:4), bar(IndAll(1:step:gen,:), 'group');
    set(gca, 'xlim', [0 (size(IndAll(1:step:gen,:),1)+1)])
    title(sprintf(['Best Individuals Variable Values (last ', num2str(size(IndAll,1)), ' Generations)']));
    xlabel(['Generation / ',num2str(step)]), ylabel('Value of Variables');
    
    drawnow;

end