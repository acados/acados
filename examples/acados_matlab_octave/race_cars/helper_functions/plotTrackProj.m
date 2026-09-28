%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%
% Author: Daniel Kloeser
% Ported by Thomas Jespersen (thomasj@tkjelectronics.dk), TKJ Electronics
%

function plotTrackProj(simX, filename)
    % load track
    s=simX(:,1);
    n=simX(:,2);
    alpha=simX(:,3);
    v=simX(:,4);
    distance=0.12;
    % transform data
    [x, y, ~, ~] = transformProj2Orig(s, n, alpha, v,filename);
    % plot racetrack map

    % Setup plot
    ylim([-1.75, 0.35]);
    xlim([-1.1, 1.6]);    
    ylabel('y[m]');
    xlabel('x[m]');

    % Plot center line
    [Sref,Xref,Yref,Psiref,~]=getTrack(filename);
    plot(Xref,Yref,'k--');

    % Draw Trackboundaries
    Xboundleft=Xref-distance*sin(Psiref);
    Yboundleft=Yref+distance*cos(Psiref);
    Xboundright=Xref+distance*sin(Psiref);
    Yboundright=Yref-distance*cos(Psiref);
    hold on;
    plot(Xboundleft,Yboundleft,'k','LineWidth',1);
    plot(Xboundright,Yboundright,'k','LineWidth',1);
    plot(x,y, '-b');

%   Not ported yet:
%     % Draw driven trajectory
%     heatmap = plt.scatter(x,y, c=v, cmap=cm.rainbow, edgecolor='none', marker='o')
%     cbar = plt.colorbar(heatmap, fraction=0.035)
%     cbar.set_label("velocity in [m/s]")
%     ax = plt.gca()
%     ax.set_aspect('equal', 'box')
% 
%     % Put markers for s values
%     xi=np.zeros(9)
%     yi=np.zeros(9)
%     xi1=np.zeros(9)
%     yi1=np.zeros(9)
%     xi2=np.zeros(9)
%     yi2=np.zeros(9)
%     for i in range(int(Sref[-1]) + 1):
%         try
%             k = list(Sref).index(i + min(abs(Sref - i)))
%         catch
%             k = list(Sref).index(i - min(abs(Sref - i)))
%         end
%         [_,nrefi,_,_]=transformOrig2Proj(Xref[k],Yref[k],Psiref[k],0)
%         [xi[i],yi[i],_,_]=transformProj2Orig(Sref[k],nrefi+0.24,0,0)
%         % plt.text(xi[i], yi[i], f'{i}m', fontsize=12,horizontalalignment='center',verticalalignment='center')
%         plt.text(xi[i], yi[i], '{}m'.format(i), fontsize=12,horizontalalignment='center',verticalalignment='center')
%         [xi1[i],yi1[i],_,_]=transformProj2Orig(Sref[k],nrefi+0.12,0,0)
%         [xi2[i],yi2[i],_,_]=transformProj2Orig(Sref[k],nrefi+0.15,0,0)
%         plt.plot([xi1[i],xi2[i]],[yi1[i],yi2[i]],color='black')
%     end
    hold off;
end