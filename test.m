plot(1:7, arrayfun(@(N) distSpanned(.5, .02, N), 1:7), 'Color', 'red')
hold
dists = 0:.01:3.5
plot(arrayfun(@(L) numAgents(.5, .02, L), dists), dists, 'Color', 'green')
