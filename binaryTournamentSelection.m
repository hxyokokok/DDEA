function idx = binaryTournamentSelection(expectation,N)

% binary tournament
tournamentGroup = randi(size(expectation,1),N,2);
tournamentResult = expectation(tournamentGroup(:,1)) > expectation(tournamentGroup(:,2));
idx = tournamentGroup(:,1) .* tournamentResult + tournamentGroup(:,2) .* ~tournamentResult;