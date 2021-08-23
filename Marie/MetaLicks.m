[RTs, JuiceLicks] =findRT(JuiceLicks);
AllLicks = MakeAllLicks(JuiceLicks);
NoJuice = ExtractNoJuice(JuiceLicks);
FirstLicksEpochs = FindLickingEpochs(AllLicks);
FirstJuice = ExtractFirstJuice(JuiceLicks);
save FirstLicksEpochs.txt FirstLicksEpochs -ascii -double
save AllLicks.txt AllLicks -ascii -double
save RTs.txt RTs -ascii -double
save NoJuice.txt NoJuice -ascii -double
save FirstJuice.txt FirstJuice -ascii -double