echo sim
root  -l -q sim_complete.C"(1000)">> sim.log 2>&1 
echo digi
root  -l -q digi_complete_boxgen.C >> digi.log 2>&1 
echo reco
root  -l -q reco_complete_boxgen.C >> reco.log 2>&1
echo pid
root  -l -q pid_complete_boxgen.C >> pid.log 2>&1 
#echo ana
#root  -l -q ana_complete.C >> ana.log 2>&1
echo END

