echo sim
root  -l -q sim_complete.C"(10000,4.6)" >> sim.log 2>&1 
echo digi
root  -l -q digi_complete.C >> digi.log 2>&1 
echo reco
root  -l -q recoideal_complete.C >> reco.log 2>&1
echo pid
root  -l -q pid_complete.C >> pid.log 2>&1 
#echo ana
#root  -l -q /home/ikp1/puetz/panda/myscripts/analysis/AnalysisTaskXi1820.C"(4.6)" >> ana.log 2>&1
echo END

