echo sim
root  -l -q sim_complete.C"(1000,4.6)" >> sim.log 2>&1 
echo digi
root  -l -q digi_complete.C >> digi.log 2>&1 
echo reco
root  -l -q recoideal_complete.C >> reco.log 2>&1
echo pid
root  -l -q pid_complete.C >> pid.log 2>&1 
echo END

