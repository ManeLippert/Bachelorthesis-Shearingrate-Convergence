read_perf('nonspec_scaling2b','pa_*','pa_64_mu16sp2s2',1)
read_perf('nonspec_scaling2b','ca_*','ca_64_mu8sp2s4',1)

%legend('Med ES, No Coll','linear','Med ES Coll','linear')

read_perf('nonspec_scaling3','s*','s1024_mu16sp2s32',1)
read_perf('nonspec_scaling3_clean_coll','ac*','ac_2048_mu8sp2s32x4',1)
read_perf('nonspec_scaling3_clean','a*','a2_1024_mu16sp2s32',1)

%legend('Large ES, No Coll','linear','Large EM, Coll','linear','Large EM, No Coll','linear')
