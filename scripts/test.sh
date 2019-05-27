p=../dist/Debug/GNU-Linux/sequence_alignment100
csp=--memo_caching_strategy
klp=--memo_key_length
vlp=--memo_value_length
saifp=--sequence_align_instance_fname



cs=non_memo
#saif=../../misc_phd/input/lcs_instances/10-random-0-15-1357.lcs
#saif=1-test-0-0-0.lcs
#saif=5-test-0-4-1.lcs
saif=5-test-0-4-1.lcs
kl=8
vl=8

$p $csp $cs $klp $kl $vlp $vl $saifp $saif
saif=5-test-0-4-2.lcs
$p $csp $cs $klp $kl $vlp $vl $saifp $saif
saif=5-test-0-4-3.lcs
$p $csp $cs $klp $kl $vlp $vl $saifp $saif
saif=5-test-0-4-4.lcs
$p $csp $cs $klp $kl $vlp $vl $saifp $saif
