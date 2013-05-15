PRO specfit_libpars,lib_mu,lib_tau_v,lib_alpha_bc,lib_alpha,vers_sps=vers_sps

if vers_sps EQ 'cb07' then begin
   lib_mu=0.3
   lib_tau_v=1.0
   lib_alpha_bc=1.3
   lib_alpha=0.7
endif
if vers_sps EQ 'bc03' then begin
   lib_mu=0.3
   lib_tau_v=1.0
   lib_alpha_bc=1.3
   lib_alpha=0.7
endif
if vers_sps EQ 'fsps' then begin
   lib_mu=0.3
   lib_tau_v=1.0
   lib_alpha_bc=0.7
   lib_alpha=0.7
endif

end
