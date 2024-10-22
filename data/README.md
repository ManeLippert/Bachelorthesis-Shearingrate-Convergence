# Data Structure

This folder contains the data for my Bachelor Thesis which were simulated with gkw with adiabtic electrons approximation of the type S6 nearby the finite heatflux threshold.

## Data 

Markdown version of [data.csv](./data.csv)

| **case** | **rlt** | **boxsize** | **time** | **timestep** | **Ns** | **Nvpar** | **Nmu** | **dtim** | **krhomax** | **Nmod** | **Nx** | **ikx** | **finit** | **stepsize_all** | **evo_start**                           | **evo_end**                              | **sel_start**    | **sel_end**      | **error_index** | **converge** | **mode** | **backup** | **path**                                                        |
|----------|---------|-------------|----------|--------------|--------|-----------|---------|----------|-------------|----------|--------|---------|-----------|------------------|-----------------------------------------|------------------------------------------|------------------|------------------|-----------------|--------------|----------|------------|-----------------------------------------------------------------|
| S6       | 6.0     | 1x1         | 6000.0   | 10000.0      | 12     | 16        | 6       | 0.02     | 1.4         | 21       | 83     | 5       | cosine2   | 1000             | "[500, 3000, 5000]"                     | "[2000, 4000, 5500]"                     |                  |                  |                 | False        | 0        | True       | /data/S6_rlt6.0/boxsize1x1/Ns12/Nvpar16/Nmu6                    |
| S6       | 6.0     | 1x1         | 12000.0  | 20000.0      | 12     | 32        | 6       | 0.02     | 1.4         | 21       | 83     | 5       | cosine2   | 1000             | "[500, 3000, 5000, 6000,  8000, 11000]" | "[2000, 4000, 5500, 7000, 10000, 12000]" |                  |                  |                 | False        | 0        | True       | /data/S6_rlt6.0/boxsize1x1/Ns12/Nvpar32/Nmu6                    |
| S6       | 6.0     | 1x1         | 6000.0   | 10000.0      | 12     | 48        | 9       | 0.02     | 1.4         | 21       | 83     | 5       | cosine2   | 1000             | "[500, 3000, 5000]"                     | "[2000, 4000, 5500]"                     |                  |                  |                 | False        | 0        | True       | /data/S6_rlt6.0/boxsize1x1/Ns12/Nvpar48/Nmu9                    |
| S6       | 6.0     | 1x1         | 6000.0   | 10000.0      | 12     | 64        | 9       | 0.02     | 1.4         | 21       | 83     | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | False        | 0        | True       | /data/S6_rlt6.0/boxsize1x1/Ns12/Nvpar64/Nmu9                    |
| S6       | 6.0     | 1x1         | 12000.0  | 20000.0      | 16     | 16        | 9       | 0.02     | 1.4         | 21       | 83     | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | False        | 0        | True       | /data/S6_rlt6.0/boxsize1x1/Ns16/Nvpar16/Nmu9                    |
| S6       | 6.0     | 1x1         | 12000.0  | 20000.0      | 16     | 32        | 9       | 0.02     | 1.4         | 21       | 83     | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | False        | 0        | True       | /data/S6_rlt6.0/boxsize1x1/Ns16/Nvpar32/Nmu9                    |
| S6       | 6.0     | 1x1         | 6000.0   | 10000.0      | 16     | 48        | 6       | 0.02     | 1.4         | 21       | 83     | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | False        | 0        | True       | /data/S6_rlt6.0/boxsize1x1/Ns16/Nvpar48/Nmu6                    |
| S6       | 6.0     | 1x1         | 6000.0   | 10000.0      | 16     | 48        | 9       | 0.02     | 1.4         | 21       | 83     | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | True         | 1        | True       | /data/S6_rlt6.0/boxsize1x1/Ns16/Nvpar48/Nmu9                    |
| S6       | 6.0     | 1x1         | 7125.0   | 10000.0      | 16     | 48        | 9       | 0.025    | 1.4         | 21       | 83     | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | True         | 1        | True       | /data/S6_rlt6.0/boxsize1x1/Ns16/Nvpar48/Nmu9/dtim0.025          |
| S6       | 6.0     | 1x1         | 6000.0   | 10000.0      | 16     | 48        | 9       | 0.02     | 0.7         | 11       | 83     | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | True         | 1        | True       | /data/S6_rlt6.0/boxsize1x1/Ns16/Nvpar48/Nmu9/krhomax0.70/Nmod11 |
| S6       | 6.0     | 1x1         | 6000.0   | 10000.0      | 16     | 48        | 9       | 0.02     | 1.4         | 21       | 43     | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | False        | 0        | True       | /data/S6_rlt6.0/boxsize1x1/Ns16/Nvpar48/Nmu9/Nx43               |
| S6       | 6.0     | 1x1         | 6000.0   | 10000.0      | 16     | 48        | 9       | 0.02     | 1.4         | 21       | 63     | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | False        | 0        | True       | /data/S6_rlt6.0/boxsize1x1/Ns16/Nvpar48/Nmu9/Nx63               |
| S6       | 6.0     | 1x1         | 12000.0  | 20000.0      | 16     | 64        | 6       | 0.02     | 1.4         | 21       | 83     | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | True         | 1        | True       | /data/S6_rlt6.0/boxsize1x1/Ns16/Nvpar64/Nmu6                    |
| S6       | 6.0     | 1x1         | 12000.0  | 20000.0      | 16     | 64        | 9       | 0.02     | 1.4         | 21       | 83     | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | True         | 1        | True       | /data/S6_rlt6.0/boxsize1x1/Ns16/Nvpar64/Nmu9                    |
| S6       | 6.0     | 1.5x1.5     | 6000.0   | 10000.0      | 16     | 48        | 9       | 0.02     | 1.4         | 31       | 125    | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     | [7000]           | [8000]           |                 | True         | 2        | False      | /data/S6_rlt6.0/boxsize1.5x1.5/Ns16/Nvpar48/Nmu9                |
| S6       | 6.0     | 2x1         | 18000.0  | 30000.0      | 16     | 48        | 9       | 0.02     | 1.4         | 21       | 165    | 10      | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | True         | 2        | True       | /data/S6_rlt6.0/boxsize2x1/Ns16/Nvpar48/Nmu9                    |
| S6       | 6.0     | 2x2         | 18000.0  | 30000.0      | 16     | 48        | 9       | 0.02     | 1.4         | 41       | 165    | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | True         | 2        | True       | /data/S6_rlt6.0/boxsize2x2/Ns16/Nvpar48/Nmu9                    |
| S6       | 6.0     | 2.5x2.5     | 3402.1   | 5672         | 16     | 48        | 9       | 0.02     | 1.4         | 51       | 207    | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     | [2000]           | [3000]           |                 | True         | 3        | False      | /data/S6_rlt6.0/boxsize2.5x2.5/Ns16/Nvpar48/Nmu9                |
| S6       | 6.0     | 3x1         | 48000.0  | 80000.0      | 16     | 48        | 9       | 0.02     | 1.4         | 21       | 247    | 15      | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     | "[10000, 28000]" | "[12000, 31000]" |                 | True         | 3        | True       | /data/S6_rlt6.0/boxsize3x1/Ns16/Nvpar48/Nmu9                    |
| S6       | 6.0     | 3x1.5       | 16840.3  | 28068        | 16     | 48        | 9       | 0.02     | 1.4         | 31       | 247    | 10      | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     | [11000]          | [12000]          |                 | True         | 4        | True       | /data/S6_rlt6.0/boxsize3x1.5/Ns16/Nvpar48/Nmu9                  |
| S6       | 6.0     | 3x1.5       | 14170.0  | 23618        | 16     | 48        | 9       | 0.02     | 1.4         | 31       | 247    | 10      | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     | [11000]          | [12000]          | 23618           | True         | 3        | True       | /data/S6_rlt6.0/boxsize3x1.5/Ns16/Nvpar48/Nmu9/Broken           |
| S6       | 6.0     | 3x2.5       | 6000.0   | 10000.0      | 16     | 48        | 9       | 0.02     | 1.4         | 51       | 247    | 6       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     | "[10000, 28000]" | "[12000, 31000]" |                 | True         | 3 & 4    | True       | /data/S6_rlt6.0/boxsize3x2.5/Ns16/Nvpar48/Nmu9                  |
| S6       | 6.0     | 3x3         | 7847.0   | 13083.0      | 16     | 48        | 9       | 0.02     | 1.4         | 61       | 247    | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     | [4000]           | [7000]           |                 | True         | 4        | True       | /data/S6_rlt6.0/boxsize3x3/Ns16/Nvpar48/Nmu9                    |
| S6       | 6.0     | 3x3         | 3454.8   | 5758.0       | 16     | 48        | 9       | 0.02     | 1.4         | 61       | 247    | 5       | cosine5   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     | [2000]           | [3000]           | "[1840, 1845]"  | True         | 4        | False      | /data/S6_rlt6.0/boxsize3x3/Ns16/Nvpar48/Nmu9/cosine5            |
| S6       | 6.0     | 3x3         | 5808.0   | 9680.0       | 16     | 48        | 9       | 0.02     | 1.4         | 61       | 247    | 5       | noise     | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     | [4000]           | [5000]           | "[1840, 1842]"  | True         | 3        | False      | /data/S6_rlt6.0/boxsize3x3/Ns16/Nvpar48/Nmu9/noise              |
| S6       | 6.0     | 3x5         | 4769.0   | 7958.0       | 16     | 48        | 9       | 0.02     | 1.4         | 101      | 247    | 3       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | True         | 4        | True       | /data/S6_rlt6.0/boxsize3x5/Ns16/Nvpar48/Nmu9                    |
| S6       | 6.0     | 4x1         | 30000.0  | 50000.0      | 16     | 48        | 9       | 0.02     | 1.4         | 21       | 329    | 20      | cosine2   | 1000             | "[24000, 26000]"                        | "[25500, 28000]"                         |                  |                  | 47084           | True         | 4        | True       | /data/S6_rlt6.0/boxsize4x1/Ns16/Nvpar48/Nmu9                    |
| S6       | 6.2     | 2x2         | 6000.0   | 10000.0      | 16     | 64        | 9       | 0.02     | 1.4         | 41       | 165    | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | True         | 2        | True       | /data/S6_rlt6.2/boxsize2x2/Ns16/Nvpar64/Nmu9                    |
| S6       | 6.2     | 3x3         | 14682.7  | 24475.0      | 16     | 48        | 9       | 0.02     | 1.4         | 61       | 247    | 5       | cosine2   | 2000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  | "[14473,16132]" | True         | 3        | True       | /data/S6_rlt6.2/boxsize3x3/Ns16/Nvpar48/Nmu9                    |
| S6       | 6.3     | 1x1         | 6000.0   | 10000.0      | 16     | 64        | 9       | 0.02     | 1.4         | 21       | 83     | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | False        | 0        | True       | /data/S6_rlt6.3/boxsize1x1/Ns16/Nvpar64/Nmu9                    |
| S6       | 6.4     | 3x3         | 12000.0  | 20000.0      | 16     | 48        | 9       | 0.02     | 1.4         | 61       | 247    | 5       | cosine2   | 1000             | "[500, 1500, 4000]"                     | "[1000, 2000, 5500]"                     |                  |                  |                 | False        | 0        | True       | /data/S6_rlt6.4/boxsize3x3/Ns16/Nvpar48/Nmu9                    |


### Important

* (*) marks files that are not included on GitHub and are only included on the NAS of TPV in Bayreuth
* (#) marks files that can be deleted to free space
* (+) marks files with additional diagnostics

## Folder Structure

The data is sorted by changed parameters in specific Folders. 
```
data
│
├── S6_rlt6.0
│  │
│  ├── boxsize1x1
│  │  │ 
│  │  ├── Ns12 (#)
│  │  │  │
│  │  │  ├── Nvpar16
│  │  │  │  └── Nmu6
│  │  │  │
│  │  │  ├── Nvpar32
│  │  │  │  └── Nmu6
│  │  │  │
│  │  │  ├── Nvpar48
│  │  │  │  └── Nmu9
│  │  │  │
│  │  │  └── Nvpar64
│  │  │     └── Nmu9
│  │  │
│  │  └── Ns16
│  │     │
│  │     ├── Nvpar16 (#)
│  │     │  └── Nmu9
│  │     │
│  │     ├── Nvpar32 (#)
│  │     │  └── Nmu9
│  │     │
│  │     ├── Nvpar48
│  │     │  ├── Nmu6 (#)
│  │     │  │
│  │     │  └── Nmu9
│  │     │     ├── dtim0.025
│  │     │     │
│  │     │     ├── krhomax0.70 (#)
│  │     │     │  └── Nmod11
│  │     │     │
│  │     │     ├── Nx43 (#)
│  │     │     └── Nx63 (#)
│  │     │
│  │     └── Nvpar64
│  │        ├── Nmu6 
│  │        └── Nmu9
│  │
│  ├── boxsize1.5x1.5
│  │  └── Ns16
│  │     └── Nvpar48
│  │        └── Nmu9
│  │
│  ├── boxsize2x1
│  │  └── Ns16
│  │     └── Nvpar48
│  │        └── Nmu9
│  │
│  ├── boxsize2x2
│  │  └── Ns16
│  │     └── Nvpar48
│  │        └── Nmu9
│  │
│  ├── boxsize2.5x2.5
│  │  └── Ns16
│  │     └── Nvpar48
│  │        └── Nmu9
│  │
│  ├── boxsize3x1
│  │  └── Ns16
│  │     └── Nvpar48
│  │        └── Nmu9
│  │
│  ├── boxsize3x1.5
│  │  └── Ns16
│  │     └── Nvpar48
│  │        └── Nmu9
│  │           └── Broken (#)
│  │
│  ├── boxsize3x2.5
│  │  └── Ns16
│  │     └── Nvpar48
│  │        └── Nmu9
│  │
│  ├── boxsize3x3
│  │  └── Ns16
│  │     └── Nvpar48
│  │        └── Nmu9
│  │           ├── cosine5 (+)
│  │           └── noise   (+)
│  │
│  ├── boxsize3x5
│  │  └── Ns16
│  │     └── Nvpar48
│  │        └── Nmu9
│  │
│  └── boxsize4x1
│     └── Ns16
│        └── Nvpar48
│           └── Nmu9
│
├── S6_rlt6.2
│  │
│  ├── boxsize2x2
│  │  └── Ns16
│  │     └── Nvpar64
│  │        └── Nmu9
│  │
│  └── boxsize3x3
│     └── Ns16
│        └── Nvpar48
│           └── Nmu9
│ 
├── S6_rlt6.3
│  └── boxsize1x1
│     └── Ns16
│        └── Nvpar64
│           └── Nmu9
│
└── S6_rlt6.4
   └── boxsize3x3
      └── Ns16
         └── Nvpar48
            └── Nmu9

```
   
### Files 

* `FDS`: restart file needed to restart run with gkw from current timestep (*)
* `FDS.dat`: database file of `FDS` with current status of the simulation after run with gkw
* `input.dat`: input database for simultaion with all parameters
* `input.out`: input output file (#)
* `kx_connect.dat`: parallel boundary connections database
* `output.dat`: output database of gkw
* `par.dat`: curvature and Coriolis functions database
* `perfloop.dat`: (no description found)
* `perform.dat`: performance database
* `status.txt`: status file of specific run generated with python script `\python\monitor_job.py` (not in every folder included)(#)
* `data.h5`: trimmed down version of `gkwdata.h5` for lower space consumption generated with shell script `\hdf5\hdf5_extract.sh` with added evaluated data (*), additional diagnostics marked with (+)
   <details><summary>Containing</summary>
   <p>

   ```
   data.h5
   │
   ├── diagnostic
   │  │
   │  ├── diagnos_fields 
   │  │  └── phi
   │  │
   │  ├── diagnos_fluxes 
   │  │  └── eflux_species01
   │  │
   │  ├── diagnos_grid
   │  │  └── lxn
   │  │
   │  ├── diagnos_moments
   │  │  └── xs_kyzero_dens (+), xs_kyzero_ene_par (+), xs_kyzero_ene_perp, (+)
   │  │
   │  └── diagnos_growth_freq 
   │     └── time
   │
   ├── grid
   │  └── xphi
   │
   └── evaluation
      └── derivative_stepsize, second_derivative_phi, zonalflow_potential, shearing_rate, shearing_rate_maximum,
          'derivative_dens' (+), 'second_derivative_dens' (+), 
          'derivative_energy_perp' (+), 'second_derivative_energy_perp' (+), 'derivative_energy_par' (+), 'second_derivative_energy_par' (+),
          'derivative_zonalflow_potential' (+), 'derivative_zonalflow_potential' (+)
   ```

   </p>
   </details>

* `gkwdata.h5`: raw data of the simulation (*) (#), additional diagnostics marked with (+)
   <details><summary>Containing</summary>
   <p>

   ```
   gkwdata.h5
   │
   ├── diagnostic
   │  │
   │  ├── diagnos_fields
   │  │  └── kxspec, kxvort, kyspec, kyvort, phi, spc
   │  │
   │  ├── diagnos_fluxes 
   │  │  └── EFlesr0001, eflux_species01, eflux_spectra, eflux_sup, eflux_xspec, flmgr01, pflux_species01, 
   │  │     pflux_spectra, pflux_sup, pflux_xspec, vflux_species01, vflux_spectra, vflux_xspec
   │  │
   │  ├── diagnos_grid
   │  │  └── intmu, intvp, lxn, lyn, mode_label, mphi, mphiw3, mrad_G, mrad_l, sgrid
   │  │
   │  ├── diagnos_growth_freq
   │  │  └── frequencies, frequencies_all_modes, growth, growth_rates_all_modes, time
   │  │
   │  ├── diagnos_mode_struct
   │  │  └── parallel
   │  │
   │  ├── diagnos_moments
   │  │  └── den01, den_spectra, ene01, ene_spectra, xs_kyzero_dens (+), xs_kyzero_ene_par (+), xs_kyzero_ene_perp, (+)
   │  │
   │  └── diagnos_rad
   │     └── prof_back
   │  
   ├── geom
   │  └── Bref, Bt_frac, D_eps, D_s, D_zeta, E_eps_s, E_eps_zeta, E_zeta_s, F, G, H_eps, H_s, H_zeta, I_eps, I_s, 
   │     I_zeta, J, Jacobian, K, NS, R, R0, Rref, Z, beta_eq, betaprime_eq, bmax, bmin, bn, eps, g_eps_eps,
   │     g_eps_s, g_eps_zeta, g_s_s, g_zeta_s, g_zeta_zeta, jfunh, jfunl, krnorm, kthnorm, lfun, poloidal_angle,
   │     q, s_grid, shat
   │
   ├── grid
   │  └── file_count, krho, krho_extended, krloc, kxrh, kzeta, time_fine, vperp, vpgr, xgr, xphi, yphi
   │
   ├── input
   │  │
   │  ├── collisions
   │  │  └── coll_freq, cons_type, en_scatter, ene_conservation, freq_override, friction_coll, lorentz, mass_conserve, 
   │  │     mom_conservation, nref, pitch_angle, rref, selfcollcon, tref, zeff
   │  │
   │  ├── control
   │  │  └── auto_restart, collisions, disp_par, disp_vp, disp_x, disp_y, dt_min, dtim, fac_dtim_est, flux_tube, 
   │  │     fluxtol, gamatol, ifluxtol, io_format, io_legacy, io_testdata, iperform_set, irun, laverage_dist_over_time,
   │  │     lflapv, lpar_vel_nl, lrestart_new_grid, ltrapping_arakawa, lverbose, matrix_format, max_gr, max_sec, 
   │  │     max_seconds, meth, method, min_gr, naverage, ncqtol, ndump_ts, neoclassics, nl_dtim_est, nlapar, nlbpar, 
   │  │     nlphi, non_linear, normalize_per_toroidal_mode, normalized, ntime, order_of_the_radial_scheme, 
   │  │     order_of_the_scheme, order_of_the_zf_scheme, parallel_boundary_conditions, radial_boundary_conditions, 
   │  │     read_file, restart_file_version, shift_metric, silent, spectral_radius, testing, uniform_mu_grid, vp_trap,
   │  │     zonal_adiabatic
   │  │
   │  ├── diagnostic
   │  │  └── apa3d, apc3d, bpa3d, bpc3d, cross_phase_timetrace, den3d, ene3d, field_fsa_kyzero, flux3d, imod_corr, 
   │  │     imod_field, kykxs_apar, kykxs_bpar, kykxs_j0_moments, kykxs_j1_moments, kykxs_moments, kykxs_phi, 
   │  │     lamplitudes, lcalc_corr, lcalc_energetics, lcalc_fluxes, lcalc_jdote, lcalc_jdote_fs, 
   │  │     lcalc_kinenergy_trappas, lcalc_stresses, lcalc_tot_energy, lcalc_zfshear, lcalc_zonal_evo, 
   │  │     lfinal_output, lfinvel, lfluxes_detail, lfluxes_em_spectra, lfluxes_spectra, lfluxes_vspace, 
   │  │     lfluxes_vspace_bpar, lfluxes_vspace_em, lfluxes_vspace_phi, lfrequencies, lgrowth_rates, lisl_follow,
   │  │     lmode_energy, lmpi_broken_io, lnonlin_transfer,lnonlin_transfer_fsa, lparallel_output, 
   │  │     lphi_diagnostics, lrad_entropy, lrad_field, lrad_kpar, lrad_moment, lrad_tint, lradial_entropy, 
   │  │     lradial_profile, lrotate_parallel, lvelspace_output, lwrite_output1, nflush_ts, 
   │  │     nmodepoints, nonlin_transfer_interval, npointsvel, out3d_interval, parallel_output_timestamps, phi3d, 
   │  │     psi_velspace, screen_output, spc3d, xs_kyzero_current, xs_kyzero_current2, xs_kyzero_dens, xs_kyzero_ene,
   │  │     xs_kyzero_ene_par, xs_kyzero_ene_perp, xs_kyzero_phi_ga_deltaf, xs_kyzero_phi_ga_fm, xs_phi, xy_apar, 
   │  │     xy_bpar, xy_current, xy_current2, xy_dens, xy_estep, xy_fluxes, xy_fluxes_bi, xy_fluxes_bpar, 
   │  │     xy_fluxes_em, xy_fluxes_fsa, xy_fluxes_k, xy_fluxes_p, xy_fluxes_v, xy_phi, xy_slice_ipar, xy_spec, 
   │  │     xy_temp, xy_vort, zeta_velspace, zevo_detail, zevo_xy_spec, zonal_scale_3d
   │  │
   │  ├── eiv_integration
   │  │  └── freq, growthrate, mat_vec_routine, max_iterations, nr_column_vec, number_eigenvalues, tolerance, 
   │  │     type_extraction, type_solver, which_eigenvalues
   │  │
   │  ├── finite_rho_parallel
   │  │  └── lflux_rhostar, lnonlinear_rhostar, ltrapdf_rhostar, lvdgrad_phi_fm_rhostar, lvdgradf_rhostar, 
   │  │     lve_grad_fm_rhostar, s_average
   │  │
   │  ├── geom
   │  │  └── N_shape, R0_loc, Zmil, beta_rota_miller, beta_rota_miller_type, c, c_prime, curv_effect, dRmil, dZmil, 
   │  │     delta, eps, eps_type, eqfile, geom_type, gradp, gradp_type, kappa, prof_type, q, qprof_coef, s, s_prime, 
   │  │     sdelta, shat, signB, signJ, skappa, square, ssquare
   │  │
   │  ├── grid
   │  │  └── lx, mumax, n_mu_grid, n_procs_mu, n_procs_s, n_procs_sp, n_procs_vpar, n_procs_x, n_s_grid, n_trapped, 
   │  │     n_vpar_grid, n_x_grid, nmod, non_blocking_vpar, nperiod, number_of_species, nx, psih, psil, vpmax
   │  │
   │  ├── gyroaverage
   │  │  └── blending_order, consistent_long_wave, gyro_average_electrons, gyro_average_ions, mk_gyro_av_hermitian, 
   │  │     n_gav_bound_ex, n_points_ga, orb_polarize, parallel_mod, use_conj
   │  │
   │  ├── header
   │  │  └── compiler, gkw_executable_name, gkw_version, number_of_processors, timestamp
   │  │
   │  ├── krook
   │  │  └── bwidth, gamkpre, gammab, gammak, krook_option, nlbound, nlkrook
   │  │
   │  ├── linear_term_switches
   │  │  └── apply_on_imod, idisp, lampere, lbpar, lg2f_correction, lneo_equil, lneo_rad, lneo_trap, lneoclassical,
   │  │     lneorotsource, lpoisson, lpoisson_zf, ltrapdf, lvd_grad_phi_fm, lvdgradf, lve_grad_fm, lvpar_grad_df, 
   │  │     lvpgrphi, neo_equil_parse_sp_seq
   │  │
   │  ├── mode
   │  │  └── chin, ikxspace, kr_type, krhomax, krrho, kthrho, mode_box, n_spacing, no_drive_of, no_transfer_from, 
   │  │     no_transfer_to, rkxspace
   │  │
   │  ├── rotation
   │  │  └── cf_drift, cf_qncheck, cf_trap, cf_upphi, cf_upsrc, coriolis, shear_profile, shear_rate, t_shear_begin, 
   │  │     toroidal_shear, vcor
   │  │
   │  ├── source_time
   │  │  └── dsfr, gauss_source_median, gauss_source_stdev, mod_freq, source_profile, source_time_ampl, 
   │  │     source_wave_number
   │  │
   │  ├── spcgeneral
   │  │  └── Ls, adiabatic_electrons, amp_imod, amp_imod_imag, amp_init, amp_zon, amp_zon_imag, beta, beta_ref, 
   │  │     beta_type, betaprime_ref, betaprime_type, energetic_particles, finit, finit_imod, icrh_params, 
   │  │     imod_init, init_coef, isl_mode, isl_rot_freq, isl_shear, lfinit_radial_dirichlet, mode_persist, 
   │  │     n_quench, psi_0, quench_modes, quench_switch, rhostar, tear_zero_epar, tearingmode, tm_start, 
   │  │     vpar_mean, wstar
   │  │
   │  ├── species01
   │  │  └── background, dens, dens_prof_coef, dens_prof_type, mass, param, rln, rlt, rlt_gauss, temp, temp_prof_coef,
   │  │     temp_prof_type, uprim, z
   │  │
   │  └── species02
   │     └── background, dens, dens_prof_coef, dens_prof_type, mass, param, rln, rlt, rlt_gauss, temp, temp_prof_coef,
   │        temp_prof_type, uprim, z
   │
   └── restart
      └── dtim, nt_complete, nt_remain, time_complete
   ``` 
   
   </p>
   </details>
