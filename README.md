# Size Convergence of the ExB Staircase Pattern in Flux Tube Simulations of Ion Temperature Gradient Driven Turbulence


![alt text](/pictures/Comparison/Boxsize/S6_rlt6.0_boxsize1-2-3-4x1-1.5-2-2.5-3-5_Ns16_Nvpar48_Nmu9_wexb_comparison.png)



## Content

1. [Introduction](#introduction)
2. [Abstract](#abstract)
3. [Journal](#journal)
4. [Literature](#literature)



## Introduction

This repository is focused on my work for my Bachelor Thesis about the topic of the size convergence of ExB Staircase Pattern with the box size. This Thesis is based on the works of Rath,F. and Peeters,A. G. and Buchholz,R. and Grosshauser,S. R. and Migliano,P. and Weikl,A. and Strintzi,D.

* [GKW-Code](https://bitbucket.org/gkw/gkw/wiki/Home)
* [Bachelor-Thesis](bachelorthesis/bachelorthesis.pdf)
* [Brief-Communication](briefcommunication/briefcommunication.pdf)

## Abstract

Ion temperature gradient driven turbulence (ITG) close to marginal stability exhibits zonal flow pattern formation on mesoscales, so-called $E\times B$ staircase structures. Such pattern formation has been observed in local gradient-driven flux-tube simulations as well as global gradient-driven and global flux-driven studies.

To reduce the computational effort for the simulations lower input parameter of GKW (Gyro Kinetic Workshop) were tested to find the optimum of minimum resolution for the performed simulations.

For convenience, a ```python``` script [```slurm_monitor.py```](python/slurm_monitor.py) was written to monitor the simulation on the ```btrzx1```-cluster and start/restart until the completion criterion is fulfilled.

Furthermore, it is shown by multiple box size convergence scans that a mesoscale pattern  size of $\sim 57-76\,\rho$ is inherent to ITG driven turbulence with Cyclone Base Case parameters in the local limit. This outcome also implies that a typical scale for avalanche-like transport is inherent to ITG driven turbulence.

## Journal
The work on the thesis is documented in from of a journal and to keep track of all changes [Source Control from GitHub](https://github.com/ManeLippert/Bachelorthesis-ZonalFlows/commits/main) was used.

<details><summary>Journal</summary>
<p>

* <details><summary>2022</summary>
  <p>

  * <details><summary>March</summary>
    <p>

    * <details><summary>24.03.2022 &nbsp; Starting Meeting</summary>
      <p>

      # Starting Meeting

      #### Thursday 24.03.2022 from 14:00 to 14:25 with Florian Rath and Arthur Peeters

      ### Discussion how to begin the work for bachelor thesis:

      * Start with reproduction of result in [[1]](/literature/Peeters%2C%20Rath%2C%20Buchholz%20-%20Gradient-driven%20flux-tube%20simulations%20of%20ion%20temperature%20gradient%20turbulence%20close%20to%20the%20non-linear%20threshold%20(Paper%2C%202016).pdf) with help of [gkw](/gkw/)
      * Because of the long runtime of the code firstly we will look only in one direction in the velocity space
      * After that small steps in all directions for better understanding of the structure and to find a minimal resolution for the best results 
      * Furthermore increase box size and search for convergence of the wavelength in zonal flows
      * There will be interpretation needed to clarify simplification steps in code

      ### Thesis
      * Work in English or German > will do it in English
      * Continues writing is better than everything in the end

      </p>
      </details>

    </p>
    </details>

  * <details><summary>April</summary>
    <p>

    * <details><summary>07.04.2022 &nbsp; Kurs "Schreiben einer MINT-Arbeit"</summary>
      </p>

      # Kurs "Schreiben einer MINT-Arbeit"

      #### Dienstag 07.04.2022 von 9:00 bis 15:00

      ## Inhalt
      * [Feststellung des Schreibtyps](#feststellung-des-schreibtyps)
      * [Störfaktorem](#störfaktoren)
      * [Phasen des Schreibprozesses](#phasen-des-schreibprozesses)
      * [Fragestellung/Forschungsfrage](#fragestellungforschungsfrage)
      * [Gliederung](#gliederung)
      * [Materialen](#materialen)
      * [Rohtext](#rohtext)
      * [Wissenschaftlicher Schreibstil](#wissenschaftlicher-schreibstil)
      * [Illustrationen](#illustrationen)
      * [Zitieren](#zitieren)
      * [Beleg im Text](#beleg-im-text)
      * [Methoden der Organisation und Planung](#methoden-zur-organisation-und-planung)

      ## Feststellung des Schreibtyps

      ![FragenSchreibtyp1](/pictures/HowToMINT/Schreibtypentest-1.png)
      ![FragenSchreibtyp2](/pictures/HowToMINT/Schreibtypentest-2.png)
      ![FragenSchreibtypAuswertung1](/pictures/HowToMINT/Schreibtypentest-3.png)
      ![FragenSchreibtypAuswertung2](/pictures/HowToMINT/Schreibtypentest-4.png)

      ## Störfaktoren

      * **Zeitdiebe** &rarr; Prokrastination am Handy?
      * **Schreiborte** &rarr; Feststellen wo die besten Schreiborte für einen sind &rarr; Draußen bei schönen Wetter
      * **Schreibzeiten** &rarr; Morgen, Nachmittags oder Abends &rarr; Nachmittags oder Abends

      ## Phasen des Schreibprozesses
      1. Orientierung und Planung
      2. Strukturieren, gliedern, forschen/lesen
      3. Material auswerten, Rohfassung schreiben
      4. Überarbent und Feedback einholen
      5. Schlusskorrektur und Abgabe

      ## Fragestellung/Forschungsfrage

      Grenzt Thema ein und leitet fokussiert durch die Arbeit

      ![Forschungsfrage1](/pictures/HowToMINT/AB1_Forschungsfrage-1.png)
      ![Forschungsfrage2](/pictures/HowToMINT/AB1_Forschungsfrage-2.png)

      ![ForschungsfrageHandout](/pictures/HowToMINT/Handout_Forschungsfrage.png)


      ## Gliederung

      * **Einleitung** &rarr; Hinführung, Problemstellung. Fragestellung (thematisieren), Methodik, Aufbau, Hauptergebnisse
      * **Methoden** &rarr; Zustandekommen der Ergebnisse, Grund für Glaubwürdigkeit (Auch Materialen)
      * **Ergebnisse** &rarr; Ausformulierung und Darstellung
      * **Diskussion** &rarr; Bezug auf Ergebnisse, dann breiter Fokus (Rückbezug zur Problemstellung)

      ## Materialen
      Quellen und Literatur frühzeitig dokumentieren (auch Anmerkungen möglich)

      ## Rohtext
      * Erstefassung eines Textes
      * Noch ungeschliffen
      * Macht as den Gedanken etwas Konkretes
      * Nimmt den Druck alles beim ersten Schreiben perfekt zu machen
      * Liefert Grundlage für weitere Schritte
      * Mehrfache Überarbeitungen machen den Rohtext zu einen abgereiften Text

      ## Wissenschaftlicher Schreibstil

      * Sachlich und Neutral
      * Logische Argumentation und Aufbau (roter Faden) &rarr; Forschungsfrage
      * Überprüfbarkeit und Nachvollziehbarkeit (Zitation)
      * Korrekte Verwendung von Fachbegriffen
      * Einheitlichkeit

      <br />

      ![Schreibstil](/pictures/HowToMINT/AB2_Schreibstil_%C3%9Cbung.png)

      ## Illustrationen
      ![Illu](/pictures/HowToMINT/Handout_Illustrations.png)

      ## Zitieren

      ### Faustregel
      1. Überhaupt zitieren
      2. Einheitlich zitieren
      3. Vorgaben beachten

      <br />

      Es gibt aber nicht den einen Zitierstil. Dieser kann sich von Fach zu Fach ändern.

      ### **WICHTIG**
      * Nachprüfbarkeit und Nachvollziehbarkeit
      * Einwandfreies zitieren &rarr; Ausdruck für wissenschaftliche Sorgfalt
      * Nachweis über über eigenständige Leitung &rarr; Trennung der Aussagen
      * Lesbarkeit &rarr; Mehr wissenschaftliche Form

      ### 1. Wörtliches/Direktes Zitat
      * Wörtliche Übernahme von Textpassagen, Sätzen, Satzteilen und Ausdrücken
      * Beginnt und endet mit Anführungszeichen
      * Längere Zitate werden i.d.R. eingerückt
      * Buchstabliche Genauigkeit 
      * Evtl. kursive Schrift, kleinere Schriftart, Absatz mit Einrückung und einzeiliger Abstand

      ### 2. Paraphrase/Indirektes Zitat
      * Sinngemäße Übernahme fremder Gedanken/Aussagen mit eigenen Worten
      * Ohne Anführungszeichen
      * Umfang muss eindeutig erkennbar sein 
      * Eventuell Zusatz "vgl."

      ### Beleg im Text
      &rarr; Verweis wird in Klammern hinter dem Zitat angefügt, gefolgt von einem Punkt: 

      &nbsp;  &nbsp; &nbsp;.....(Vgl. Eco, 2010, S.204). (**Vor dem Punkt**)

      &rarr; Wenn Autoren explizit erwähnt wurden, folgt die Quelle direkt hinter dem Namen: 

      &nbsp;  &nbsp; &nbsp;.....Eco (2010, S.204)

      &rarr; Verweis mit Fußnote. Jede Fußnote beginnt mit einem Großbuchstaben und endet mit einem Punkt. Zahl der Fußnote folgt hinter dem Punkt

      &nbsp;  &nbsp; &nbsp;.....xyz.³

      ___
      &nbsp;  &nbsp; &nbsp;³Vgl. Eco, 2010, S.204.

      ## Methoden zur Organisation und Planung

      ![Orga1](/pictures/HowToMINT/Methodenhandout_WS%20Orga%20und%20Planen-1.png)
      ![Orga2](/pictures/HowToMINT/Methodenhandout_WS%20Orga%20und%20Planen-2.png)
      ![Orga3](/pictures/HowToMINT/Methodenhandout_WS%20Orga%20und%20Planen-3.png)

      </p>
      </details>

    </p>
    </details>

  * <details><summary>May</summary>
    <p>

    * <details><summary>05.05.2022 &nbsp; Start with Bachelor Work</summary>
      <p>

      # Start with Bachelor Work

      #### Thursday 24.03.2022 from 14:00 to 14:27 with Florian Rath and Arthur Peeters

      ### Discussion on how to run the code:

      #### Login:

      * Login on local machine through ```x2go``` because ```ssh``` is too slow. 
      * When someone uses login through ```ssh``` the command line is shrunk down to a limited amount of executables that results in no ```make``` command. To get full access to the command line one has too ```ssh``` to ```bpptx```

      #### Cluster:

      * ```btrzx1``` is easier to run code 
      * ```btrzx3``` could cause problems with the nodes but is more efficient than ```btrzx1```

      Run code first on ```btrzx1``` with [```bashrc_btrzx1```](/gkw/run_btrzx1/bashrc_btrzx1) (loads all modules for ```GKW```) with jobmanager ```SLURM``` (started with ```sbatch```) and jobscript [```jobscript_btrzx1_simple```](/gkw/run_btrzx1/jobscript_btrzx1_simple).

      #### Sync Files:

      From local to remote machine
      ```
      scp -r Bachelorthesis-ZonalFlows/gkw/ user@btrzx1-1.rz.uni-bayreuth.de:gkw/
      ```
      From remote to local
      ```
      scp -r user@btrzx1-1.rz.uni-bayreuth.de:gkw/ Bachelorthesis-ZonalFlows/gkw/ 
      ```

      on Linux account just use ```git``` protocol

      ### What to do first:

      * Use test cases with adiabatic electrons
      * Work with spectral and non-spectral (cheaper, but steps in heat production not reproducible) and compare the time duration
      * In [paper](/literature/Peeters%2C%20Rath%2C%20Buchholz%20-%20Comparison%20of%20gradient%20and%20flux%20driven%20gyro-%0Akinetic%20turbulent%20transport%20(Paper%2C%202016).pdf) they used spectral 
      * Compare spectral outcome with [paper](/literature/Peeters%2C%20Rath%2C%20Buchholz%20-%20Comparison%20of%20gradient%20and%20flux%20driven%20gyro-%0Akinetic%20turbulent%20transport%20(Paper%2C%202016).pdf)
      * Verify the decrease of turbulence and heat flux on work point (condition of this bachelor thesis)

      </p>
      </details>

    * <details><summary>10.05.2022 &nbsp; First Day in the Office in Bayreuth</summary>
      <p>

      # First Day in the Office in Bayreuth

      #### Thusday 10.05.2022 from 10:00 to 17:30

      ### First Run with gkw
      For the first run I used the [input.dat.minimum](https://github.com/ManeLippert/Bachelorthesis-ZonalFlows/blob/main/gkw/doc/input.dat.minimum) that gaves me the examination files in the ```~/gkw/run``` directory. For futher examination I will use ```python``` on my local machine.

      ### Discussion with Florian Rath

      * Run ```gkw``` with configuration (S6) from [[1]](/literature/Peeters%2C%20Rath%2C%20Buchholz%20-%20Gradient-driven%20flux-tube%20simulations%20of%20ion%20temperature%20gradient%20turbulence%20close%20to%20the%20non-linear%20threshold%20(Paper%2C%202016).pdf) page 2

      Use [```cyclone```](https://github.com/ManeLippert/Bachelorthesis-ZonalFlows/blob/main/gkw/doc/input/cyclone) as basis ```input.dat``` and change parameter according (S6)

      * Save data as ```hdf5``` (8 times more compact than ```ASCII```). ```python``` can read files easily 

      * As diagnostic run ```xy_phi``` to get data from [[1]](/literature/Peeters%2C%20Rath%2C%20Buchholz%20-%20Gradient-driven%20flux-tube%20simulations%20of%20ion%20temperature%20gradient%20turbulence%20close%20to%20the%20non-linear%20threshold%20(Paper%2C%202016).pdf) page 8 pictures

      ```
      !------------------------------------------------------------------------------------------------------------------------
      &CONTROL
      zonal_adiabatic = .true.,               !If zonal flows corrections included for adiabiatic electrons       (default = F)

      order_of_the_zf_scheme = 'sixth_order'  !Use a different finite-differences scheme for (default = order_of_the_scheme)

      D      = disp_par = 1.0                 !(Hyper) dissipation coefficient for parallel derivatives.          (default=0.2)
      D_vpar = disp_vp  = 0.2                 !(Hyper) dissipation coefficient for parallel velocity space        (default=0.2)
      D_x    = disp_x   = 0.1                 !(Hyper) dissipation coefficient in perpendicular x direction       (default=0.0)
      D_y    = disp_y   = 0.1                 !(Hyper) dissipation coefficient in perpendicular y direction       (default=0.0)

      io_format = 'hdf5'                      ! Use 'ascii' to output all data as formatted text files      (default = 'mixed')
                                              !     'binary' to output all data as unformatted binary files
                                              !     'mixed' to output some binary and mostly text files
                                              !     'hdf5' to output a single HDF5 file (needs compilation with HDF5 libraries)
                                              !     'hdf5+ascii' to output a single HDF5 file and duplicate 1D and 2D data to
                                              !         formatted text files.
                                              !     'none' to output no data at all.
      /
      !------------------------------------------------------------------------------------------------------------------------
      &GRIDSIZE
      N_m    = NMOD        = 21               !Number of binormal modes - do not interact for linear runs
      N_x    = NX          = 83               !Number of radial wave vectors / points: needs to be an odd number for spectral
      N_s    = N_s_grid    = 16               !Number of grid points along the field line
      N_vpar = n_vpar_grid = 64               !Number of grid points for parallel velocity (must be even)
      N_mu   = N_mu_grid   = 9                !Total number of magnetic moment grid points
      /
      !------------------------------------------------------------------------------------------------------------------------
      &MODE
      mode_box = .true.,                      !Determines if there is a 2D grid of ky,kx. if true use nperiod = 1 (default = F)
                                              !If nperiod = 1 and mode box = .true. the kx modes will be coupled.
      krhomax = 1.4,                          !For mode_box, this is the maximum k_theta rho_i (ky) on the grid.(default = 0.0)
                                              !For nmod>1, modes are equidistantly spaced from 0.0 to to krhomax.
                                              !k_perp is evaluated on the low field side of the outboard midplane.
                                              !rho_i evaluated on the flux surface at the major radius of the magnetic axis.
                                              !Note that other codes may normalise the thermal velocity differently
                                              !which can correspond to  gkw k_theta that are a factor sqrt(2) greater.
      /
      !------------------------------------------------------------------------------------------------------------------------
      &SPECIES
      rlt = 6.0
      /
      !------------------------------------------------------------------------------------------------------------------------
      &GEOM
      GEOM_TYPE = 'circ'                      !Switch for the metric: 's-alpha', 'circ', 'miller', 'fourier' or 'chease'   
                                              !(default = 's-alpha')
      /
      !------------------------------------------------------------------------------------------------------------------------
      &DIAGNOSTIC
      xy_phi = .true.                         ! Electrostatic potential in perpendicular plane at LFS midplane    (default = T)
      /
      !------------------------------------------------------------------------------------------------------------------------
      &LINEAR_TERM_SWITCHES                   
      v_d = idisp = 1                         !Select between dissipation schemes in finite differences 
      /
      ```
      </p>
      </details>

    * <details><summary>11.05.2022 &nbsp; Run for Standard Resolution 6th order (S6)</summary>
      <p>

      # Run for Standard Resolution 6th order (S6)

      #### Wednesday 11.05.2022 9:45 to 13:30

      ### New Input file

      [```input_S6_rtl6.dat```](../data/S6_rlt6.0/Nsgrid16_Nvpargrid64_Nmugrid9/input.dat)

      On ```btrzx1``` the maximal available processors are 32 so that you have to determine additional values. Furthermore ```gkw``` needs time to write files and the maximal runtime should be 15min less than the ```walltime```. On ```btrzx1``` the ```walltime``` is set to 24h (maximum duration). Lastly I set the parameter for the timesteps for writing checkpoint files in ```ndump_ts```.

      #### Conditions:
      * ```N_procs_mu``` < ```N_mu_grid```
      * ```N_procs_vpar``` * ```N_procs_s``` != 32
      * ```max_seconds``` = ```walltime``` - 900


      ```
      !------------------------------------------------------------------------------------------------------------------------
      &CONTROL
      zonal_adiabatic = .true.,               !If zonal flows corrections included for adiabiatic electrons       (default = F)

      order_of_the_zf_scheme = 'sixth_order'  !Use a different finite-differences scheme for (default = order_of_the_scheme)

      D      = disp_par = 1.0                 !(Hyper) dissipation coefficient for parallel derivatives.          (default=0.2)
      D_vpar = disp_vp  = 0.2                 !(Hyper) dissipation coefficient for parallel velocity space        (default=0.2)
      D_x    = disp_x   = 0.1                 !(Hyper) dissipation coefficient in perpendicular x direction       (default=0.0)
      D_y    = disp_y   = 0.1                 !(Hyper) dissipation coefficient in perpendicular y direction       (default=0.0)

      io_format = 'hdf5'                      ! Use 'ascii' to output all data as formatted text files      (default = 'mixed')
                                              !     'binary' to output all data as unformatted binary files
                                              !     'mixed' to output some binary and mostly text files
                                              !     'hdf5' to output a single HDF5 file (needs compilation with HDF5 libraries)
                                              !     'hdf5+ascii' to output a single HDF5 file and duplicate 1D and 2D data to
                                              !         formatted text files.
                                              !     'none' to output no data at all.

      ndump_ts=500                   !Number of large timesteps between writing of checkpoint DMP files    

      max_seconds = 85500            ! 24h = 86400s 15min = 900s -> 85500
      /
      !------------------------------------------------------------------------------------------------------------------------
      &GRIDSIZE
      N_m    = NMOD        = 21               !Number of binormal modes - do not interact for linear runs
      N_x    = NX          = 83               !Number of radial wave vectors / points: needs to be an odd number for spectral
      N_s    = N_s_grid    = 16               !Number of grid points along the field line
      N_vpar = n_vpar_grid = 64               !Number of grid points for parallel velocity (must be even)
      N_mu   = N_mu_grid   = 9                !Total number of magnetic moment grid points

      N_procs_mu   = 3                        !As above, but for mu                                              
      N_procs_vpar = 8                        !As above, but for vpar (>1 only works if vp_trap = 0)             
      N_procs_s    = 4                        !As above, but for s
      /
      !------------------------------------------------------------------------------------------------------------------------
      &MODE
      mode_box = .true.,                      !Determines if there is a 2D grid of ky,kx. if true use nperiod = 1 (default = F)
                                              !If nperiod = 1 and mode box = .true. the kx modes will be coupled.
      krhomax = 1.4,                          !For mode_box, this is the maximum k_theta rho_i (ky) on the grid.(default = 0.0)
                                              !For nmod>1, modes are equidistantly spaced from 0.0 to to krhomax.
                                              !k_perp is evaluated on the low field side of the outboard midplane.
                                              !rho_i evaluated on the flux surface at the major radius of the magnetic axis.
                                              !Note that other codes may normalise the thermal velocity differently
                                              !which can correspond to  gkw k_theta that are a factor sqrt(2) greater.
      /
      !------------------------------------------------------------------------------------------------------------------------
      &SPECIES
      rlt = 6.0
      /
      !------------------------------------------------------------------------------------------------------------------------
      &GEOM
      GEOM_TYPE = 'circ'                      !Switch for the metric: 's-alpha', 'circ', 'miller', 'fourier' or 'chease'   
                                              !(default = 's-alpha')
      /
      !------------------------------------------------------------------------------------------------------------------------
      &DIAGNOSTIC
      xy_phi = .true.                         ! Electrostatic potential in perpendicular plane at LFS midplane    (default = T)
      /
      !------------------------------------------------------------------------------------------------------------------------
      &LINEAR_TERM_SWITCHES                   
      v_d = idisp = 1                         !Select between dissipation schemes in finite differences 
      /
      ```

      In the input file is also more options that are provided by Florian Rath. For example the option that ```gkw``` will automatically write restart files and additional DIAGNOSTICS.

      ### Jobscript

      #### Conditions:
      * ```SBATCH --nodes=N_procs_mu*N_procs_vpar*N_procs_s``` = 96
      * ```SBATCH --ntasks-per-node=32```
      * ```SBATCH --nodes=3``` = 32 * 3 = 96
      * ```SBATCH --time=0-24:00:00```

      [```jobscript_btrzx1_S6```](../gkw/btrzx1/jobscript_btrzx1_S6)

      </p>
      </details>

    * <details><summary>12.05.2022 &nbsp; Discussion about Resolution & Run for (S6) with rtl=6.3</summary>
      <p>

      # Discussion about Resolution

      #### Thursday 24.03.2022 from 14:00 to 14:25 with Florian Rath and Arthur Peeters

      ### Minimum Values

      WIP so we will try to find the best minimum resolution

      * ```N_s_grid``` = 12
      * ```N_vpar_grid``` = 16 or 32
      * ```N_mu_grid``` = 6

      Numeric dissipation gains with smaller scales of resolution that could cause the **lost** of zonal flows

      ### ```Python``` Program

      * Write ```python``` program to evaluate the ```xy_phi``` diagnostics and symbolize 'Scherrrate' and heat flux
      * Learn how to evaluate ```h5``` files

      # Run for (S6) with rtl=6.3

      [```input_S6_rtl6.3.dat```](../data/S6_rtl6.3/input.dat)

      [```jobscript_btrzx1_S6```](../data/S6_rtl6.3/jobscript_btrzx1)

      </p>
      </details>

    * <details><summary>16.05.2022 &nbsp; Writing of useful shell scripts</summary>
      <p>

      # Writing of useful shell scripts

      #### Monday 16.05.2022 13:15 to 23:00

      ## Shell Scripts
      * [```ssh_btrzx1```](../ssh/ssh_btrzx1.sh) turns automatically the vpn connection on and connects to ```btrzx1-1.rz.uni-bayreuth.de```

      * [```ssh_copy```](../ssh/ssh_copy.sh) useful copy script to copy files from remote to local or in the other direction

      </p>
      </details>

    * <details><summary>20.05.2022 &nbsp; Discussion about evaluation of the shearing rate $\omega_{\mathrm{E \times B}}$</summary>
      <p>

      # Discussion about evaluation of the shearing rate

      #### Friday 20.05.2022 12:00 to 12:15 with Florian Rath and Arthur Peeters

      ## Coordinate

      The coordinate ```x``` is in the ```h5```-file marked as ```xphi``` and is the radial coordinate

      ## Derivative

      The derivative is periodic which means that at the start point $f_0$ the other two points for derivative would be $f_{N}$ and $f_1$ and at the end point $f_{N}$ the other two points would be $f_{N-1}$ and $f_0$.\
      \
      That concludes to the formula:\
      \
      Start: $\frac{f_1 - 2 \cdot f_0 + f_N}{h^2}$\
      \
      Middle: $\frac{f_{i+1} - 2 \cdot f_i + f_{i-1}}{h^2}$\
      \
      End: $\frac{f_{0} - 2 \cdot f_N + f_{N-1}}{h^2}$

      ## Additional Diagnostic

      Use fourier spetrum as additional diagnostic to evaluate the shearing rate $\omega_{\mathrm{E \times B}}$ like in Fig 5a in [[1]](https://doi.org/10.1063/1.4961231)

      </p>
      </details>

    </p>
    </details>

  * <details><summary>June</summary>
    <p>

    * <details><summary>08.06.2022 &nbsp; Resolution, Folder Structure & Comparison of Resolution</summary>
      <p>

      # Resolution, Folder Structure and Comparison of Resolution

      #### Wednesday 08.06.2022

      ## Resolution

      Best resolution: 

      ```Nsgrid = 16```, ```Nvpar = 48```, ```Nmugrid = 9```

      Possible Variations: 

      * ```krhomax = 0.70 | nmod = 11```
      * ```krhomax = 1.05 | nmod = 16```
      * ```nx = 63```, ```nx = 43```

      ## New Folder Structure

      Every change in ```input.dat``` gets it own folder and the evaluation notebook write changes in picture name. Furthermore the notebook will write with python the picture folder.

      ## Comparison of Resolution

      * ```Nsgrid = 12/16``` | ```Nvpargrid = 64```, ```Nmugrid = 9``` | ```Nvpargrid = 48```, ```Nmugrid = 9```
      * ```Nvpargrid = 64/48/32/16``` | ```Nsgrid = 16```, ```Nmugrid = 9``` 
      * ```Nmugrid = 6/9``` | ```Nvpargrid = 64```, ```Nsgrid = 16``` | ```Nvpargrid = 48```, ```Nsgrid = 16```

      </p>
      </details>

    * <details><summary>09.06.2022 &nbsp; Meeting to increase Boxsize radially</summary>
      <p>

      # Meeting to increase Boxsize radially

      #### Thursday 09.06.2022 14:00 to 14:30 with Florian Rath and Arthur Peeters

      ## Change Timestep
        Set ```dtim = 0.02``` to ```dtim = 0.025``` and compare outcome with $\delta t$. The graph should decrease vor ```dtim```.
      ```dtim``` is a timestep measured with gkw.

      ## Final Resolution

      ```Nsgrid = 16```, ```Nvpar = 48```, ```Nmugrid = 9```

      ## Increase Boxsize radially

      Change following variables according to increase factor $N$:

      * ```ikx_space_N``` $= 5 * N$
      * ```nx_N``` $= [($ ```nx_1``` $-1 ) * N ] +1$ 

      Boxsize 1x1: `nx_1` = 83, &nbsp; `ikx_space_1` = 5\
      Boxsize 2x1: `nx_2` = 165, `ikx_space_2` = 10\
      Boxsize 3x1: `nx_3` = 247, `ikx_space_3` = 15\
      Boxsize 4x1: `nx_4` = 329, `ikx_space_4` = 20

      </p>
      </details>

    * <details><summary>15.07.2022 - 29.07.2022 &nbsp; Work on Restart Script & Run for increased Boxsize</summary>
      </p>

      # Work on Restart Script

      Lots of work for the rest of the month gone into the development of the restart script [`slurm_monitor.py`](/python/slurm_monitor.py) to 
      tackle the problem of dealing everyday with restarts of the code due to some wall time of the cluster btrzx1.
      For that the script is developed as python3 script that only needs built in python modules to ensure running on every system. 
      The script itself looks in a specific time interval if the job is running, pending or needs to be started und das this routine until a defined
      timestep is reached all by load the output of SLURM Job Manager with `squeue` and analyse the output.
      The Core build could be adopted for diffent jobmanager as well the script is build variable enough to ensure the changing of the inportant values.

      # Run of increased boxsize

      The Rest of the time was waiting for the simulation for Boxsize 4x1 to be complete
      </p>
      </details>

    </p>
    </details>

  * <details><summary>July</summary>
    <p>

    * <details><summary>05.07.2022 &nbsp; Meeting to increase Boxsize binormal </summary>
      <p>

      # Meeting to increase Boxsize binormal

      #### Thursday 05.07.2022 14:00 to 14:30 with Florian Rath and Arthur Peeters

      ## Increase Boxsize binormal

      Change following variables according to increase factor $N$:
      * ```ikx_space_N``` $= 5$
      * ```nx_N``` $= [($ ```nx_1``` $-1 ) * N ] +1$ 
      * ```nmod_N``` $= [($ ```nmod_1```$ -1) * N ] +1$

      Boxsize 1x1: `nx_1` = 83, &nbsp; `nmod_1` = 21\
      Boxsize 2x2: `nx_2` = 165, `nmod_2` = 41\
      Boxsize 3x3: `nx_3` = 247, `nmod_3` = 61\
      Boxsize 4x4: `nx_4` = 329, `nmod_4` = 81

      </p>
      </details>

    * <details><summary>06.07.2022 - 29.07.2022 &nbsp; Problems with hdf5-file & Further work on restart script</summary>
      <p>

      ## Problems with hdf5-file

      hdf5 files have to be closed every time you are done with processing data. Otherwise the file gets curupted and the data gets lost because
      only the programm that opens the hdf5 file can close it. This behaviour results in lots of problems on the server because of the storage limit
      on btrzx1 GKW got stopped and the file remained open.

      ## Further Work on restart script

      Because of that the restart script now features a backup option to safe data between successful runs und can restore it after error.
      As additonal the restart script now can write the job name into the jobscript file, has timestaps for each new status update, 
      writes outputs in `status.txt` and sends mails at the start and the end of on total run.

      </p>
      </details>

    </p>
    </details>

  * <details><summary>August</summary>
    <p>

    * <details><summary>06.08.2022 - 16.08.2022 &nbsp; Evaluate Data</summary>
      <p> 

      ## Evaluate Data
      To make sure every simulations has no turbulence a fourier plot of fourer mode 1 to 5 (in Plots $k_1$ to $k_5$) will in the time domain be made. 
      It has shown that the mode with a value of $\omega_{\mathrm{E \times B}, max} \sim 0.20$ is also the wavelength thats converges 
      with the boxsize. So if the mode $k_3$ is at $\omega_{\mathrm{E \times B}, max} \sim 0.20$ we know when the other modes are nearly zero
      that in the boxsize the 3 times wavelength converges with the boxsize.

      ## Results

      Boxsize 1x1: $k_1$\
      Boxsize 2x1: $k_2$\
      Boxsize 2x2: $k_2$\
      Boxsize 3x1: $k_3$\
      Boxsize 3x3: $k_4$\
      Boxsize 4x1: $k_4$

      Note that the boxsize 3x3 the fourier mode is $k_4$ has the value $0.20$. So this could be inconsitent with the other results for the Xx1
      boxsizes

      </p>
      </details>


    * <details><summary>17.08.2022 &nbsp; Meeting about Boxsize 3x3 & Further Work</summary>
      <p> 

      # Meeting about Boxsize 3x3

      #### Thursday 17.08.2022 14:15 to 14:45 with Florian Rath and Arthur Peeters

      The wavelength is not well defined (in german 'scharf') so because of the results of boxsize 3x3 with the fourier mode $k_4$ the actually result is
      not cruial because the question of the thesis is if the wavelength does converge at all so if $k_4$ or $k_3$ is the stabilizing fourier mode
      is not from intrested. However it would be consistent with the results of Xx1 if 3x3 have had the mode $k_4$. 

      Although the runs are very long the result that the stairscase structure fully developes is remarkable so the underlying process of 
      turbulence that gets stablize through zonal flows holds for even longer runs and hint to a mechanism of toridial plasma.

      The results are very good because boxsize was chosen really well that gets confirmed with the formation of the staircase structure and
      the convergence of the wavelength.

      # Further Work

      To determined the results even more Florian suggests to run more nearby the finite heatflux threshold. To recall all simulations were run with an 
      gradient length $R/L_T$ of $6.0$ because of the run at $6.3$ was not stabilize quick enough as $6.0$. The last simulations sould bei run with an
      gradient length of $6.2$ and a boxsize of 2x2.

      </p>
      </details> 

    * <details><summary>18.08.2022 - 24.08.2022 &nbsp; Result for $R/L_T = 6.2$ & Plots for Thesis </summary>
      <p>

      # Result for $R/L_T = 6.2$

      The simulation for $R/L_T = 6.2$ does stabilize very quick like 1x1, 2x2 and 3x3 so the wavelength does indeed converge with the boxsize. As an 
      view in the future maybe long simulations can even converge even faster if the boxsize in radially and binormal is suitable chosen.

      # Plots for Thesis

      When displaying plots of the wavelength over a time interval it is enough to only display the intresting time intervals e.g. instabil, semi-stabil
      and stabil in addition to that show time intervals where it seems the turbulence is stabilized but the staircase structure is not fully 
      developed.

      </p>
      </details>

    * <details><summary>25.08.2022 - 31.08.2022 &nbsp; Saving data to NAS tp5-peeters </summary>
      <p>

      # Saving data to NAS tp5-peeters

      Cloned repository and copied data to NAS of tp5-peeters. For more informations read [README-DATA](/data/README.md).

      </p>
      </details>

  * <details><summary>September</summary>
    <p>      

    * <details><summary>01.09.2022 &nbsp; Meeting about Publication in Physiccs of Plasma & Layout Bachelor Thesis</summary>
      <p>   

      # Meeting about Publication in Physiccs of Plasma & Layout Bachelor Thesis

      #### Thursday 01.09.2022 14:00 to 14:30 with Florian Rath and Arthur Peeters

      # Layout Bachelor Thesis

      For the bachelor thesis is no other regulation than that of the examination office of the MPI in Bayreuth.

      # Publication in Physiccs of Plasma

      Results of the bachelor thesis should be publicated as brief communication in physics of plasma for that it is cruial to know the
      layout und the rules of the journal.

      </p>
      </details>  

    * <details><summary>02.09.2022 - 30.09.2022 &nbsp; Work on Publication Layout </summary>
      <p>

      # Work on Publication Layout

      Publication should not be longer than $3500$ words and not longer than 4 Pages to ensure that [wordcount.tex](/breifcommunication/wordcount.tex)
      count the words and pages. 

      For Plots only include the necessary plots that are comparison in of the boxsizes in Xx1, XxY, $R/L_T$ and 2x1 between 2x2 and 3x1 between 3x3
      with plots of the turbulence to show simulation is stabil and one example plot with $\omega_{\mathrm{E \times B, max}}$ and
      the corresponding fourier modes.
      No Plots needed for instabil and semi-stabil of wavelength plots in publication only stabil and intresting semi-stabil wavelength plots. 

      </p>
      </details>

    </p>
    </details>

  * <details><summary>October</summary>
    <p>

    * <details><summary>01.10.2022 - 19.10.2022 &nbsp; Work on Publication Plots </summary>
      <p>

      # Work on Publication Plots

      The plots are mostly generated with subplots in matplotlib. To achive a clean look most of work the work gone into programming this plots.
      Especially the plots for the different wavelength took the longest time because it is a overlay off four subplots an idea of myself but not easy 
      to realise. For better visibility the staircase structure of got shifted to achive an overlay of every staircase.

      </p>
      </details>

    * <details><summary>20.10.2022 &nbsp; Meeting about Publication Plots </summary>
      <p>
    
      # Meeting about Publication Plots

      #### Thursday 20.10.2022 14:00 to 14:30 with Arthur Peeters

      # Feedback of Wavelength Plot

      The plot itself is beautiful the only thig that should get added are the linar growth rate $\gamma$ from [[2]](https://doi.org/10.1063/1.4952621).
      The shift of the staircase structure is physically possible because toridial plasma has it symmetry in the rotation so such shift will not 
      affect the nature of the plasma but should certainly be addressed in publication.
      For more space the comparison of gradient length $R/L_T$ will be excluded from publication

      From now on the start of the writing porcess can start. 

      </p>
      </details>
  
    * <details><summary>21.10.2022 - 31.10.2022 &nbsp; Reading Paper [1] & Search References </summary>
      <p>

      # Reading Paper [1] & Search References

      Read paper [[1]](https://doi.org/10.1063/1.4961231) again and extract references from the important section and import them to 
      [references.bib](/briefcommunication/refernces.bib). Briefcommunication will reference paper [[1]](https://doi.org/10.1063/1.4961231) as the basis work.

      </p>
      </details>

    </p>
    </details>

  * <details><summary>November</summary>
    <p>

    * <details><summary>01.11.2022 - 10.11.2022 &nbsp; Reading Paper [...] </summary>
      <p>

      Additional paper found with references
        
      </p>
      </details>
      
    * <details><summary>11.10.2022 &nbsp; Meeting about progess in Publication </summary>
      <p>

      # Meeting about progress in Publication

      #### Friday 11.11.2022 10:15 to 10:45 with Florian Rath

      The shift of wavelength is possible but could get in trouble at the surface of the plasma and at the point the 
      plasma gets connected after one complete round because of the boundary condition. But it will be enough to tell for better 
      visibility the staircase sturcture gets shifted maybe later a explaination would be needed.

      Plan to write publication in the next two weeks.

    * <details><summary>12.11.2022 - 18.11.2022 &nbsp; Restart Script Professonial </summary>
      <p>

      # Restart Script Professonial

      Restart script got overwrite. Now the script creates the jobscript and the status file by itself. Included a parser to give arguments directly
      through the command line and a progressbar with job info that gets updated every 5 seconds. 
      The script can now be run with `nohup` or `screen` ducumentation is included in helpers message.

      </p>
      </details>

    * <details><summary>19.11.2022 - 30.11.2022 &nbsp; Writing Brief Communication and Add Ons in Shell Scripts</summary>
      <p>

      #Writing Brief Communication and Add Ons in Shell Scripts
      ## Writing Brief Communication

      The main focus lay down on making graphics for the brief communication and writing section after section based on Rath2021 und Peeters2016 und Rath2016

      ## Shell Scripts
      
      The copy script got an parser so the use from the command line gets easier and all scripts that needs vpn connection got an updated vpn command for MacOS.

      </p>
      </details>

    </p>
    </details>

  * <details><summary>Dezember</summary>
    <p>

    * <details><summary>06.12.2022 - 29.12.2022 &nbsp; Draft complete Brief Communicationand Issuses with Restart Script</summary>

      # Draft complete for brief communication

      First draft of brif communication was completed (29.12.2022) and After that correction were made on grammar, spelling and graphics. 

      ## Restart Script

      The restart script has multiple prolems:

      * After pending status check the write output to status file stopped
        -> Fix was to only write ones to output file not 12 times in a row like before because the buffer got to fast filled.
      * Parser gets now options for frametype of table, control over sleep time 
      * Change the documentation multiple times 
      * Delete line function rewrite with open to get statusfile correctly closed
      * Send mail function does throw errors again because of whitespaces in subject 
        -> replaced whitspaces with underlines

    </p>
    </details>
    
  </p>
  </details>

* <details><summary>2023</summary>
  <p>

  * <details><summary>January</summary>
    <p>

    * <details><summary>02.01.2023 - 10.01.2023 &nbsp; Corrections Brief Communication and Rerun of box size 4x1</summary>
      <p>

      # Corrections Brief Communication and Rerun of box size 4x1

      ## Correction of Brief Communication

      Thanks to Dominik Müller, Anna-Maria Pleyer and my Sister Cornelia Lippert for reading my first and providing feedback. The corrected version was send to Prof Arthur Peeters and Florian Rath.

      ## Rerun box size 4x1

      A rerun of boxsize 4x1 was made with the goal to get the repetition of the staircase structures aligned with the boxsize. A rerun was necessary because the data file got currupted and can not be fixed.
      The results yield that even after long time intervals of subdued turbulence that the staircase structures got not better aligned with the box size.

      </p>
      </details>

    * <details><summary>16.01.2023 &nbsp; Meetng about Draft of Brief Communcation </summary>
      <p>

      # Meeting about Draft of Brief Communcation

      #### Monday 16.01.2023 14:00 to 14:45 with Florian Rath and Arthur Peeters

      The First draft was good but the focus as an continuation of peeters2016 is an problem. The brief communication should be an paper on its own.

      The case that for 3x3 the staircase structures repeats itself 4 times requires an boxsize scan in binormal direction for that an scan of 3x1.5, 3x2.5 and 3x5 will be made additionally as well a scan for 3x3 for $R/L_T$ = 6.2, 6.4 to lengthen the time of turbulence and to approach the heat flux threshold for verification.

      </p>
      </details>

    * <details><summary>30.01.2023 &nbsp; Presentation Style in LaTeX </summary>
      <p>
      
      # Presentation Style in LaTeX

      To make an presentation about my bachelor work the decision were made in favor of latex because of the cross plattform compability of pdfs.
      As style sheet will be used sleek theme which is an 16:9 variant of HSMR by Benjamin Weiss.

      It could be that the use of powerpoint is necessary to add animations.

      </p>
      </details>
    
    </p>
    </details>

  * <details><summary>February</summary>
    <p>

    * <details><summary>04.02.2023 - 26.02.2023 &nbsp; Corrections of Brief Communication, New Fetaures for Restart Script, Juypter Notebook Problems and Binormal box scan </summary>
      <p>
      
      # Corrections of Brief Communication, New Fetaures for Restart Script, Juypter Notebook Problems and Binormal box scan
      ## Corrections of Brief Communication

      The brief communication has get its first correction from Florian Rath which got accepted and minor things were changed afterwards. 
      Graphics got reworked for new variable names and the box size plot reworked for an big box size plot including radial, isotropic and binormal box size plot.

      ## New features restart script

      The restart script gets new features:
      * Kill option for nohup process to not kill of the false process with multiple user using the script
      * Script continues writing to status file and does not rewrite it
      * New backup locations to chose from praser
      * Increased refresh rate to 300
      * Additional check level to ensure successful run
      * Reset Simulation with dump files. Thanks for Florian Rath to provide the function
      * New jobStatusInfo header to get even informations when slumr `squeue` has no output
      * Script can now be run from everythere so no need to copy it every time again into simultaion folder
      * Check if h5 file is closed before restart and if `FDS.dat` and `gkwdata.h5` has same modified timestamp

      ## Juypter Notebook Problems

      After an update of python and visual studio code the juypter notebooks stop working and lost connection to the server after a image was produced. Because of that the evaluation got rewritten in python file in a new folder to be certian that the evaluation can go on.

      Additional an datasheet was created in `csv` to have an main file with all informations.

      ## Evaluation

      Every simulation converges except for $R/L_T$ = 6.4 which was anticipate. After multiple errors 3x1.5 has to rerun the old file is under the folder Broken.

      Results:
      * 3x1.5 -> Convergence $k_4$
      * 3x2.5 -> Convergence $k_3$, $k_4$
      * 3x5 -> Convergence $k_4$

      * rlt = 6.2 -> Convergence $k_3$
      * rlt = 6.4 -> Turbulent -> Consistent with results of peeters2016

      </p>
      </details>

    </p>
    </details>

  * <details><summary>March</summary>
    <p>

    * <details><summary>19.03.2023 - 22.03.2023 &nbsp; Writing new version of brief communication </summary>
      </p>
      
      # Writing new version of brief communication

      The breif commnication gets additional section for binormal box scan and minor issues were corrected (_ref, rearrangments, colors in plots, captions, name of quantities and formula and commands for quantites that repeats very often).

      </p>
      </details>

    * <details><summary>23.03.2023 - 26.03.2023 &nbsp; New features for Restart script get implemented </summary>
      </p>

      # New features for Restart script get implemented

      As stated befor the restart script gets new features that got implemented into the script itself.
      Additional to that the reset function needs `h5py`, `pandas`and `numpy` modules installed. For that, a function to automatically install modules were made. 

      The code itself got rewritten at some parts because of obsolete code and other isssues were fixed look into [#24](https://github.com/ManeLippert/Bachelorthesis-Shearingrate-Convergence/issues/24) for more.

      </p>
      </details>

    </p>
    </details>

  * <details><summary>April</summary>
    <p>

    * <details><summary>04.04.2023 - 06.04.2023 &nbsp; Submit Brief Communication to "Physics of Plasma" </summary>
      </p>

      # Submit Brief Communication to "Physics of Plasma"

      Brief communication got the last corrections from Arthur Peeters (wavelength -> size/radial size) and Florian Rath and submitted to AIPs "Physics of Plasma"
      
      </p>
      </details>

    * <details><summary>09.04.2023 - 09.05.2023 &nbsp; Writing Bachelor Thesis and Backup Data </summary>
      </p>

      # Writing Bachelor Thesis and Backup Data

      ## Writing Bachelor Thesis

      This time the focus lay down to write the bachelor thesis. For that, changes for brief communication plots were made to ensure readablity and the chapter pages style were changed to a new modern style (cover page will follow up).

      Additionally the error indexs for the 6.2/3x3 simulation were found and written into datasheet. 

      ## Backup Data

      Data got uploaded to NAS and the git repository syncronized as well the server folder gets an clean up.
      
      </p>
      </details>

    </p>
    </details>

  * <details><summary>May</summary>
    <p>

    * <details><summary>17.05.2023 &nbsp; Meeting about Brief Communication Review </summary>
      </p>

      # Meeting about Brief Communication Review

      #### Wednesday 17.05.2023 9:00 to 9:45 with Florian Rath

      Feedback of Referees was good the brief communication will be corrected accordingly for that additional simulation will be performed:

      - Two simulations with box size ```1.5x1.5``` and ```2.5x2.5```
      - Two simulations with different initial conditions ```noise``` and and ```cosine5``` (default = ```cosine2```) for box size ```3x3```
      - Additionally the diagnostics ```xy_kyzero_dens```, ```xy_kyzero_ene_par``` and ```xy_kyzero_ene_perp``` to investigate the influence of the pressure gradient on the shearing rate 

      </p>
      </details>

    </p>
    </details>

  * <details><summary>June</summary>
    <p>

    * <details><summary>07.06.2023 &nbsp; Submission of revised Brief Communication </summary>
      </p>

      # Submission of revised Brief Communication

      The revised brief communictaion got submitted to Physics of Plasma without issues with the corresponding response to the referees

      </p>
      </details>

    * <details><summary>08.06.2023 - 30.06.2023 &nbsp; Correction Bachelor Thesis, Publish Repository and Presentation </summary>
      </p>

      # Correction Bachelor Thesis, Publish Repository and Presentation

      The rest of the month the last finishing touches on my Bachelor Thesis were made and the repository published as well as my presentation prepared. 

      The last data backup were made as well. This is probably the last entrance of this journal.

      </p>
      </details>

  </p>
  </details>

</p>
</details>

## Literature

<details><summary>Literature</summary>
</p>

[1] 2018 nohup. URL https://wiki.ubuntuusers.de/nohup/ – Accessed: 2023-04-15.

[2] 2021 Screen. URL https://wiki.ubuntuusers.de/Screen/ – Accessed: 2023-04-15.

[3] Barton, Justin E., Wehner, William P., Schuster, Eugenio, Felici, Federico & Sauter, Olivier 2015 Simultaneous closed-loop control of the current profile and the electron temperature profile in the tcv tokamak.

[4] Beer, M.A. 1994 Gyrofluid models of turbulent transport in tokamaks. PhD thesis, Princeton University.

[5] Biglari, H., Diamond, P. H. & Terry, P. W. 1990 Influence of sheared poloidal rotation on edge turbulence. Phys. Fluids B: Plasma Physics 2 (1), 1–4.

[6] Brizard, A. J. & Hahm, T. S. 2007 Foundations of nonlinear gyrokinetic theory. Rev. Mod. Phys. 79, 421–468.

[7] Burrell, K. H. 1997 Effects of E×B velocity shear and magnetic shear on tur-
bulence and transport in magnetic confinement devices. Physics of Plasmas 4 (5),1499–1518.

[8] Cary, John R. 1981 Lie transform perturbation theory for Hamiltonian systems.Physics Reports 79 (2), 129–159.

[9] Cary, John R & Littlejohn, Robert G 1983 Noncanonical Hamiltonian mechanics and its application to magnetic field line flow. Annals of Physics 151 (1), 1–34.

[10] Casson, F.J. 2011 Turbulent transport in rotating tokamak plasmas. PhD thesis, University of Warwick.

[11] Coppi, B., Rosenbluth, M. N. & Sagdeev, R. Z. 1967 Instabilities due to temperature gradients in complex magnetic field configurations. The Physics of Fluids 10 (3), 582–587.

[12] Cowley, S. C., Kulsrud, R. M. & Sudan, R. 1991 Considerations of ion‐temperature‐gradient‐driven turbulence. Physics of Fluids B: Plasma Physics 3 (10), 2767–2782.

[13] Dannert, T. 2005 Gyrokinetische Simulation von Plasmaturbulenz mit gefangenen Teilchen und Elektromagnetischen Effekten. PhD thesis, Technische Universtät München.

[14] Diamond, P. H., Itoh, S.-I., Itoh, K. & Hahm, T. S. 2005 Zonal flows in plasma—a review. Plasma Phys. Controlled Fusion 47, R35.

[15] Diamond, P. H. & Kim, Y.‐B. 1991 Theory of mean poloidal flow generation by turbulence. Physics of Fluids B: Plasma Physics 3 (7), 1626–1633.

[16] Dif-Pradalier, G., Diamond, P. H., Grandgirard, V., Sarazin, Y., Abiteboul, J., Garbet, X., Ghendrih, Ph., Strugarek, A., Ku, S. & Chang, C. S. 2010 On the validity of the local diffusive paradigm in turbulent plasma transport. Phys. Rev. E 82, 025401.

[17] Dif-Pradalier, G., Hornung, G., Ghendrih, Ph., Sarazin, Y., Clairet, F., Vermare, L., Diamond, P. H., Abiteboul, J., Cartier-Michaud, T., Ehrlacher, C., Estève, D., Garbet, X., Grandgirard, V., Gürcan, Ö. D., Hennequin, P., Kosuga, Y., Latu, G., Maget, P., Morel, P., Norscini, C., Sabot, R. & Storelli, A. 2015 Finding the elusive E×B staircase in magnetized plasmas. Phys. Rev. Lett. 114, 085004.

[18] Dimits, A. M., Bateman, G., Beer, M. A., Cohen, B. I., Dorland, W., Hammett, G. W., Kim, C., Kinsey, J. E., Kotschenreuther, M., Kritz, A. H., Lao, L. L., Mandrekas, J., Nevins, W. M., Parker, S. E., Redd, A. J., Shumaker, D. E., Sydora, R. & Weiland, J. 2000 Comparisons and physics basis of tokamak transport models and turbulence simulations. Phys. of Plasmas 7 (3), 969–983.

[19] Dubin, Daniel H. E., Krommes, John A., Oberman, C. & Lee, W. W. 1983 Nonlinear gyrokinetic equations. The Physics of Fluids 26 (12), 3524–3535.

[20] Garbet, X., Idomura, Y., Villard, L. & Watanabe, T. H. 2010 Gyrokinetic simulations of turbulent transport. Nuclear Fusion 50 (4).

[21] Hahm, T. S. & Burrell, K. H. 1995 Flow shear induced fluctuation suppression in finite aspect ratio shaped tokamak plasma. Physics of Plasmas 2 (5), 1648–1651.

[22] Hamada, S. 1958 Kakuyugo Kenkyu 1, 542.

[23] Hammett, Greg 2009 The Ion Temperature Gradient (ITG) Instability. CM- PD/CMSO Winter School, UCLA, 1/09/2009.

[24] Hasegawa, Akira & Mima, Kunioki 1978 Pseudo‐three‐dimensional turbulence in magnetized nonuniform plasma. The Physics of Fluids 21 (1), 87–92.

[25] H.Isliker, Pisokas, Th., Strintzi, D. & Vlahos, L. 2010 A self-organized criticality model for ion temperature gradient mode driven turbulence in confined plasma. Physics of Plasmas 17.

[26] Horton, W. 1999 Drift waves and transport. Rev. Mod. Phys. 71, 735–778.

[27] Idomura, Y., Urano, H., Aiba, N. & Tokuda, S. 2009 Study of ion turbulent transport and profile formations using global gyrokinetic full- f vlasov simulation. Nuclear Fusion 49 (6), 065029.

[28] Kim, Y. J., Imadera, K., Kishimoto, Y. & Hahm, T. S. 2022 Transport events and E×B staircase in flux-driven gyrokinetic simulation of ion temperature gradient turbulence. Journal of the Korean Physical Society 81, 636.

[29] Kishimoto, Y., Imadera, K., Ishizawa, A., Wang, W. & Li, J. Q. 2023 Characteristics of constrained turbulent transport in flux-driven toroidal plasmas. Philosophical Transactions of the Royal Society A: Mathematical, Physical and Engineering Sciences 381 (2242), 20210231.

[30] Kosuga, Y., Diamond, P. H. & Gürcan, Ö. D. 2013 How the propagation of heat-flux modulations triggers e×b flow pattern formation. Phys. Rev. Lett. 110, 105002.

[31] Krommes, John A. 2012 The Gyrokinetic Description of Microturbulence in Magnetized Plasmas. Annual Review of Fluid Mechanics.

[32] Krommes, John A. & Kim, Chang-Bae 2000 Interactions of disparate scales in drift-wave turbulence. Phys. Rev. E 62, 8508–8539.

[33] Lippert, M. 2022 torque_monitor.py. URL https://github.com/ ManeLippert/Bachelorthesis-Shearingrate-Convergence/blob/ main/python/torque_monitor.py – Accessed: 2023-04-14.

[34] Lippert, M. & Rath, F. 2023 slurm_monitor.py. URL https://bitbucket.org/gkw/gkw/src/develop/python/slurm_monitor.py – Accessed: 2023-04-12.

[35] Lippert, M., Rath, F. & Peeters, A. G. 2023 Size convergence of the E×B staircase pattern in flux tube simulations of ion temperature gradient driven turbulence. Physics of Plasmas 7 (3), 969–983.

[36] Maeyama, S., Ishizawa, A., Watanabe, T.-H., Nakata, M., Miyato, N., Yagi, M. & Idomura, Y. 2014 Comparison between kinetic-ballooning-mode-driven turbulence and ion-temperature-gradient-driven turbulence. Physics of Plasmas 21 (5), 052301.

[37] Makwana, K. D., Terry, P. W., Pueschel, M. J. & Hatch, D. R. 2014 Subdominant modes in zonal-flow-regulated turbulence. Phys. Rev. Lett. 112, 095002.

[38] McMillan, B. F., Jolliet, S., Tran, T. M., Villard, L., Bottino, A. & Angelino, P. 2009 Avalanchelike bursts in global gyrokinetic simulations. Physics of Plasmas 16 (2), 022310.

[39] Mittendorf, J., Schobert, B. & Müller, D. 2023 Rmhd-code. URL https://bitbucket.org/astro_bayreuth/rmhdcode – Accessed: 2023-04-14.

[40] Müller, D. 2023 Numerical simulations of exor events in protoplanetary disks: Numerical stability and growth of ring structures in the surface density. Bachelorthesis, University of Bayreuth.

[41] Nakata, M., Watanabe, T.-H. & Sugama, H. 2012 Nonlinear entropy transfer via zonal flows in gyrokinetic plasma turbulence. Physics of Plasmas 19, 022303.

[42] Peeters, A. G., Camenen, Y., Casson, F. J., Hornsby, W. A., Snodin, A. P., Strintzi, D. & Szepesi, G. 2009 The nonlinear gyro-kinetic flux tube code gkw. Comput. Phys. Commun. 180, 2650.

[43] Peeters, A. G., Rath, F., Buchholz, R., Camenen, Y., Candy, J., Casson, F. J., Grosshauser, S. R., Hornsby, W. A., Strintzi, D. & Weikl, A. 2016 Gradient-driven flux-tube simulations of ion temperature gradient turbulence close to the non-linear threshold. Physics of Plasmas 23 (8), 082517.

[44] Pueschel, M. J., Kammerer, M. & Jenko, F. 2008 Gyrokinetic turbulence simulations at high plasma beta. Physics of Plasmas 15 (10), 102310.

[45] Rath, F., Peeters, A. G., Buchholz, R., Grosshauser, S. R., Migliano, P., Weikl, A. & Strintzi, D. 2016 Comparison of gradient and flux driven gyro-kinetic turbulent transport. Physics of Plasmas 23 (5), 052309.

[46] Rath, F., Peeters, A. G. & Weikl, A. 2021 Analysis of zonal flow pattern formation and the modification of staircase states by electron dynamics in gyrokinetic near marginal turbulence. Physics of Plasmas 28 (7), 072305.

[47] Rudakov, L.I. & Sagdeev, R.Z. 1961 On the instability of a nonuniform rarefied plasma in a strong magnetic field. Dokl. Akad. Nauk. SSSR 138 (3), 581–583.

[48] Schelter, Dr.rer.nat. Ingo 2016 btrzx2 (2016). URL https://www.bzhpc.uni-bayreuth.de/de/keylab/Cluster/btrzx2_page/index.html – Accessed: 2023-04-14.

[49] Schelter, Dr.rer.nat. Ingo 2020 btrzx1 (2020). URL https://www.bzhpc.uni-bayreuth.de/de/keylab/Cluster/btrzx1_page/index.html – Accessed: 2023-04-12.

[50] Seiferling, F., Peeters, A. G., Grosshauser, S. R., Rath, F. & Weikl, A. 2019 The interplay of an external torque and e×b structure formation in tokamak plasmas. Physics of Plasmas 26 (10), 102306.

[51] Seo, Janghoon, Jhang, Hogun & Kwon, Jae-Min 2022 Effects of light impurities on zonal flow activities and turbulent thermal transport. Physics of Plasmas 29 (5), 052502.

[52] Stroth, Ulrich 2011 Plasmaphysik. Wiesbaden: Viewg+Teubner.

[53] Villard, L, Angelino, P, Bottino, A, Brunner, S, Jolliet, S, McMillan, B F, Tran, T M & Vernay, T 2013 Global gyrokinetic ion temperature gradient turbulence simulations of iter. Plasma Physics and Controlled Fusion 55 (7), 074017.

[54] Waltz, R. E., Dewar, R. L. & Garbet, X. 1998 Theory and simulation of rotational shear stabilization of turbulence. Physics of Plasmas 5 (5), 1784–1792.

[55] Waltz, R. E., Kerbel, G. D. & Milovich, J. 1994 Toroidal gyro-landau fluid modelturbulence simulations in a nonlinearballooning mode representation with radial modes. Physics of Plasmas 1, 2229.

[56] Wang, W., Kishimoto, Y., Imadera, K., Liu, H.R., Li, J.Q., Yagi, M. & Wang, Z.X. 2020 Statistical study for itg turbulent transport in flux-driven tokamak plasmas based on global gyro-kinetic simulation. Nuclear Fusion 60 (6), 066010.

[57] Weikl, A., Peeters, A. G., Rath, F., Grosshauser, S. R., Buchholz, R., Hornsby, W. A., Seiferling, F. & Strintzi, D. 2017 Ion temperature gradient turbulence close to the finite heat flux threshold. Physics of Plasmas 24 (10), 102317.

[58] Wesson, John 2011 Tokamaks. Oxford: Oxford University Press.

[59] Whelan, G. G., Pueschel, M. J. & Terry, P. W. 2018 Nonlinear electromagnetic stabilization of plasma microturbulence. Phys. Rev. Lett. 120, 175002.

[60] Whelan, G. G., Pueschel, M. J., Terry, P. W., Citrin, J., McKinney, I. J., Guttenfelder, W. & Doerk, H. 2019 Saturation and nonlinear electromagnetic stabilization of itg turbulence. Physics of Plasmas 26 (8), 082302.

[61] W.M.Newins, J.Candy, S.Cowley, T.Dannert, A.Dimits, W.Dorland, C.Estrada-Mila, G.W.Hammet, F.Jenko, M.J.Pueschel & D.E.Shumaker 2006 Characterizing electron temperature gradient turbulence via numerical simulations. Physics of Plasmas 13.

</p>
</details>