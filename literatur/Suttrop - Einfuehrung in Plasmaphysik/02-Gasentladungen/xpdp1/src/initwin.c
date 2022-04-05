#include "pdp1.h"
#include "xgrafix.h"
#include <string.h>

void XGWrite(void *, int, int, FILE *, char *); 
void XGRead(void *, int, int, FILE *, char *);


/****************************************************************/

void InitWindows(void)
{
  register int i, isp;
  char buffer[50];

  /*********************************************/
  /* Set up each window structure              */
  
  if(nsp) {
    if(nsp > 1)
      {
	XGSet2D("linlin", "Time", "s1_dens(t)", "closed", 100, 400, 1.0, 1.0,
		True, True, 0.0, 0.0, 0.0, 0.0);
	XGCurve(t_array, s1n_hist, &hist_hi, 4);
	XGSet2D("linlin", "Time", "s2_dens(t)", "closed", 100, 400, 1.0, 1.0,
		True, True, 0.0, 0.0, 0.0, 0.0);
	XGCurve(t_array, s2n_hist, &hist_hi, 4);
	XGSet2D("linlin", "Time", "s1+s2_dens(t)", "closed", 100, 400, 1.0, 1.0,
		True, True, 0.0, 0.0, 0.0, 0.0);
	XGCurve(t_array, s1s2n_hist, &hist_hi, 4);
      }

    XGSet2D("linlin", "Time", "Number(t)", "closed", 100, 400, 1.0, 1.0,
						True, True, 0.0, 0.0, 0.0, 0.0);
    for (isp=0; isp<nsp; isp++) XGCurve(t_array, np_hist[isp], &hist_hi, isp);

    XGSet2D("linlin", "Time", "Bulk Velocity(t)", "closed", 100, 400, 1.0, 1.0,
						True, True, 0.0, 0.0, 0.0, 0.0);
    for (isp=0; isp<nsp; isp++) XGCurve(t_array, vel_hist[isp], &hist_hi, isp);    
    
    /*********************************************/

/*		XGSet2D("linlin", "Time", "Trapped(t)", "closed", 100, 400, 1.0, 1.0,
						True, True, 0.0, 0.0, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++) XGCurve(t_array, np_trapped[isp], &hist_hi, isp);
    */
    /*********************************************/
		/*
		XGSet2D("linlin", "Time", "UnTrapped(t)", "closed", 100, 400, 1.0, 1.0,
						True, True, 0.0, 0.0, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++) XGCurve(t_array, np_untrapped[isp], &hist_hi, isp);
    */
    /*********************************************/   
 /*
		XGSet2D("linlin", "Time", "UnTrapped energy (t) \\[eV\\]", "closed", 100, 400, 1.0, 1.0,
						True, True, 0.0, 0.0, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++) XGCurve(t_array, TE_untrapped[isp], &hist_hi, isp);
   */ 
    /*********************************************/
/*
		XGSet2D("linlin", "Time", "Trapped energy (t) \\[eV\\]", "closed", 100, 400, 1.0, 1.0,
						True, True, 0.0, 0.0, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++) XGCurve(t_array, TE_trapped[isp], &hist_hi, isp);
  */  
    /*********************************************/
/*		
		XGSet2D("linlin", "Time", "Total Particle energy (t) \\[eV\\]", "closed", 100, 400, 1.0, 1.0,
						True, True, 0.0, 0.0, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++) XGCurve(t_array, TE_particle[isp], &hist_hi, isp);
  */  
    /*********************************************/
    for (isp=0; isp<nsp; isp++) {
      sprintf(buffer, "Vx-X Phase Space %d", isp+1);
      
      XGSet2D("linlin", "X", strdup(buffer), "closed", 300, 10, dx, dx/dt,
							False, False, 0.0, length, vscale*(v0[isp][1]+5*vt[isp][1]), 
							vscale*(v0[isp][0]+5*vt[isp][0]));
      
      XGScat2D(x[isp], vx[isp], &np[isp], isp);
    }
    /*********************************************/
		
		if (nsp>1){
			sprintf(buffer, "Vx-X Phase Space all species");
			
			XGSet2D("linlin", "X", strdup(buffer), "closed", 300, 10, dx, dx/dt,
							False, False, 0.0, length, vscale*(v0[isp][1]+5*vt[isp][1]), 
							vscale*(v0[isp][0]+5*vt[isp][0]));
			for (isp=0; isp<nsp; isp++) 
				XGScat2D(x[isp], vx[isp], &np[isp], isp);
		}
    
    /*********************************************/
		for (isp=0; isp<nsp; isp++) {
      sprintf(buffer, "Vy-X Phase Space %d", isp+1);
      
      XGSet2D("linlin", "X", strdup(buffer), "closed", 300, 10, dx, dx/dt,
							False, False, 0.0, length, vscale*(v0y[isp]-5*vty[isp]),
							vscale*(v0y[isp]+5*vty[isp]));
      
      XGScat2D(x[isp], vy[isp], &np[isp], isp);
    }
    /*********************************************/
		for (isp=0; isp<nsp; isp++) {
      sprintf(buffer, "Vz-X Phase Space %d", isp+1);
      
      XGSet2D("linlin", "X", strdup(buffer), "closed", 300, 10, dx, dx/dt,
							False, False, 0.0, length, vscale*(v0z[isp]-5*vtz[isp]), 
							vscale*(v0z[isp]+5*vtz[isp]));
      
      XGScat2D(x[isp], vz[isp], &np[isp], isp);
    }
    /*********************************************/
		
    for (isp=0; isp<nsp; isp++) {
			sprintf(buffer, "Vy-Vx Phase Space %d", isp+1);
			XGSet2D("linlin", "Vx", strdup(buffer), "closed", 400, 10, dx/dt, dx/dt,
							False, False, vscale*(v0[isp][1]+5*vt[isp][1]), 
							vscale*(v0[isp][0]+5*vt[isp][0]), vscale*(v0y[isp]-5*vty[isp]),
							vscale*(v0y[isp]+5*vty[isp]));
			
			XGScat2D(vx[isp], vy[isp], &np[isp], isp);
    }
    /*********************************************/
		for (isp=0; isp<nsp; isp++) {
			sprintf(buffer, "Vz-Vx Phase Space %d", isp+1);
			XGSet2D("linlin", "Vx", strdup(buffer), "closed", 400, 10, dx/dt, dx/dt,
							False, False, vscale*(v0[isp][1]+5*vt[isp][1]), 
							vscale*(v0[isp][0]+5*vt[isp][0]), vscale*(v0z[isp]-5*vtz[isp]),
							vscale*(v0z[isp]+5*vtz[isp]));
			
			XGScat2D(vx[isp], vz[isp], &np[isp], isp);
    }
    /*********************************************/
		for (isp=0; isp<nsp; isp++) {
			sprintf(buffer, "Vy-Vz Phase Space %d", isp+1);
			XGSet2D("linlin", "Vz", strdup(buffer), "closed", 400, 10, dx/dt, dx/dt,
							False, False, vscale*(v0z[isp]-5*vtz[isp]),
							vscale*(v0z[isp]+5*vtz[isp]), vscale*(v0y[isp]-5*vty[isp]),
							vscale*(v0y[isp]+5*vty[isp]));
			
			XGScat2D(vz[isp], vy[isp], &np[isp], isp);
    }
    /*********************************************/

    XGSet2D("linlin", "X", "N(x)", "closed", 600, 500, dx, nc2p/area/dx,
						False, True, 0.0, length, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++) XGCurve(x_grid, sp_n[isp], &ng, isp);

    /*********************************************/
    
    XGSet2D("linlin", "X", "Ux(x)", "closed", 650, 50, dx, 0.5*dx/dt,
						False, True, 0.0, length, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++) XGCurve(x_grid, sp_u_x_show[isp], &ng, isp);

    /*********************************************/
		
		XGSet2D("linlin", "X", "Uy(x)", "closed", 650, 50, dx, 0.5*dx/dt,
						False, True, 0.0, length, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++) XGCurve(x_grid, sp_u_y_show[isp], &ng, isp);

    /*********************************************/

		  XGSet2D("linlin", "X", "Uz(x)", "closed", 650, 50, dx, 0.5*dx/dt,
						False, True, 0.0, length, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++) XGCurve(x_grid, sp_u_z_show[isp], &ng, isp);

    /*********************************************/
    XGSet2D("linlin", "X", "KE(x) \\[eV\\]", "closed", 650, 50, dx, 0.25,
						False, True, 0.0, length, 0.0, 0.0);

    for (isp=0; isp<nsp; isp++) XGCurve(x_grid, sp_ke_show[isp], &ng, isp);

    /*********************************************/

    XGSet2D("linlin", "X", "KE(x) in X \\[eV\\]", "closed", 650, 50, dx, 0.25,
						False, True, 0.0, length, 0.0, 0.0);

    for (isp=0; isp<nsp; isp++) XGCurve(x_grid, sp_ke_x_show[isp], &ng, isp);

    /*********************************************/ 

		XGSet2D("linlin", "X", "KE(x) in Y \\[eV\\]", "closed", 650, 50, dx, 0.25,
						False, True, 0.0, length, 0.0, 0.0);

    for (isp=0; isp<nsp; isp++) XGCurve(x_grid, sp_ke_y_show[isp], &ng, isp);

    /*********************************************/ 
		
		XGSet2D("linlin", "X", "KE(x) in Z \\[eV\\]", "closed", 650, 50, dx, 0.25,
						False, True, 0.0, length, 0.0, 0.0);

    for (isp=0; isp<nsp; isp++) XGCurve(x_grid, sp_ke_z_show[isp], &ng, isp);

    /*********************************************/ 
   
    XGSet2D("linlin", "X", "Jx(x)", "closed", 650, 50, dx, 1.0,
						False, True, 0.0, length, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++) XGCurve(x_grid, sp_j_x_show[isp], &ng, isp);

    /*********************************************/

		XGSet2D("linlin", "X", "Jy(x)", "closed", 650, 50, dx, 1.0,
						False, True, 0.0, length, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++) XGCurve(x_grid, sp_j_y_show[isp], &ng, isp);

    /*********************************************/

		XGSet2D("linlin", "X", "Jz(x)", "closed", 650, 50, dx, 1.0,
						False, True, 0.0, length, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++) XGCurve(x_grid, sp_j_z_show[isp], &ng, isp);

    /*********************************************/
    
    XGSet2D("linlin", "X", "J-E(x)", "closed", 650, 50, dx, 1.0,
						False, True, 0.0, length, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++) XGCurve(x_grid, jdote_show[isp], &ng, isp);
  }
  /*********************************************/
  
  XGSet2D("linlin", "X", "Rho(x)", "closed", 400, 100, dx, 1.0,
					False, True, 0.0, length, 0.0, 0.0);
  
  XGCurve(x_grid, rho, &ng, 2);
  
  /*********************************************/
  
  XGSet2D("linlin", "X", "E field(x)", "closed", 500, 10, dx, 1.0,
					True, True, 0.0, length, 0.0, 0.0);
  
  XGCurve(x_grid, e, &ng, 3);
  
  /*********************************************/

  XGSet2D("linlin", "X", "Potential(x)", "closed", 400, 300, dx, 1.0,
					False, True, 0.0, length, 0.0, 0.0);
  
  XGCurve(x_grid, phi, &ng, 4);
  
  /*********************************************/

  XGSet2D("linlin", "Time", "Mid Potential(t)", "closed", 500, 400, 1.0, 1.0,
					True, True, 0.0, 0.0, 0.0, 0.0);
  
  XGCurve(t_array, com_phi_hist[1], &hist_hi, 4);

  /*********************************************/
  
  XGSet2D("linlin", "Time", "LHS Potential(t)", "closed", 700, 400, 1.0, 1.0,
					True, True, 0.0, 0.0, 0.0, 0.0);
  
  XGCurve(t_array, com_phi_hist[0], &hist_hi, 4);

  /*********************************************/

  XGSet2D("linlin", "Time", "Current(t)", "closed", 600, 400, 1.0, 1.0,
					True, True, 0.0, 0.0, 0.0, 0.0);
  
  XGCurve(t_array, com_cur_hist, &hist_hi, 4);

  /*********************************************/

  XGSet2D("linlin", "I(t)", "V(t)", "closed", 600, 400, 1.0, 1.0,
					True, True, 0.0, 0.0, 0.0, 0.0);
  
  XGCurve(com_cur_hist, com_phi_hist[0], &hist_hi, 4);

  /*********************************************/

  XGSet2D("linlin", "Time", "Power(t)", "closed", 550, 450, 1.0, 1.0,
					True, True, 0.0, 0.0, 0.0, 0.0);
  
  XGCurve(t_array, com_pow_hist, &hist_hi, 4);

  /*********************************************/
  
  XGSet2D("linlin", "Time", "s1(t)", "closed", 600, 400, 1.0, 1.0,
	  True, True, 0.0, 0.0, 0.0, 0.0);
  
  XGCurve(t_array, s1_hist, &hist_hi, 4);
  
  /*********************************************/
  
  XGSet2D("linlin", "Time", "s2(t)", "closed", 600, 400, 1.0, 1.0,
	  True, True, 0.0, 0.0, 0.0, 0.0);
  
  XGCurve(t_array, s2_hist, &hist_hi, 4);
  
  /*********************************************/
  
  XGSet2D("linlin", "Time", "(s1 + s2)(t)", "closed", 600, 400, 1.0, 1.0,
	  True, True, 0.0, 0.0, 0.0, 0.0);
  
  XGCurve(t_array, s1s2_hist, &hist_hi, 4);

 /*********************************************/
  
  XGSet2D("linlin", "Time", "s1_E=0_(t)", "closed", 600, 400, 1.0, 1.0,
	  True, True, 0.0, 0.0, 0.0, 0.0);
  
  XGCurve(t_array, s1e0, &hist_hi, 4);

 /*********************************************/
  
  XGSet2D("linlin", "Time", "s2_E=0_(t)", "closed", 600, 400, 1.0, 1.0,
	  True, True, 0.0, 0.0, 0.0, 0.0);
  
  XGCurve(t_array, s2e0, &hist_hi, 4);

 /*********************************************/
  
  XGSet2D("linlin", "Time", "Vc_E=0(t)", "closed", 600, 400, 1.0, 1.0,
	  True, True, 0.0, 0.0, 0.0, 0.0);
  
  XGCurve(t_array, v2e0, &hist_hi, 4);

 /*********************************************/

  XGSet2D("linlin", "Time", "Sigma(t)", "closed", 500, 500, 1.0, 1.0,
					True, True, 0.0, 0.0, 0.0, 0.0);
  
  XGCurve(t_array, wall_sigma_hist, &hist_hi, 3);
  
  /*********************************************/
  
  XGSet2D("linlog", "Time", "Field Energy(t) \\[J\\]", "closed", 650, 400, 1.0, 1.0,
					True, True, 0.0, 0.0, 0.0, 0.0);
  
  XGCurve(t_array, ese_hist, &hist_hi, 4);
  
  /*********************************************/
  
  if(nsp) {
    XGSet2D("linlin", "Time", "Ave Kinetic Energy(t) \\[eV\\]", "closed",
						500, 200, 1.0, 0.25, True, True, 0.0, 0.0, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++)
      XGCurve(t_array, kes_hist[isp], &hist_hi, isp);
    
    /*********************************************/

		XGSet2D("linlin", "Time", "Ave Kinetic Energy(t) in X \\[eV\\]", "closed",
						500, 200, 1.0, 0.25, True, True, 0.0, 0.0, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++)
      XGCurve(t_array, kes_x_hist[isp], &hist_hi, isp);
    
    /*********************************************/

		XGSet2D("linlin", "Time", "Ave Kinetic Energy(t) in Y \\[eV\\]", "closed",
						500, 200, 1.0, 0.25, True, True, 0.0, 0.0, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++)
      XGCurve(t_array, kes_y_hist[isp], &hist_hi, isp);
    
    /*********************************************/

		XGSet2D("linlin", "Time", "Ave Kinetic Energy(t) in Z \\[eV\\]", "closed",
						500, 200, 1.0, 0.25, True, True, 0.0, 0.0, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++)
      XGCurve(t_array, kes_z_hist[isp], &hist_hi, isp);
    
    /*********************************************/

    XGSet2D("linlin", "Time", "Wall Charge(t)", "closed", 600, 200,
						1.0, 1.0, True, True, 0.0, 0.0, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++)
      XGCurve(t_array, jwall_hist[isp], &hist_hi, isp);
    
    /*********************************************/
    
    for(isp=0; isp<nsp; isp++) {
      sprintf(buffer, "F(E) %d", isp+1);
      
      XGSet2D("linlin", "E", strdup(buffer), "closed", 650, 350,
							Escale[isp], nc2p, False, True, e_array[isp][0]*Escale[isp],
							e_array[isp][nbin[isp]-1]*Escale[isp], 0.0, 0.0);
      
      XGCurve(e_array[isp], fe[isp], &nbin[isp], isp);
    }

		/*********************************************/

    for(isp=0; isp<nsp; isp++) {
      sprintf(buffer, "Local F(E) mid %d", isp+1);
      
      XGSet2D("linlog", "E", strdup(buffer), "closed", 650, 350,
							Escale[isp], nc2p, False, True, emin_mid[isp]*Escale[isp],
							e_mid_array[isp][nbin_mid[isp]-1]*Escale[isp], 0.0, 0.0);
      
      XGCurve(e_mid_array[isp], sp_fe_show[isp], &nbin_mid[isp], isp);
    }

		/*********************************************/

		for(isp=0; isp<nsp; isp++) {
      sprintf(buffer, "F(E) Mid %d", isp+1);
      
      XGSet2D("linlog", "E", strdup(buffer), "closed", 650, 350,
							Escale[isp], nc2p, False, True, emin_mid[isp]*Escale[isp],
							e_mid_array[isp][nbin_mid[isp]-1]*Escale[isp], 0.0, 0.0);
      
      XGCurve(e_mid_array[isp], fe_mid_show[isp], &nbin_mid[isp], isp);
    }
    
    /*********************************************/
    
    XGSet2D("linlin", "Theta", "F(theta)", "closed", 600, 350,
						180/M_PI, 1.0, False, True, 0.0, 90.0, 0.0, 0.0);
    
    for (isp=0; isp<nsp; isp++)
      XGCurve(th_array[isp], ftheta[isp], &nbin[isp], isp);
    
    /*********************************************/

		
		for(isp=0; isp<nsp; isp++)
			if (nvxbin[isp]){
				sprintf(buffer, "v_x dist species %d", isp+1);
				
				XGSet3D("linlinlin", "X", "velocity", strdup(buffer), 0.0, 0.0, "closed",
								700, 175, dx, 1, 1.0, True, True, True, 0, 1.0, 0, 1.0, 0, 1.0);
				
				XGSurf(x_grid, vx_array[isp], vx_dist[isp], &ng, &nvxbin[isp], 4); 
			}
		
		/*********************************************/

		for(isp=0; isp<nsp; isp++)
			if (nvybin[isp]){
				sprintf(buffer, "v_y dist species %d", isp+1);
				
				XGSet3D("linlinlin", "X", "velocity", strdup(buffer), 0.0, 0.0, "closed",
								700, 175, dx, 1, 1.0, True, True, True, 0, 1.0, 0, 1.0, 0, 1.0);
				
				XGSurf(x_grid, vy_array[isp], vy_dist[isp], &ng, &nvybin[isp], 4); 
			}
		
		/*********************************************/
		
		for(isp=0; isp<nsp; isp++)
			if (nvzbin[isp]){
				sprintf(buffer, "v_z dist species %d", isp+1);
				
				XGSet3D("linlinlin", "X", "velocity", strdup(buffer), 0.0, 0.0, "closed",
								700, 175, dx, 1, 1.0, True, True, True, 0, 1.0, 0, 1.0, 0, 1.0);
				
				XGSurf(x_grid, vz_array[isp], vz_dist[isp], &ng, &nvzbin[isp], 4); 
			}
		
    /*********************************************/
		
		for(i=0; i<ndiag; i++) {
      XGSet2D("linlin", "X", rate_title[i], "closed",
							100, 500, dx, nc2p/(dx*area*dt), False, True, 0.0, length, 0.0, 0.0);
      
      XGCurve(x_grid, rate_show[i], &ng, 3);
    }
  }
  /*********************************************/
  
  if(n_ave) {
    if(nsp) {
      XGSet2D("linlin", "X", "Time Ave n(x)", "closed",
							600, 500, dx, nc2p/area/dx, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, sp_n_ave_show[isp], &ng, isp);
      
      /*********************************************/
      
      XGSet2D("linlin", "X", "Time Ave Ux(x)", "closed",
							600, 500, dx, 0.5*dx/dt, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, sp_u_x_ave_show[isp], &ng, isp);
      
      /*********************************************/

			  XGSet2D("linlin", "X", "Time Ave Uy(x)", "closed",
							600, 500, dx, 0.5*dx/dt, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, sp_u_y_ave_show[isp], &ng, isp);
      
      /*********************************************/

			  XGSet2D("linlin", "X", "Time Ave Uz(x)", "closed",
							600, 500, dx, 0.5*dx/dt, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, sp_u_z_ave_show[isp], &ng, isp);
      
      /*********************************************/
      
      XGSet2D("linlin", "X", "Time Ave Jx(x)", "closed",
							600, 500, dx, 1.0, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, sp_j_x_ave_show[isp], &ng, isp);
      
      /*********************************************/
      
      XGSet2D("linlin", "X", "Time Ave Jy(x)", "closed",
							600, 500, dx, 1.0, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, sp_j_y_ave_show[isp], &ng, isp);
      
      /*********************************************/
			      
      XGSet2D("linlin", "X", "Time Ave Jz(x)", "closed",
							600, 500, dx, 1.0, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, sp_j_z_ave_show[isp], &ng, isp);
      
      /*********************************************/
      
      XGSet2D("linlin", "X", "Time Ave KE(x)", "closed",
							600, 500, dx, 0.25, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, sp_ke_ave_show[isp], &ng, isp);

			/*********************************************/

			XGSet2D("linlin", "X", "Time Ave KEx(x)", "closed",
							600, 500, dx, 0.25, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, sp_ke_x_ave_show[isp], &ng, isp);

			/*********************************************/

			XGSet2D("linlin", "X", "Time Ave KEy(x)", "closed",
							600, 500, dx, 0.25, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, sp_ke_y_ave_show[isp], &ng, isp);

			/*********************************************/

			XGSet2D("linlin", "X", "Time Ave KEz(x)", "closed",
							600, 500, dx, 0.25, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, sp_ke_z_ave_show[isp], &ng, isp);

			/*********************************************/
      
      XGSet2D("linlin", "X", "Time Ave T(x)", "closed",
							600, 500, dx, 0.5, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, T_ave_show[isp], &ng, isp);    

			/*********************************************/
      
      XGSet2D("linlin", "X", "Time Ave Tx(x)", "closed",
							600, 500, dx, 0.5, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, Tx_ave_show[isp], &ng, isp);    

			/*********************************************/
      
      XGSet2D("linlin", "X", "Time Ave Ty(x)", "closed",
							600, 500, dx, 0.5, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, Ty_ave_show[isp], &ng, isp);    

			/*********************************************/
      
      XGSet2D("linlin", "X", "Time Ave Tz(x)", "closed",
							600, 500, dx, 0.5, False, True, 0.0, length, 0.0, 0.0);
      
      for (isp=0; isp<nsp; isp++)
				XGCurve(x_grid, Tz_ave_show[isp], &ng, isp);    
		}
		/*********************************************/
    
    XGSet2D("linlin", "X", "Time Ave E field(x)", "closed",
						10, 500, dx, 1.0, False, True, 0.0, length, 0.0, 0.0);
    
    XGCurve(x_grid, e_ave_show, &ng, 3);
    
    /*********************************************/
    
    XGSet2D("linlin", "X", "Time Ave Potential(x)", "closed",
						400, 300, dx, 1.0, False, True, 0.0, length, 0.0, 0.0);
    
    XGCurve(x_grid, phi_ave_show, &ng, 4);
    
    /*********************************************/
	}

	if (nfft){
		/*********************************************/
    XGSet2D("linlin", "Time", "Local Mid Potential(t)", "closed",
						300, 550, 1.0, 1.0, True, True, 0.0, 0.0, 0.0, 0.0);
    
    XGCurve(Local_t_array, phi_hist[1], &thist_hi, 3);
    
    /*********************************************/
    
    XGSet2D("linlin", "Time", "Local LHS Potential(t)", "closed",
						300, 550, 1.0, 1.0, True, True, 0.0, 0.0, 0.0, 0.0);
    
    XGCurve(Local_t_array, phi_hist[0], &thist_hi, 3);
    
    /*********************************************/
    
    XGSet2D("linlin", "Time", "Local Current(t)", "closed",
						300, 550, 1.0, 1.0, True, True, 0.0, 0.0, 0.0, 0.0);
    
    XGCurve(Local_t_array, cur_hist, &thist_hi, 3);
    
    /*********************************************/
    
    XGSet2D("linlin", "Time", "Local Power(t)", "closed",
						300, 550, 1.0, 1.0, True, True, 0.0, 0.0, 0.0, 0.0);
    
    XGCurve(Local_t_array, pow_hist, &thist_hi, 3);

    /*********************************************/
    
    XGSet2D("linlog", "Freq", "Mag of Mid V(f)", "closed",
						100, 250, 1.0, 1.0, True, True, 0.0, 0.0, 0.0, 0.0);
    
    XGCurve(f_array, mphi_fft, &freq_hi, 4);

    /*********************************************/

    XGSet2D("linlin", "Freq", "Phase of Mid V(f)", "closed",
						400, 250, 1.0, 1.0, True, False, 1.0, 1.0, -180.0, 180.0);
    
    XGCurve(f_array, mphi_fft+freq_hi, &freq_hi, 4);
    
    /*********************************************/

    XGSet2D("linlog", "Freq", "Mag of LHS V(f)", "closed",
						100, 250, 1.0, 1.0, True, True, 0.0, 0.0, 0.0, 0.0);
    
    XGCurve(f_array, phi_fft, &freq_hi, 4);

    /*********************************************/

    XGSet2D("linlin", "Freq", "Phase of LHS V(f)", "closed",
						400, 250, 1.0, 1.0, True, False, 1.0, 1.0, -180.0, 180.0);
    
    XGCurve(f_array, phi_fft+freq_hi, &freq_hi, 4);
    
    /*********************************************/

    XGSet2D("linlog", "Freq", "Mag of I(f)", "closed",
						100, 250, 1.0, 1.0, True, True, 0.0, 0.0, 0.0, 0.0);
    
    XGCurve(f_array, cur_fft, &freq_hi, 4);

    /*********************************************/

    XGSet2D("linlin", "Freq", "Phase of I(f)", "closed",
						400, 250, 1.0, 1.0, True, False, 1.0, 1.0, -180.0, 180.0);
    
    XGCurve(f_array, cur_fft+freq_hi, &freq_hi, 4);
    
    /*********************************************/

    XGSet2D("linlog", "Freq", "Mag of Power(f)", "closed",
						100, 250, 1.0, 1.0, True, True, 0.0, 0.0, 0.0, 0.0);
    
    XGCurve(f_array, pow_fft, &freq_hi, 4);

    /*********************************************/

    XGSet2D("linlin", "Freq", "Phase of Power(f)", "closed",
						400, 250, 1.0, 1.0, True, False, 1.0, 1.0, -180.0, 180.0);
    
    XGCurve(f_array, pow_fft+freq_hi, &freq_hi, 4);
  }
}

/***************************************************************/
/* Dumping the current state of the system into a binary file. */
/* Note: the binary file is written in the format used by IBM  */
/* and DEC compatible machines (low bytes followed by high     */
/* bytes). For this Dump() and Restore() call the functions    */
/* sun_read() and sun_write().                                 */

char Revision[]={'1','.','1','1'};

void Dump(char *filename)
{
  register int i, isp;
  SCALAR ftemp;
  FILE *DMPFile;
  
  if ((DMPFile = fopen(filename, "w")) == NULL) {
    puts("Dump: open failed");
    return;
  }

  XGWrite(Revision, 1, 4, DMPFile, "char");
  ftemp = t;
  XGWrite(&ftemp, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
  XGWrite(&oldsigma, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
  XGWrite(&extq, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
  XGWrite(&extq_1, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
  XGWrite(&extq_2, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
  XGWrite(&extq_3, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
  XGWrite(&nsp, sizeof(SCALAR), 1, DMPFile, "int");
  XGWrite(jwall, sizeof(SCALAR), nsp, DMPFile, SCALAR_CHAR);
  XGWrite(np, sizeof(SCALAR), nsp, DMPFile, "int");
  
  for (isp=0; isp<nsp; isp++) {
    for (i=0; i<np[isp]; i++) {
      ftemp = x[isp][i]/xnc;
      XGWrite(&ftemp, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
      ftemp = vscale*vx[isp][i];
      XGWrite(&ftemp, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
      ftemp = vscale*vy[isp][i];
      XGWrite(&ftemp, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
      ftemp = vscale*vz[isp][i];
      XGWrite(&ftemp, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
    }
  }
  fclose(DMPFile);
}

/****************************************************************/

void Restore(char *filename)
{
  char Rev[5];
  register int i, isp;
  SCALAR ftemp;
  FILE *DMPFile;

  if ((DMPFile = fopen(filename, "r+b")) == NULL) {
    puts("Dump: open failed");
    return;
  }
  XGRead(Rev, 1, 4, DMPFile, "char");
  for (i=0; i<4; i++)
    if (Rev[i]!=Revision[i]) {
      puts("Incompatible dump file version");
      exit(1);
    }
  XGRead(&ftemp, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
  t = ftemp;
  XGRead(&sigma, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
  XGRead(&extq, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
  XGRead(&extq_1, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
  XGRead(&extq_2, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
  XGRead(&extq_3, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
  XGRead(&nsp, sizeof(SCALAR), 1, DMPFile, "int");
  XGRead(jwall, sizeof(SCALAR), nsp, DMPFile, SCALAR_CHAR);
  XGRead(np, sizeof(SCALAR), nsp, DMPFile, "int");

  for (isp=0; isp<nsp; isp++) {
    for (i=0; i<np[isp]; i++) {
      XGRead(&ftemp, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
      x[isp][i] = ftemp*xnc;
      XGRead(&ftemp, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
      vx[isp][i]= ftemp/vscale;
      XGRead(&ftemp, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
      vy[isp][i]= ftemp/vscale;
      XGRead(&ftemp, sizeof(SCALAR), 1, DMPFile, SCALAR_CHAR);
      vz[isp][i]= ftemp/vscale;
    }
  }
  fclose(DMPFile);
}

/****************************************************************/

void Quit(void)
{
  
}

/****************************************************************/



