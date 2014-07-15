/******************************************************************************
  Module Name : bulk.c
  Module Date : 04/24/2014
  Module Auth : Yonggang Li

  Description : Bulk or multi-layer samples.

  Others :
      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#include "bulk.h"

/*=============================================================================
  Function Name : read_csg_shape
  Description   : Read bulk shape parameters form Config file.

  Inputs  : char* file_name
  Outputs : no.

  Notes : no.
=============================================================================*/
int read_bulk_shape (char *file_name) {
	int i;

    shape_file = fopen (file_name, "rt");  /* TargetCompositionFileName */
    if (shape_file == NULL) {
        printf ("Error: File '%s' cannot be opend for reading.\n",
                 TargetCompositionFileName);
    	return -1000;
    }

    fscanf (shape_file, "%3i", &max_no_layers);
#ifdef MPI_PRALLEL
    /* MPI=============================================== */
    if (my_node == ROOT) printf ("\nMax number of layers: %i\n", max_no_layers);
    /* MPI=============================================== */
#else
    printf ("\nMax number of layers: %i\n", max_no_layers);
#endif
    if (max_no_layers >= 20) {
	    printf ("Error: Number of layers (%2i) in bulk materials is no more than 20!",
                 max_no_layers);
        return -1100;
    }

    sect[0] = 0.0;  /* vacuum */
    if (no_substrate) {  /* thin multi-layers without substrate */
        sect[max_no_layers+1] = 0.0;
    }
    else {  /* 1/l0 of the half-infinite substrate */
        sect[max_no_layers+1] = list_of_materials[0].rev_atomic_distance;
    }
    sect0 = list_of_materials[0].rev_atomic_distance;
    material_of_layer[0] = 0;  /* vacuum */
    material_of_layer[max_no_layers+1] = 1;  /* half-infinite substrate */
    thick_of_layer[0] = sub_surf_z;
    for (i=1; i<=max_no_layers; i++) {
    	fscanf (shape_file, "%i",  &material_of_layer[i]);
        fscanf (shape_file, "%lf", &thick_of_layer[i]);

        if (material_of_layer[i] == 0) {
        	sect[i] = 0.0;
        }
        else {
            sect[i] = list_of_materials[material_of_layer[i]-1].rev_atomic_distance;
        }

#ifdef MPI_PRALLEL
        /* MPI=============================================== */
        if (my_node == ROOT) {
            printf ("Layer %i:\n", i);
            printf ("    Material  : %i\n", material_of_layer[i]);
            printf ("    Thick (nm): %g\n", thick_of_layer[i]);
        }
        /* MPI=============================================== */
#else
        printf ("Layer %i:\n", i);
        printf ("    Material  : %i\n", material_of_layer[i]);
        printf ("    Thick (nm): %g\n", thick_of_layer[i]);
#endif
        thick_of_layer[i] += thick_of_layer[i-1];  /* distance from surface */
    }

    fclose (shape_file);

    return 0;
}

/*=============================================================================
  Function Name : freepath_bulk
  Description   : Free flight path of ion for bulk or multi-layer samples.

  Inputs  : int proj_Z, int proj_M
            float *energy
            double z, double vz
            double flight_length
            int *current_material_index
            int *ion_left_target
  Outputs : float *energy
            double *s
            int *current_material_index
            int *ion_left_target
            float *stopping
            float *straggling

  Notes :
      Function call:
        fetm - float check_surf_bind (int current_material_index, int projZ);
        bulk - double dist_to_surf (double z, double z0, double vz);
        transport - void e_energy_loss (int projZ, int projM, float projE,
                        int mater_type, double s, float *stopping,
                        float *straggling);
=============================================================================*/
void freepath_bulk (int projZ, int projM, float *energy, double z, double vz,
                    double flight_length, double *s, int *current_material_index,
                    int *ion_left_target, float *stopping, float *straggling) {
    int    i, j;
    float  old_surf_bind, current_surf_bind;
    double s_const = 0.0, s_temp, si = 0.0;
    double T0, TT[20];

    if (flight_length_type == 1) {  /* 1 - constant */
        s_const = flight_length;
    }
    else {  /* 0 or default - random */
        s_const = log_list[ran_log_list];  /* free flight path */
        if (++ran_log_list >= MAXLOGLIST) ran_log_list = 0;
    }
    if (*current_material_index == 0) {
        old_surf_bind = 0.0;
    }
    else {
        old_surf_bind = check_surf_bind (*current_material_index, projZ);
    }

    *stopping = 0.0;  *straggling = 0.0;
    s_temp    = 0.0;  *s = 0.0;

    if (max_no_layers == 0) {  /* half-infinite substrate */
        /* distance from the ion to the surface */
        T0 = dist_to_surf (z, sub_surf_z, vz);

        /* no interaction, go into the bulk */
        if (vz>1.e-35 || (z>sub_surf_z && fabs (vz)<1.e-35)) {
            if (flight_length_type == 1) {  /* 1 - constant */
                si = s_const;
            }
            else {  /* 0 or default - random */
                si = s_const / sect0;
            }

            if (T0 > 0.0) {
			    *s = T0 + si;
			    *current_material_index = 1;
			    current_surf_bind = check_surf_bind (*current_material_index, projZ);
			    *energy += current_surf_bind;
		    }
		    else {
			    *s = si;
		    }
            *ion_left_target = 1;
        }
        else {
            if (T0 > s_const) {
                if (flight_length_type == 1) {  /* 1 - constant */
                    si = s_const;
                }
                else {  /* 0 or default - random */
                    si = s_const / sect0;
                }

            	*ion_left_target = 1;
            }
            else {  /* escape */
            	si = T0;

            	*ion_left_target = 0;
            	*energy -= old_surf_bind;
            }
            *s = si;
        }
        e_energy_loss (projZ, projM, *energy, *current_material_index, si, stopping, straggling);
    }
    else {  /* multi-layers */
        /* determine current layer */
        j = 0;
        if (vz > 0.0) {
            for (i=0; i<=max_no_layers; i++) {
            	if (z < thick_of_layer[i]) {
            		j = i;
            		break;
            	}
            }
            if (i > max_no_layers) j = max_no_layers + 1;
        }
        else {
        	for (i=max_no_layers; i>=0; i--) {
        		if (z > thick_of_layer[i]) {
        			j = i;
        			break;
        		}
        	}
        }

        if (vz>1.e-35 || (z>sub_surf_z && fabs (vz)<1.e-35)) {
            for (i=j; i<=max_no_layers; i++) {
            	if (i == j) {
            		TT[i] = dist_to_surf (z, thick_of_layer[i], vz);
            	}
            	else {
            	    TT[i] = dist_to_surf (z, thick_of_layer[i], vz) - TT[i-1];
            	}

            	if (flight_length_type == 1) {  /* 1 - constant */
                    if (sect[i] > 0.0) s_temp += TT[i];
                }
                else {  /* 0 or default - random */
                    s_temp += TT[i] * sect[i];
                }

                if (s_temp >= s_const) {
                	if (flight_length_type == 1) {  /* 1 - constant */
                        si = s_const;
                    }
                    else {  /* 0 or default - random */
                        si = s_const / sect[i];
                    }

                    *s += si;
                    *current_material_index = material_of_layer[i];
                    *ion_left_target = 1;
                    e_energy_loss (projZ, projM, *energy, *current_material_index, si, stopping, straggling);
                    return;
                }
                else {
                	si = TT[i];
                	if (*current_material_index > 0)
                	    e_energy_loss (projZ, projM, *energy, *current_material_index, si, stopping, straggling);
                	*current_material_index = material_of_layer[i+1];
                	*ion_left_target = 1;
                	current_surf_bind = check_surf_bind (*current_material_index, projZ);
                    *energy += current_surf_bind - old_surf_bind;
                }
                *s += si;
            }

            if (flight_length_type == 1) {  /* 1 - constant */
            	si = (s_const - s_temp);
            }
            else {
            	si = (s_const - s_temp) / sect0;
            }
            *s += si;
            *current_material_index = 1;
            *ion_left_target = 1;
            e_energy_loss (projZ, projM, *energy, *current_material_index, si, stopping, straggling);
        }
        else {
            for (i=j; i>=0; i--) {
            	if (i == j) {
            		TT[i] = dist_to_surf (z, thick_of_layer[i], vz);
            	}
            	else {
            	    TT[i] = dist_to_surf (z, thick_of_layer[i], vz) - TT[i+1];
            	}
            	if (flight_length_type == 1) {  /* 1 - constant */
                    if (sect[i] > 0.0) s_temp += TT[i];
                }
                else {  /* 0 or default - random */
                    s_temp += TT[i] * sect[i+1];
                }

                if (s_temp >= s_const) {
                	if (flight_length_type == 1) {  /* 1 - constant */
                        si = s_const;
                    }
                    else {  /* 0 or default - random */
                        si = s_const / sect[i+1];
                    }

                    *s += si;
                    *current_material_index = material_of_layer[i+1];
                    *ion_left_target = 1;
                    e_energy_loss (projZ, projM, *energy, *current_material_index, si, stopping, straggling);
                    return;
                }
                else {
                	si = TT[i];
                	if (*current_material_index > 0)
                	    e_energy_loss (projZ, projM, *energy, *current_material_index, si, stopping, straggling);
                	*current_material_index = material_of_layer[i];
                	*ion_left_target = 1;
                	current_surf_bind = check_surf_bind (*current_material_index, projZ);
                    *energy = *energy - old_surf_bind + current_surf_bind;
                }
                *s += si;
            }
            *ion_left_target = 0;
            current_surf_bind = check_surf_bind (*current_material_index, projZ);
            *energy -= current_surf_bind;
        }
    }

    return;
}

/*=============================================================================
  Function Name : freepath_bulk
  Description   : Check if projectile is still in the bulk/layer.

  Inputs  : double z, double vz
            double flight_length
  Outputs : no.

  Notes :
      Function call:
        bulk - double dist_to_surf (double z, double z0, double vz);
=============================================================================*/
int check_target_type_bulk (double z, double vz, double flight_length) {
    int    i, j;
    double s = 0.0, T0;

    if (max_no_layers == 0) {  /* half-infinite substrate */
        if (vz > 0.0) {  /* no interaction, go into the bulk */
            return 1;
        }
        else {
            /* distance from the ion to the surface */
        	T0 = dist_to_surf (z, sub_surf_z, vz);
            if (T0 > flight_length) {
                return 1;
            }
            else {  /* escape */
                return 0;
            }
        }
    }
    else {  /* multi-layers */
        j = 0;
        if (vz > 0.0) {
            for (i=0; i<=max_no_layers; i++) {
            	if (z < thick_of_layer[i]) {
            		j = i;
            		break;
            	}
            }
            if (i > max_no_layers) j = max_no_layers + 1;
        }
        else {
        	for (i=max_no_layers; i>=0; i--) {
        		if (z > thick_of_layer[i]) {
        			j = i;
        			break;
        		}
        	}
        }

        s = dist_to_surf (z, thick_of_layer[j], vz);
        if (s > flight_length) {
            return 1;
        }
        else {
            return 0;
        }
    }

    return 0;
}

/*=============================================================================
  Function Name : dist_to_surf
  Description   : Calculate the distance from ion to the substrate surface.

  Inputs  :
            float x
            float vx
  Outputs : double T0

  Notes : no.
=============================================================================*/
double dist_to_surf (double z, double z0, double vz) {
    double T0;

    T0 = -1.e10;
    if (fabs (vz) > 1.0e-20) {  /* not parallel to the surface of substrate */
        T0 = - (z - z0) / vz;
    }
    else {  /* parallel to the surface of substrate */
        T0 = 1.e10;
    }

    return T0;
}
