/******************************************************************************
  Module Name : mpimod.h
  Module Date : 04/22/2014
  Module Auth : Yonggang Li

  Description : MPI parallel module.

  Others :
      Revision History:
      Date    Rel Ver.    Notes
******************************************************************************/
#include "mpimod.h"

/*=============================================================================
  Function Name : mpi_seed
  Description   : MPI random.

  Inputs  : no.
  Outputs : no.

  Notes : no.
=============================================================================*/
int mpi_seed (int seed) {
    return seed += my_node;
}

/*=============================================================================
  Function Name : mpi_init
  Description   : Arrays used for massage transfer.

  Inputs  : no.
  Outputs : no.

  Notes : no.
=============================================================================*/
int mpi_init_array (void) {

	  reduce_int_data    = (int*) calloc (cell_count, sizeof (int));
	  reduce_double_data = (double*) calloc (cell_count, sizeof (double));

	  if (reduce_int_data    == NULL) return -7001;  /* cannot allocate memory */
	  if (reduce_double_data == NULL) return -7002;  /* cannot allocate memory */

	  return 0;
}

/*=============================================================================
  Function Name : mpi_distribute
  Description   : Distribute jobs to nodes.

  Inputs  : no.
  Outputs : no.

  Notes : no.
=============================================================================*/
void mpi_distribute (void) {
    int num_ion_per_node, num_ion_last_node;

    num_ion_per_node = max_no_ions / num_nodes;
    num_ion_last_node = num_ion_per_node + (max_no_ions % num_nodes);
    if (my_node == num_nodes-1) {
        max_no_ions_node = num_ion_last_node;
    }
    else {
        max_no_ions_node = num_ion_per_node;
    }

	return;
}

/*=============================================================================
  Function Name : mpi_distribute
  Description   : Reduced final output data.

  Inputs  : no.
  Outputs : no.

  Notes : no.
=============================================================================*/
void mpi_reduce_data (void) {
    int i, j;

    MPI_Reduce (target_implanted_ions, reduce_int_data, cell_count, MPI_INT, MPI_SUM,
                ROOT, MPI_COMM_WORLD);
    if (my_node == ROOT) copy_int_array (reduce_int_data, target_implanted_ions, cell_count);
    MPI_Reduce (target_replacing_ions , reduce_int_data, cell_count, MPI_INT, MPI_SUM,
                ROOT, MPI_COMM_WORLD);
    if (my_node == ROOT) copy_int_array (reduce_int_data, target_replacing_ions, cell_count);
    MPI_Reduce (target_total_displacements, reduce_int_data, cell_count, MPI_INT, MPI_SUM,
                ROOT, MPI_COMM_WORLD);
    if (my_node == ROOT) copy_int_array (reduce_int_data, target_total_displacements, cell_count);
    MPI_Reduce (target_total_interstitials, reduce_int_data, cell_count, MPI_INT, MPI_SUM,
                ROOT, MPI_COMM_WORLD);
    if (my_node == ROOT) copy_int_array (reduce_int_data, target_total_interstitials, cell_count);
    MPI_Reduce (target_total_vacancies, reduce_int_data, cell_count, MPI_INT, MPI_SUM,
                ROOT, MPI_COMM_WORLD);
    if (my_node == ROOT) copy_int_array (reduce_int_data, target_total_vacancies, cell_count);
    MPI_Reduce (target_total_replacements, reduce_int_data, cell_count, MPI_INT, MPI_SUM,
                ROOT, MPI_COMM_WORLD);
    if (my_node == ROOT) copy_int_array (reduce_int_data, target_total_replacements, cell_count);

	  if (detailed_sputtering == 1) {  /* for detailed calucation of sputtering, we need these */
        MPI_Reduce (target_total_sputtered, reduce_int_data, cell_count, MPI_INT, MPI_SUM,
                    ROOT, MPI_COMM_WORLD);
        if (my_node == ROOT) copy_int_array (reduce_int_data, target_total_sputtered, cell_count);
    }

    if (store_energy_deposit == 1) {  /* for detailed storing of deposited energy, we need these */
        MPI_Reduce (target_energy_phonons, reduce_double_data, cell_count, MPI_DOUBLE,
                    MPI_SUM, ROOT, MPI_COMM_WORLD);
        if (my_node == ROOT) copy_double_array (reduce_double_data, target_energy_phonons, cell_count);
        MPI_Reduce (target_energy_electrons, reduce_double_data, cell_count, MPI_DOUBLE,
                    MPI_SUM, ROOT, MPI_COMM_WORLD);
        if (my_node == ROOT) copy_double_array (reduce_double_data, target_energy_electrons, cell_count);
    }

    if (single_ion_sputter_yields == 1) {
        MPI_Reduce (sputter_yield_histogram, reduce_int_data3, MAX_SPUTTERED+1, MPI_INT, MPI_SUM,
                    ROOT, MPI_COMM_WORLD);
        if (my_node == ROOT) copy_int_array (reduce_int_data3, sputter_yield_histogram, MAX_SPUTTERED+1);
    }

    if (detailed_sputtering == 1) {
        MPI_Reduce (total_sputter_counter, reduce_int_data2, 8, MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);
        if (my_node == ROOT) copy_int_array (reduce_int_data2, total_sputter_counter, 8);
        MPI_Reduce (leaving_ions, reduce_int_data2, 8, MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);
        if (my_node == ROOT) copy_int_array (reduce_int_data2, leaving_ions, 8);
    }

    /* go through materials, init arrays for interstitials and vacancies */
    for (i=0; i<number_of_materials; i++) {
        /* go through elements, make arrays for each one */
        for (j=0; j<list_of_materials[i].element_count; j++) {
            MPI_Reduce (list_of_materials[i].target_implanted_recoils_int[j], reduce_int_data, cell_count,
                        MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);
            if (my_node == ROOT)
                copy_int_array (reduce_int_data, list_of_materials[i].target_implanted_recoils_int[j], cell_count);
            MPI_Reduce (list_of_materials[i].target_implanted_recoils_repl[j], reduce_int_data, cell_count,
                        MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);
            if (my_node == ROOT)
                copy_int_array (reduce_int_data, list_of_materials[i].target_implanted_recoils_repl[j], cell_count);
            MPI_Reduce (list_of_materials[i].target_elemental_vacancies[j], reduce_int_data, cell_count,
                        MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);
            if (my_node == ROOT)
                copy_int_array (reduce_int_data, list_of_materials[i].target_elemental_vacancies[j], cell_count);
            MPI_Reduce (list_of_materials[i].target_elemental_disp[j], reduce_int_data, cell_count,
                        MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);
            if (my_node == ROOT)
                copy_int_array (reduce_int_data, list_of_materials[i].target_elemental_disp[j], cell_count);
            //MPI_Reduce (list_of_materials[i].elemental_leaving_recoils[j], reduce_int_data, cell_count,
            //            MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);
            //if (my_node == ROOT)
            //    copy_int_array (reduce_int_data, list_of_materials[i].elemental_leaving_recoils[j], cell_count);
        }
        MPI_Reduce (list_of_materials[i].sputter_counter, reduce_int_data4, MAX_EL_PER_MAT*8, MPI_INT, MPI_SUM,
                    ROOT, MPI_COMM_WORLD);
        if (my_node == ROOT)
            copy_int_array (reduce_int_data4, list_of_materials[i].sputter_counter, MAX_EL_PER_MAT*8);
    }

    free (reduce_int_data);
    free (reduce_double_data);

    return;
}
