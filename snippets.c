#include <petscsys.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscviewerhdf5.h>
#include <petscksp.h>


void print_vectors(PetscScalar *data, PetscInt *indices, PetscInt *indptr, PetscInt nrows, PetscInt ncols)
{
    /////
    // Assuming indptr, indices, and data are already loaded and available
    PetscInt i;
    PetscInt num_entries = 10; // Number of entries to print from the start and end

    PetscPrintf(PETSC_COMM_WORLD, "rows %D\n", nrows);

    PetscPrintf(PETSC_COMM_WORLD, "cols %D\n", ncols);
    // Print the first 10 entries of indptr, indices, and data
    PetscPrintf(PETSC_COMM_WORLD, "First 10 entries of indptr:\n");
    for (i = 0; i < num_entries; i++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "%D\n", indptr[i]);
    }

    PetscPrintf(PETSC_COMM_WORLD, "First 10 entries of indices:\n");
    for (i = 0; i < num_entries; i++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "%D\n", indices[i]);
    }

    PetscPrintf(PETSC_COMM_WORLD, "First 10 entries of data:\n");
    for (i = 0; i < num_entries; i++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "%g\n", PetscRealPart(data[i]));
    }

    // Print the last 10 entries of indptr, indices, and data
    PetscInt indptr_length, indices_length, data_length;
    // These values should be known or computed based on your application logic.
    // Example assignments for illustration:
    indptr_length = 77518;    // You must set this value appropriately
    indices_length = 3081123; // You must set this value appropriately
    data_length = 3081123;    // You must set this value appropriately

    PetscPrintf(PETSC_COMM_WORLD, "Last 10 entries of indptr:\n");
    for (i = indptr_length - num_entries; i < indptr_length; i++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "%D\n", indptr[i]);
    }

    PetscPrintf(PETSC_COMM_WORLD, "Last 10 entries of indices:\n");
    for (i = indices_length - num_entries; i < indices_length; i++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "%D\n", indices[i]);
    }

    PetscPrintf(PETSC_COMM_WORLD, "Last 10 entries of data:\n");
    for (i = data_length - num_entries; i < data_length; i++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "%g\n", PetscRealPart(data[i]));
    }

    ////
}