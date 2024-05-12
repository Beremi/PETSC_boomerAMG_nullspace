#include <petscsys.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscviewerhdf5.h>
#include <petscksp.h>

int main(int argc, char **args)
{
    PetscInitialize(&argc, &args, NULL, NULL);

    Mat A, U;                                            // Matrix
    Vec b, indices_vec, indptr_vec, data_vec, shape_vec; // Vectors
    PetscViewer viewer;                                  // HDF5 viewer
    PetscInt nrows, ncols;
    PetscScalar *data, *indices_f, *indptr_f;
    PetscInt *indices, *indptr, num_indices, num_indptr, i, j;
    PetscErrorCode ierr;

    // Open HDF5 file
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, "matrix_and_vector2.h5", FILE_MODE_READ, &viewer));

    // Load the CSR components
    // Data
    PetscCall(VecCreate(PETSC_COMM_WORLD, &data_vec));
    PetscCall(PetscObjectSetName((PetscObject)data_vec, "MyMatrix/data"));
    PetscCall(VecLoad(data_vec, viewer));
    PetscCall(VecGetArray(data_vec, &data));

    // Indices
    PetscCall(VecCreate(PETSC_COMM_WORLD, &indices_vec));
    PetscCall(PetscObjectSetName((PetscObject)indices_vec, "MyMatrix/indices"));
    PetscCall(VecLoad(indices_vec, viewer));
    /* do some ugly calculations */
    PetscCall(VecGetSize(indices_vec, &num_indices));
    PetscCall(VecGetArray(indices_vec, &indices_f));
    /* and now we have to convert the double representation to integers to pass over, argh */
    PetscCall(PetscMalloc1(num_indices, &indices));
    for (j = 0; j < num_indices; j++)
        indices[j] = (PetscInt)indices_f[j];

    // Indptr
    PetscCall(VecCreate(PETSC_COMM_WORLD, &indptr_vec));
    PetscCall(PetscObjectSetName((PetscObject)indptr_vec, "MyMatrix/indptr"));
    PetscCall(VecLoad(indptr_vec, viewer));
    PetscCall(VecGetSize(indptr_vec, &num_indptr));
    PetscCall(VecGetArray(indptr_vec, &indptr_f));
    PetscCall(PetscMalloc1(num_indptr, &indptr));
    for (j = 0; j < num_indptr; j++)
        indptr[j] = (PetscInt)indptr_f[j];

    // Shape
    PetscCall(VecCreate(PETSC_COMM_WORLD, &shape_vec));
    PetscCall(PetscObjectSetName((PetscObject)shape_vec, "MyMatrix/shape"));
    PetscCall(VecLoad(shape_vec, viewer));
    PetscScalar *shape;
    PetscCall(VecGetArray(shape_vec, &shape));
    nrows = (PetscInt)shape[0];
    ncols = (PetscInt)shape[1];

    // Print the first and last 10 entries of indptr, indices, and data
    // print_vectors(data, indices, indptr, nrows, ncols);

    // Create and assemble the matrix A
    PetscCall(MatCreateSeqAIJ(PETSC_COMM_WORLD, nrows, ncols, PETSC_DEFAULT, NULL, &A));
    PetscCall(MatSeqAIJSetPreallocationCSR(A, indptr, indices, data));
    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

    // Load the RHS vector b
    PetscCall(VecCreate(PETSC_COMM_WORLD, &b));
    PetscCall(PetscObjectSetName((PetscObject)b, "MyMatrix/rhs"));
    PetscCall(VecLoad(b, viewer));

    // Load the null space matrix U

    // Open HDF5 file
    PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_HDF5_MAT));

    // Load the dense matrix U
    PetscCall(MatCreate(PETSC_COMM_WORLD, &U));
    PetscCall(MatSetType(U, MATDENSE));
    PetscCall(PetscObjectSetName((PetscObject)U, "MyMatrix/null_space"));
    PetscCall(MatLoad(U, viewer));

    // Print the shape of U
    PetscInt U_rows, U_cols;
    MatGetSize(U, &U_cols, &U_rows);
    PetscPrintf(PETSC_COMM_WORLD, "Shape of U: %D x %D\n", U_rows, U_cols);

    // Set matrix A near null space using dense matrix U
    PetscScalar *U_data;
    PetscCall(MatDenseGetArray(U, &U_data));
    MatNullSpace matnullspace;
    Vec *nullvecs;
    PetscReal dot_product;
    PetscCall(PetscMalloc1(U_cols, &nullvecs));
    for (j = 0; j < U_rows; j++)
    {
        PetscCall(VecCreate(PETSC_COMM_WORLD, &nullvecs[j]));
        PetscCall(VecSetSizes(nullvecs[j], PETSC_DECIDE, U_cols));
        PetscCall(VecSetFromOptions(nullvecs[j]));
        PetscCall(VecPlaceArray(nullvecs[j], U_data + j * U_cols));
    }
    for (i = 0; i < U_rows; i++)
    {
        for (j = 0; j < U_rows; j++)
        {
            PetscCall(VecDot(nullvecs[j], nullvecs[i], &dot_product));
            PetscPrintf(PETSC_COMM_WORLD, "%g, ", PetscRealPart(dot_product));
        }
        PetscPrintf(PETSC_COMM_WORLD, "\n");
    }

    PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, U_rows, nullvecs, &matnullspace));
    PetscCall(MatSetNearNullSpace(A, matnullspace));

    // Use A, b, U for computation or output results
    // Create the KSP solver
    KSP ksp;
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, A, A));
    PetscCall(KSPSetType(ksp, KSPCG));                                                   // Choose the KSP solver type, e.g., KSPCG for conjugate gradient method
    PetscCall(KSPSetTolerances(ksp, 1e-3, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT)); // Set the convergence criteria
    PetscCall(KSPSetFromOptions(ksp));                                                   // Set the options from the command line
    PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));                               // Set the initial guess to be nonzero
    PetscViewerAndFormat *vf;
    PetscCall(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf));
    PetscCall(KSPMonitorSet(ksp, (PetscErrorCode(*)(KSP, PetscInt, PetscReal, void *))KSPMonitorSingularValue,
                            vf, (PetscErrorCode(*)(void **))PetscViewerAndFormatDestroy));

    // Set the preconditioner
    PC pc;
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCHYPRE));
    PetscCall(PCHYPRESetType(pc, "boomeramg"));

    PetscCall(PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_vec_interp_variant", "3"));
    PetscCall(PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_nodal_coarsen", "6"));

    PetscCall(PCSetFromOptions(pc));

    // Solve the linear system
    Vec x;
    PetscCall(VecDuplicate(b, &x));
    PetscCall(KSPSolve(ksp, b, x));

    // Print the solution
    // PetscScalar *solution;
    // PetscInt solution_length, i;
    // PetscCall(VecGetArray(x, &solution));
    // PetscCall(VecGetLocalSize(x, &solution_length));
    // PetscPrintf(PETSC_COMM_WORLD, "Solution:\n");
    // for (i = 0; i < solution_length; i++)
    // {
    //     PetscPrintf(PETSC_COMM_WORLD, "%g\n", PetscRealPart(solution[i]));
    // }

    // Cleanup

    PetscCall(KSPDestroy(&ksp));
    PetscCall(VecDestroy(&x));
    PetscCall(MatDestroy(&U));
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(PetscFinalize());
    return 0;
}
