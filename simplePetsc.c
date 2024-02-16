#include <petsc.h>

int main(int argc, char const *argv[]) {
  PetscMPIInt rank;
  PetscInt i;
  PetscReal localval, globalsum;

  PetscCall(PetscInitialize(&argc, &argv, NULL, "Compute e in parallel with PETSc.\n\n"));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

  localval = 1.0;
  for (i = 2; i < rank + 1; i++)
    localval /= i;

  PetscCallMPI(MPI_Allreduce(&localval, &globalsum, 1, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "e is about %17.15f\n", globalsum));
  PetscCall(PetscPrintf(PETSC_COMM_SELF, "rank %d did %d flops\n", rank, (rank > 0) ? rank - 1 : 0));
  PetscCall(PetscFinalize());
  return 0;
}
