#include <petsc.h>


PetscErrorCode formExact(DM da, Vec uexact){
    PetscInt i;
    PetscReal hx, x, *auexect;   
    DMDALocalInfo info;

    PetscCall(DMDAGetLocalInfo(da, &info));
    hx = 1.0/(info.mx-1);
    PetscCall(DMDAVecGetArray(da, uexact, &auexect));
    for( i = info.xs; i < info.xs+info.xm; i++){
        //x = hx*i;
        x=i;
        auexect[i] = x*x*(1-x*x);
    }
    PetscCall(DMDAVecRestoreArray(da, uexact, &auexect));
    return 0;
}

PetscErrorCode formRHS(DM da, Vec b){
    PetscInt i;
    PetscReal hx, x, *ab;
    DMDALocalInfo info;

    PetscCall(DMDAGetLocalInfo(da, &info));
    PetscCall(DMDAVecGetArray(da, b, &ab));
    
    hx = 1.0/(info.mx-1);

    for( i = info.xs; i < info.xs+info.xm; i++){
        if(i == 0 || i == info.mx-1){
            ab[i] = 0;
        }
        else{
            //x = i*hx;
            x=i;
            ab[i] = 1.0-6.0*(x*x);
        }
    } 
    PetscCall(DMDAVecRestoreArray(da, b, &ab));
    return 0;
}

PetscErrorCode formMat(DM da, Mat A){
    PetscInt i, ncols;
    PetscReal hx, x, v[3];
    MatStencil row, col[3];
    DMDALocalInfo info;

    PetscCall(DMDAGetLocalInfo(da, &info));
    hx = 1.0/(info.mx-1);

    for( i = info.xs; i < info.xs+info.xm; i++){
        row.i = i;
        col[0].i = i;
        ncols = 1;
        if (i == 0 || i == info.mx-1)       
        {
            v[0] = 1.0;
        }
        else
        {
            v[0] = -2/(i*i);
            if(i - 1 > 0){
                col[ncols].i = i-1;
                v[ncols] = 1/(i*i);
                ncols++; 
            }
            if(i+1 < info.mx-1){
                col[ncols].i = i+1;
                v[ncols] = 1/(i*i);
                ncols++;
            }
        }
        PetscCall(MatSetValuesStencil(A, 1, &row, ncols, col, v, INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
    return 0;
}

int main(int argc, char **argv) {
    DM da;
    Vec u, uexact, b;
    Mat A;
    KSP ksp;

    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
    
    //create and setup DM
    PetscCall(DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, 10, 1, 1, NULL, &da));
    PetscCall(DMSetFromOptions(da));

    //create and setup Mat
    PetscCall(DMCreateMatrix(da, &A));
    PetscCall(MatSetFromOptions(A));

    //create and setup Vecs
    PetscCall(DMCreateGlobalVector(da, &u));
    PetscCall(VecDuplicate(u, &b));
    PetscCall(VecDuplicate(u, &uexact));
    
    //Add value to Vecs and Mat
    formExact(da, uexact);
    formRHS(da, b);
    formMat(da, A);


    PetscCall(PetscFinalize());
    return 0;
}