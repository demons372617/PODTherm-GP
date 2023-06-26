#include "POD_ODE_Solver.hpp"

#include <petscmat.h>
#include <stdio.h>
#include <petscts.h>
#include <petscksp.h>
#include <unistd.h> // for sleep function

#include<sstream>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdraw.h>
#include <mpi.h>

namespace Thermal {
namespace POD {
        
/* import C, G, P_vector and temperature gradient to solve ODE.  And then save the coefficient of ppodmode*/

typedef struct {
    Mat A;                 /* RHS mat, used with IFunction interface */
    Vec vp;
} AppCtx;

// include user defined routine
extern PetscErrorCode RHSFunctionHeat(TS ts,PetscReal t,Vec X,Vec r,void *ctx);

bool ODE_Solver::init( int argc, char** argv ) {
    
    PetscInitialize(&argc,&argv,0,0);
    read_xml_config("../../POD_Para.xml");
    bool res = true;
    res &= C.read_file("../../C.csv",false);
    C.print_info("C");
    if(!res) { printf("After loading C!\n"); }
    res &= Gmatrix.read_file("../../G.csv",false);
    Gmatrix.print_info("G");
    if(!res) { printf("After loading G!\n"); }
    printf("FLP\n");
    res = res && load_file(floorplan_file ,MFloorplan,NU,4);
    if(!res) { printf("After loading FLP!\n"); }
    printf("power\n");
    res = res && load_file(ptrace_file,Pmatrix1,NT,NU); // power to be predicted
    if(!res) { printf("After loading power!\n"); }
    res &= P_lib.read_file("../../P_mat.csv",false);
    P_lib.print_info("P");
    if(!res) { printf("After loading P!\n"); }
   /* res &= temp_grad_lib.read_file("TG_mat.csv",false);
    temp_grad_lib.print_info("TG");
    if(!res) { printf("After loading TG!\n"); }
    res &= times_per_step.read_file("timings.csv",false);
    if(!res) { printf("After loading timings!\n"); }*/
    
    return res;
}

int ODE_Solver::run()
{
    // define the petsc parameter
    Mat   CP,GP,PP,PSP,RC,CG,CPSP; // context of eps  
    Vec              u,u0; // eigenvector
    MPI_Comm        comm;
    PetscInt        *dnnz,*onnz,i,j; // variable will be used in slepc
    /*****************************************************************************/
    IS             perm,iperm;
    MatFactorInfo  info;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    comm = MPI_COMM_WORLD;
    PetscErrorCode                 ierr;
    KSP                            ksp;
    PetscViewer                    viewer;
    
    /***********************************************************
     * compute the power density. Hereafter, Pmatrix is power density.
     * ***********************************************************/
    printf("PD\n");
    for(int i = 0; i < NT; i++) { //iterate all time steps
        for(int j =0; j < NU;  j++) { //iterate all functional units
            Pmatrix1[i][j] = Pmatrix1[i][j]/(MFloorplan[j][0]*MFloorplan[j][1]*thick_actl);             
            // printf("%lf\n",Pmatrix1[i][j]);
        }
    }
    /***********************************************************
     * compute the P vector.
     * ***********************************************************/
    double **PSmatrix = (double**)malloc(NT*sizeof(double*));
    if (!PSmatrix)
        return -1;
    for(int i = 0; i < NT; i++) {
        PSmatrix[i] = (double*)malloc(N_p*sizeof(double));
        if (!PSmatrix[i])
            return -1;
    }

    /**
    ::Matrix new_PS;
    new_PS.read_file("PS_4s.csv",false);
    **/
    double sum_p =0;
    for (int i = 0; i < NT; i++) {  //iterate time steps
        for (int j = 0; j < N_p; j++) { //iterate POD modes
            sum_p =0;
            for(int k =0; k < NU; k++) { //iterate function units
                sum_p =  sum_p +P_lib[j][k]*Pmatrix1[i][k];
            }
	    //std::cout << temp_grad_lib.getNRows() << "," << temp_grad_lib.getNCols() << ":" << i << "," << j << std::endl;
            //PSmatrix[i][j] = sum_p + temp_grad_lib[j][i]; // BAD
            PSmatrix[i][j] = sum_p; // GOOD
        }
    }


    printf("Starting Solver!\n");
    //sleep(60);
    for (int num_pod = 1; num_pod <= NPOD; num_pod++) {
        //std::cout << "STEP: " << num_pod << std::endl;
        PetscMalloc2(num_pod,&dnnz,num_pod,&onnz); // allocate 2 arrays of matrix
        MatCreateDense(comm,PETSC_DECIDE,PETSC_DECIDE,num_pod,num_pod,NULL,&CP);
        MatSetUp(CP);
        PetscFree2(dnnz,onnz);
        // printf("bp 0.1\n");
        MatSetOption(CP,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
        for (i=0; i<num_pod; i++) {
            for(j=0; j<num_pod;j++){
                MatSetValue(CP,i,j,C[i][j],INSERT_VALUES);
            }
        }
        MatAssemblyBegin(CP,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(CP,MAT_FINAL_ASSEMBLY);

        //MatGetOrdering(CP,MATORDERINGRCM,&perm,&iperm);
        //MatFactorInfoInitialize(&info);
        //info.fill          = 1.0;
        //MatLUFactor(CP,perm,iperm,&info);
        // printf("bp 0.2\n");

        // if(rank == 0){
        //     printf("CP view, rank %d\n", rank);
        //     MatView(CP, PETSC_VIEWER_STDOUT_(comm)); 
        // }
        
        
        PetscMalloc2(num_pod,&dnnz,num_pod,&onnz); // allocate 2 arrays of matrix
        MatCreateDense(comm,PETSC_DECIDE,PETSC_DECIDE,num_pod,num_pod,NULL,&RC); // RC is inverse matrix of C matrix
        MatSetUp(RC);
        PetscFree2(dnnz,onnz);
        MatSetOption(RC,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
        for (i=0; i<num_pod; i++) {
            for(j=0; j<num_pod;j++){
                if(i==j)
                    MatSetValue(RC,i,j,1,INSERT_VALUES);
                else
                    MatSetValue(RC,i,j,0,INSERT_VALUES);
            }
        }
        MatAssemblyBegin(RC,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(RC,MAT_FINAL_ASSEMBLY);
        

        ////////////////////////////////////////////////////////
        PetscMalloc2(num_pod,&dnnz,num_pod,&onnz); // allocate 2 arrays of matrix
        MatCreateDense(comm,PETSC_DECIDE,PETSC_DECIDE,num_pod,num_pod,NULL,&PP);
        MatSetUp(PP);
        PetscFree2(dnnz,onnz);
        MatSetOption(PP,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
        for (i=0; i<num_pod; i++) {
            for(j=0; j<num_pod;j++){
                if(i==j)
                    MatSetValue(PP,i,j,1,INSERT_VALUES);
                else
                    MatSetValue(PP,i,j,0,INSERT_VALUES);
            }
        }
        MatAssemblyBegin(PP,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(PP,MAT_FINAL_ASSEMBLY);
        // printf("bp 0.3\n");
        
        //printf("PP view, rank %d\n", rank);
        //MatView(PP, PETSC_VIEWER_STDOUT_(comm));

        
        
        KSPCreate(PETSC_COMM_WORLD,&ksp);
        KSPSetOperators(ksp,CP,CP);
        KSPSetFromOptions(ksp);
        KSPSetUp(ksp);
        // printf("bp 0\n");

        KSPMatSolve(ksp,PP,RC);
        // printf("bp 1\n");
        
        
        // MatView(RC, PETSC_VIEWER_STDOUT_(comm));

        //MatMatSolve(CP,PP,RC);

        /*-----------------------------------------------------------------------------------------------------
          save G matrix in petsc format GP
          -------------------------------------------------------------------------------------------------------*/
        PetscMalloc2(num_pod,&dnnz,num_pod,&onnz); // allocate 2 arrays of matrix
        MatCreateDense(comm,PETSC_DECIDE,PETSC_DECIDE,num_pod,num_pod,NULL,&GP);
        MatSetUp(GP);
        PetscFree2(dnnz,onnz);
        MatSetOption(GP,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
        for (i=0; i<num_pod; i++) {
            for(j=0; j<num_pod;j++){
                MatSetValue(GP,i,j,-Gmatrix[i][j],INSERT_VALUES);
            }
        }
        MatAssemblyBegin(GP,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(GP,MAT_FINAL_ASSEMBLY);
        //free Gmatrix

        /*------------------------------------------------------------------------------------------------------
          compute C inverse matrix times G matrix
          ---------------------------------------------------------------------------------------------------------*/
        // printf("bp 2\n");
        MatMatMult(RC,GP,MAT_INITIAL_MATRIX,1,&CG);
        // printf("bp 3\n");

        /*------------------------------------------------------------------------------------------------------
          save source term in petsc format matrix PSP
          ---------------------------------------------------------------------------------------------------------*/

        PetscMalloc2(num_pod,&dnnz,NT,&onnz); // allocate 2 arrays of matrix
        MatCreateDense(comm,PETSC_DECIDE,PETSC_DECIDE,num_pod,NT,NULL,&PSP);
        MatSetUp(PSP);
        PetscFree2(dnnz,onnz);
        MatSetOption(PSP,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);

        for (i=0; i<num_pod; i++) {
            for(j=0; j<NT;j++){
                MatSetValue(PSP,i,j,PSmatrix[j][i],INSERT_VALUES); // psp is transpose matrix of PSmatrix
            }
        }
        MatAssemblyBegin(PSP,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(PSP,MAT_FINAL_ASSEMBLY);

        // printf("bp 4\n");
        MatMatMult(RC,PSP,MAT_INITIAL_MATRIX,1,&CPSP);
        // printf("bp 5\n");

        /*-----------------------------------------------------------------------------------------------------
          save mat4 matrix in petsc format PHI
          -------------------------------------------------------------------------------------------------------*/
        /*-------------------------------------------------------------------------------------------
          SLOVE THE EQUATION
          -------------------------------------------------------------------------------------------------*/   
        AppCtx         appctx;                 /* user-defined application context */
        TS             ts;                     /* timestepping context */
        Vec            u_local,u0_local;                      /* approximate solution vector */
        PetscInt       time_steps_max = 100;   /* default max timesteps */
        PetscInt       u0_size, u_size_global;
        PetscReal      dt;
        PetscInt       global_size, local_size;


        VecCreate(PETSC_COMM_WORLD,&u);
        VecSetSizes(u,PETSC_DECIDE,num_pod);
        VecCreate(PETSC_COMM_WORLD,&u0);
        VecSetSizes(u0,PETSC_DECIDE,num_pod);
        VecCreate(PETSC_COMM_WORLD,&appctx.vp);
        VecSetSizes(appctx.vp,PETSC_DECIDE,num_pod);
        VecSetType(appctx.vp, VECMPI);
        VecSetType(u, VECMPI);
        VecSetType(u0, VECMPI);

        VecZeroEntries(u);
        VecZeroEntries(u0);
        

        VecZeroEntries(appctx.vp);

        MatCreate(PETSC_COMM_SELF,&appctx.A);
        MatSetSizes(appctx.A,PETSC_DECIDE,PETSC_DECIDE,num_pod,num_pod);
        MatSetFromOptions(appctx.A);
        MatSetUp(appctx.A);
        MatDuplicate(CG, MAT_COPY_VALUES,&appctx.A);

        for (int i = 0; i < NT; i ++){
	        //sampling_interval = times_per_step[i][0];

            TSCreate(PETSC_COMM_WORLD,&ts);
            TSSetProblemType(ts,TS_LINEAR);
            TSSetType(ts,TSRK);
            TSRKSetType( ts,TSRK4);
            TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP );
            TSSetRHSJacobian(ts,CG,CG,TSComputeRHSJacobianConstant,&appctx);

            TSSetMaxSteps(ts,time_steps_max);
            TSSetMaxTime(ts,sampling_interval);

            // set solution vector and initial timestepping
            dt = sampling_interval/20;
            TSSetTimeStep(ts,dt);

            MatGetColumnVector(CPSP,appctx.vp,i);

            //VecCopy(vp,u);
            TSSetRHSFunction(ts,NULL,RHSFunctionHeat,&appctx);
            TSSetSolution(ts,u0);
            double timestamp_before = MPI_Wtime();
            printf("Timestamp: %f, rank %d, front\n", timestamp_before, rank);
            TSSolve(ts,u);
            // VecView(u, PETSC_VIEWER_STDOUT_WORLD);
            double timestamp_after = MPI_Wtime();
            printf("Timestamp: %f, rank %d, back\n", timestamp_after, rank);
            VecZeroEntries(u0);
            VecCopy(u,u0);
            
            // printf("bk check, rank %d\n", rank);
            VecGetSize(u, &u_size_global);
            // printf("u global size is %d, rank %d\n", u_size_global, rank);
            VecCreate(PETSC_COMM_SELF,&u_local);
            VecSetSizes(u_local,u_size_global,u_size_global);
            VecSetType(u_local, VECSEQ);
            VecZeroEntries(u_local);
            // printf("bk check2, rank %d\n", rank);
            IS is_global;
            ISCreateStride(PETSC_COMM_WORLD, u_size_global, 0, 1, &is_global);
            // printf("bk check3, rank %d\n", rank);
            VecScatter ctx;
            VecScatterCreate(u, is_global, u_local, NULL, &ctx);
            // printf("bk check4, rank %d\n", rank);
            ierr = VecScatterBegin(ctx, u, u_local, INSERT_VALUES, SCATTER_FORWARD);
            // printf("bk check5, rank %d\n", rank);
            CHKERRQ(ierr);
            ierr = VecScatterEnd(ctx, u, u_local, INSERT_VALUES, SCATTER_FORWARD);
            // printf("bk check6, rank %d\n", rank);
            CHKERRQ(ierr);

            VecScatterDestroy(&ctx);
            ISDestroy(&is_global);

            MPI_Barrier(MPI_COMM_WORLD);
            // VecView(u_local, PETSC_VIEWER_STDOUT_SELF);
            if(rank == 0){
                for( int k = 0; k < u_size_global; k++) {
                    VecGetValues(u_local,1,&k,&CU[k][i]);
                }
            }
            
            MPI_Barrier(MPI_COMM_WORLD);
            
            
            TSDestroy(&ts);
            
        }
        
        VecDestroy(&u_local);
        if(rank == 0){
            std::cout << "WRITING CU " << num_pod << std::endl;
            char namecu[200];
            sprintf(namecu, "../../pod_result/CU%d.txt", num_pod);
            FILE *fileMAT4 = fopen(namecu,"w");

            for(int j = 0; j < NT; j++){
                for (int i = 0; i < num_pod; i++){
                    fprintf(fileMAT4,"%.16lg\t",CU[i][j]);
                }
                fprintf(fileMAT4,"\n");
            }
            fclose(fileMAT4); 
        }
        
        
        MatDestroy(&CP);
        MatDestroy(&RC);
        MatDestroy(&PP);
        MatDestroy(&GP);
        MatDestroy(&CG);
        MatDestroy(&PSP);
        MatDestroy(&CPSP);
        MatDestroy(&appctx.A);
        VecDestroy(&u);
        VecDestroy(&u0);
        VecDestroy(&appctx.vp);
    }

    // free memory
    // for (int i= 0 ; i < NC; i++){
    //     free(C[i]);
    // }
    // free(C);

    // for (int i= 0 ; i < NC; i++){
    //     free(Gmatrix[i]);
    // }
    // free(Gmatrix);

    for (int i= 0 ; i < NU; i++){
        free(MFloorplan[i]);
    }
    free(MFloorplan);

    for (int i= 0 ; i < NT; i++){
        free(Pmatrix1[i]);
    }
    free(Pmatrix1);

    // for (int i= 0 ; i < N_p; i++){
    //     free(P_lib[i]);
    // }
    // free(P_lib);

    // for (int i= 0 ; i < NT; i++){
    //     free(temp_grad_lib[i]);
    // }
    // free(temp_grad_lib);

    for (int i= 0 ; i < NT; i++){
        free(PSmatrix[i]);
    }
    free(PSmatrix);
    PetscFinalize();
    return 0;
}

PetscErrorCode  RHSFunctionHeat(TS ts,PetscReal t,Vec X,Vec r,void *ctx)
{

    // Mat A;
    // PetscFunctionBeginUser;
    // PetscCall(TSGetRHSJacobian(ts, &A, NULL, NULL, &ctx));
    // PetscCall(RHSMatrixHeat(ts, t, X, A, NULL, ctx));
    // PetscCall(MatView(A,PETSC_VIEWER_STDOUT_WORLD));
    // PetscCall(MatMult(A, X, r));

    AppCtx *appctx = (AppCtx*)ctx;     /* user-defined application context */
    // printf("bp omg\n");
    // MatView(appctx -> A, PETSC_VIEWER_STDOUT_(MPI_COMM_WORLD));
    // PetscInt appctx_rows, appctx_cols;
    // PetscInt vec_size;
    // MatGetLocalSize(appctx -> A, &appctx_rows, &appctx_cols);
    // printf("rows: %d, cols %d\n", appctx_rows, appctx_cols);
    // VecGetLocalSize(X, &vec_size);
    // printf("Vec x size: %d\n", vec_size);
    MatMult(appctx -> A,X,r);
    VecAXPY(r,1.0,appctx -> vp);

    return 0;
}

bool ODE_Solver::load_file(
    const std::string& ifname,
    double**& mat,
    const unsigned& nrows,
    const unsigned& ncols
) {
    mat = (double**)malloc(nrows*sizeof(double*));
    if(!mat) { return false; }
    for(int i = 0; i < nrows; i++){
        mat[i] = (double*)malloc(ncols*sizeof(double));
        if (!mat[i]) { return false; }
    }

    FILE* ifh = fopen(ifname.c_str(),"r");
    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            mat[i][j]=0;
            if(!fscanf(ifh , "%lf", &mat[i][j])) { break; }
            //printf("%lf\n",mat[i][j]);
        }
    }
    fclose(ifh);
    return true;
}


} // ns POD
} // ns Thermal
