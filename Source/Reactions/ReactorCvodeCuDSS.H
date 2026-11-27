#ifndef REACTORCVODECUDSS_H
#define REACTORCVODECUDSS_H

// Batched sparse direct solves with NVIDIA cuDSS, used in place of the
// deprecated cuSolverSp batched QR when PelePhysics is built with
// PELE_USE_CUDSS. All cells share the sparsity pattern of the Newton (or
// preconditioner) matrix, so they are solved as one cuDSS uniform batch.

#ifdef PELE_USE_CUDSS

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <cuda_runtime.h>
#include <cudss.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#if !defined(CUDSS_VERSION) || CUDSS_VERSION < 10000
#error "PELE_USE_CUDSS requires cuDSS 1.0.0 or newer"
#endif

namespace pele::physics::reactions::cvode {

// LU factorization and solve of a uniform batch of CSR matrices: one
// row-pointer/column-index pattern (device memory, length n+1 and nnz) and
// one device array of nnz * nbatch values, system k at offset k * nnz.
// Right-hand sides and solutions are device arrays of n * nbatch values.
class CuDSSBatchSolver
{
public:
  CuDSSBatchSolver(
    int n,
    int nnz,
    int nbatch,
    int* rowptr_d,
    int* colind_d,
    amrex::Real* values_d,
    cudssIndexBase_t base,
    cudaStream_t stream);
  ~CuDSSBatchSolver();

  CuDSSBatchSolver(const CuDSSBatchSolver&) = delete;
  CuDSSBatchSolver& operator=(const CuDSSBatchSolver&) = delete;
  CuDSSBatchSolver(CuDSSBatchSolver&&) = delete;
  CuDSSBatchSolver& operator=(CuDSSBatchSolver&&) = delete;

  // Factorize the current values (reordering and symbolic analysis run on
  // the first call only, since the pattern is fixed). Asynchronous on the
  // stream given at construction.
  cudssStatus_t factorize(amrex::Real* values_d);

  // Solve with the last factorization: x = A^{-1} b. x and b must not
  // alias. Asynchronous on the stream given at construction.
  cudssStatus_t solve(amrex::Real* x_d, amrex::Real* b_d);

private:
  cudssHandle_t m_handle = nullptr;
  cudssConfig_t m_config = nullptr;
  cudssData_t m_data = nullptr;
  cudssMatrix_t m_A = nullptr;
  cudssMatrix_t m_x = nullptr;
  cudssMatrix_t m_b = nullptr;
  bool m_analyzed = false;
};

// SUNLinearSolver wrapping CuDSSBatchSolver for a SUNMatrix_cuSparse block-CSR
// matrix and CUDA N_Vectors: the matrix is factorized in Setup and the factors
// are reused for every Solve until CVODE calls Setup again.
struct SUNLinearSolverContent_cuDSS_struct
{
  CuDSSBatchSolver* solver;
  int last_flag;
};

using SUNLinearSolverContent_cuDSS = SUNLinearSolverContent_cuDSS_struct*;

SUNLinearSolver SUNLinSol_cuDSS_batch(
  N_Vector y, SUNMatrix A, cudaStream_t stream, SUNContext sunctx);

SUNLinearSolver_Type SUNLinSolGetType_cuDSS_batch(SUNLinearSolver S);

SUNLinearSolver_ID SUNLinSolGetID_cuDSS_batch(SUNLinearSolver S);

int SUNLinSolSetup_cuDSS_batch(SUNLinearSolver S, SUNMatrix A);

int SUNLinSolSolve_cuDSS_batch(
  SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, amrex::Real tol);

sunindextype SUNLinSolLastFlag_cuDSS_batch(SUNLinearSolver S);

int SUNLinSolFree_cuDSS_batch(SUNLinearSolver S);

} // namespace pele::physics::reactions::cvode

#endif
#endif
