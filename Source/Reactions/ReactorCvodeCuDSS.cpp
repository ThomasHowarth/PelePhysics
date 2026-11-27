#include "ReactorCvodeCuDSS.H"

#ifdef PELE_USE_CUDSS

#include <type_traits>
#include <string>
#include <AMReX.H>
#include <AMReX_BLassert.H>
#include <AMReX_Print.H>
#include <AMReX_BLProfiler.H>
#include <nvector/nvector_cuda.h>
#include <sunmatrix/sunmatrix_cusparse.h>

namespace pele::physics::reactions::cvode {

namespace {

constexpr cudaDataType_t
cudssRealType()
{
  return std::is_same_v<amrex::Real, float> ? CUDA_R_32F : CUDA_R_64F;
}

void
cudssCheck(cudssStatus_t status, const char* what)
{
  if (status != CUDSS_STATUS_SUCCESS) {
    amrex::Abort(
      std::string("cuDSS: ") + what + " failed with status " +
      std::to_string(static_cast<int>(status)));
  }
}

} // namespace

CuDSSBatchSolver::CuDSSBatchSolver(
  const int n,
  const int nnz,
  const int nbatch,
  int* rowptr_d,
  int* colind_d,
  amrex::Real* values_d,
  const cudssIndexBase_t base,
  cudaStream_t stream)
{
  cudssCheck(cudssCreate(&m_handle), "cudssCreate");
  // The data object binds to this stream at analysis; later phases use it too
  cudssCheck(cudssSetStream(m_handle, stream), "cudssSetStream");
  cudssCheck(cudssConfigCreate(&m_config), "cudssConfigCreate");
  cudssCheck(cudssDataCreate(m_handle, &m_data), "cudssDataCreate");

  // One sparsity pattern shared by all cells: a uniform batch. This must be
  // set before the analysis phase.
  int ubatch_size = nbatch;
  cudssCheck(
    cudssConfigSet(
      m_config, CUDSS_CONFIG_UBATCH_SIZE, &ubatch_size, sizeof(ubatch_size)),
    "cudssConfigSet(UBATCH_SIZE)");

  cudssCheck(
    cudssMatrixCreateCsr(
      &m_A, n, n, nnz, rowptr_d, nullptr, colind_d, values_d, CUDA_R_32I,
      CUDA_R_32I, cudssRealType(), CUDSS_MTYPE_GENERAL, CUDSS_MVIEW_FULL, base),
    "cudssMatrixCreateCsr");
  // Right-hand side and solution pointers are set at each solve
  cudssCheck(
    cudssMatrixCreateDn(
      &m_x, n, 1, n, nullptr, cudssRealType(), CUDSS_LAYOUT_COL_MAJOR),
    "cudssMatrixCreateDn(x)");
  cudssCheck(
    cudssMatrixCreateDn(
      &m_b, n, 1, n, nullptr, cudssRealType(), CUDSS_LAYOUT_COL_MAJOR),
    "cudssMatrixCreateDn(b)");
}

CuDSSBatchSolver::~CuDSSBatchSolver()
{
  cudssMatrixDestroy(m_b);
  cudssMatrixDestroy(m_x);
  cudssMatrixDestroy(m_A);
  cudssDataDestroy(m_handle, m_data);
  cudssConfigDestroy(m_config);
  cudssDestroy(m_handle);
}

cudssStatus_t
CuDSSBatchSolver::factorize(amrex::Real* values_d)
{
  BL_PROFILE("Pele::CuDSSBatchSolver::factorize()");
  cudssStatus_t status = cudssMatrixSetValues(m_A, values_d);
  if (status != CUDSS_STATUS_SUCCESS) {
    return status;
  }
  if (!m_analyzed) {
    // Reordering (on the host) and symbolic factorization; synchronous
    status = cudssExecute(
      m_handle, CUDSS_PHASE_ANALYSIS, m_config, m_data, m_A, m_x, m_b);
    if (status != CUDSS_STATUS_SUCCESS) {
      return status;
    }
    m_analyzed = true;
  }
  status = cudssExecute(
    m_handle, CUDSS_PHASE_FACTORIZATION, m_config, m_data, m_A, m_x, m_b);
#ifdef AMREX_DEBUG
  if (status == CUDSS_STATUS_SUCCESS) {
    // Device-side errors surface through CUDSS_DATA_INFO; reading it
    // synchronizes, so only do so in debug builds
    int info = 0;
    size_t written = 0;
    status = cudssDataGet(
      m_handle, m_data, CUDSS_DATA_INFO, &info, sizeof(info), &written);
    if (status == CUDSS_STATUS_SUCCESS && info != 0) {
      amrex::Print() << "cuDSS factorization reported info = " << info << "\n";
    }
  }
#endif
  return status;
}

cudssStatus_t
CuDSSBatchSolver::solve(amrex::Real* x_d, amrex::Real* b_d)
{
  BL_PROFILE("Pele::CuDSSBatchSolver::solve()");
  AMREX_ASSERT(m_analyzed);
  AMREX_ASSERT(x_d != b_d);
  cudssStatus_t status = cudssMatrixSetValues(m_x, x_d);
  if (status != CUDSS_STATUS_SUCCESS) {
    return status;
  }
  status = cudssMatrixSetValues(m_b, b_d);
  if (status != CUDSS_STATUS_SUCCESS) {
    return status;
  }
  return cudssExecute(
    m_handle, CUDSS_PHASE_SOLVE, m_config, m_data, m_A, m_x, m_b);
}

#define SUN_CUDSS_CONTENT(S) ((SUNLinearSolverContent_cuDSS)((S)->content))
#define SUN_CUDSS_SOLVER(S) (SUN_CUDSS_CONTENT(S)->solver)
#define SUN_CUDSS_LASTFLAG(S) (SUN_CUDSS_CONTENT(S)->last_flag)

SUNLinearSolver
SUNLinSol_cuDSS_batch(
  N_Vector y, SUNMatrix A, cudaStream_t stream, SUNContext sunctx)
{
  if (y == nullptr || A == nullptr) {
    return nullptr;
  }
  if (N_VGetVectorID(y) != SUNDIALS_NVEC_CUDA) {
    return nullptr;
  }
  if (SUNMatGetID(A) != SUNMATRIX_CUSPARSE) {
    return nullptr;
  }
  if (SUNMatrix_cuSparse_SparseType(A) != SUNMAT_CUSPARSE_BCSR) {
    return nullptr;
  }
  if (N_VGetLength(y) != SUNMatrix_cuSparse_Columns(A)) {
    return nullptr;
  }

  SUNLinearSolver S = SUNLinSolNewEmpty(sunctx);
  if (S == nullptr) {
    return nullptr;
  }
  S->ops->gettype = SUNLinSolGetType_cuDSS_batch;
  S->ops->getid = SUNLinSolGetID_cuDSS_batch;
  S->ops->setup = SUNLinSolSetup_cuDSS_batch;
  S->ops->solve = SUNLinSolSolve_cuDSS_batch;
  S->ops->lastflag = SUNLinSolLastFlag_cuDSS_batch;
  S->ops->free = SUNLinSolFree_cuDSS_batch;

  auto* content = new SUNLinearSolverContent_cuDSS_struct;
  content->last_flag = SUN_SUCCESS;
  // SUNMatrix_cuSparse block-CSR storage is zero-based, with one pattern for
  // all blocks and the block values stored contiguously
  content->solver = new CuDSSBatchSolver(
    SUNMatrix_cuSparse_BlockRows(A), SUNMatrix_cuSparse_BlockNNZ(A),
    SUNMatrix_cuSparse_NumBlocks(A), SUNMatrix_cuSparse_IndexPointers(A),
    SUNMatrix_cuSparse_IndexValues(A), SUNMatrix_cuSparse_Data(A),
    CUDSS_BASE_ZERO, stream);
  S->content = content;

  return S;
}

SUNLinearSolver_Type
SUNLinSolGetType_cuDSS_batch(SUNLinearSolver /*S*/)
{
  return SUNLINEARSOLVER_DIRECT;
}

SUNLinearSolver_ID
SUNLinSolGetID_cuDSS_batch(SUNLinearSolver /*S*/)
{
  return SUNLINEARSOLVER_CUSTOM;
}

int
SUNLinSolSetup_cuDSS_batch(SUNLinearSolver S, SUNMatrix A)
{
  const cudssStatus_t status =
    SUN_CUDSS_SOLVER(S)->factorize(SUNMatrix_cuSparse_Data(A));
  SUN_CUDSS_LASTFLAG(S) =
    (status == CUDSS_STATUS_SUCCESS) ? SUN_SUCCESS : SUN_ERR_EXT_FAIL;
  return SUN_CUDSS_LASTFLAG(S);
}

int
SUNLinSolSolve_cuDSS_batch(
  SUNLinearSolver S,
  SUNMatrix /*A*/,
  N_Vector x,
  N_Vector b,
  amrex::Real /*tol*/)
{
  if (SUN_CUDSS_LASTFLAG(S) != SUN_SUCCESS) {
    return SUN_CUDSS_LASTFLAG(S);
  }
  const cudssStatus_t status = SUN_CUDSS_SOLVER(S)->solve(
    N_VGetDeviceArrayPointer_Cuda(x), N_VGetDeviceArrayPointer_Cuda(b));
  SUN_CUDSS_LASTFLAG(S) =
    (status == CUDSS_STATUS_SUCCESS) ? SUN_SUCCESS : SUN_ERR_EXT_FAIL;
  return SUN_CUDSS_LASTFLAG(S);
}

sunindextype
SUNLinSolLastFlag_cuDSS_batch(SUNLinearSolver S)
{
  return SUN_CUDSS_LASTFLAG(S);
}

int
SUNLinSolFree_cuDSS_batch(SUNLinearSolver S)
{
  if (S == nullptr) {
    return SUN_SUCCESS;
  }
  if (S->content != nullptr) {
    delete SUN_CUDSS_SOLVER(S);
    delete SUN_CUDSS_CONTENT(S);
    S->content = nullptr;
  }
  SUNLinSolFreeEmpty(S);
  return SUN_SUCCESS;
}

} // namespace pele::physics::reactions::cvode

#endif
