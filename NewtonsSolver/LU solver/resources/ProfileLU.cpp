#include "../headers/ProfileLU.h"

inline void LU::ProfileSolver::_Direct(const ProfileMatrix& mat, std::vector<double>& x, const std::vector<double>& F) {
   for (size_t i = 0; i < mat.Size(); i++)
   {
      size_t j = i - (mat.ia[i + 1] - mat.ia[i]);
      double sum = 0;
      for (size_t k = mat.ia[i]; k < mat.ia[i + 1]; k++, j++)
      {
         sum += F[j] * mat.al[k];
      }
      x[i] = (F[i] - sum) / mat.diag[i];
   }
}

inline void LU::ProfileSolver::_Reverse(const ProfileMatrix& mat, std::vector<double>& x, const std::vector<double>& F) {
   for (size_t i = mat.Size(); i > 0; )
   {
      --i;
      size_t j = i - (mat.ia[i + 1] - mat.ia[i]);
      for (size_t k = mat.ia[i]; k < mat.ia[i + 1]; k++, j++)
      {
         x[j] -= x[i] * mat.au[k];
      }
   }
}
