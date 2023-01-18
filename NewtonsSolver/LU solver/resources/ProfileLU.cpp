#include "../headers/ProfileLU.h"

void LU::ProfileSolver::_Direct(const ProfileMatrix& mat, std::vector<double>& x) {
   for (size_t i = 0; i < mat.Size(); i++)
   {
      size_t j = i - (mat.ia[i + 1] - mat.ia[i]);
      double sum = 0;
      for (size_t k = mat.ia[i]; k < mat.ia[i + 1]; k++, j++)
      {
         sum += x[j] * mat.al[k];
      }
      x[i] = (x[i] - sum) / mat.diag[i];
   }
}

void LU::ProfileSolver::_Reverse(const ProfileMatrix& mat, std::vector<double>& x) {
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
