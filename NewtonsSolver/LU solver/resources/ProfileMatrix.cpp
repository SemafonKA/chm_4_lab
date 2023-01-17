#include "../headers/ProfileMatrix.h"

inline void ProfileMatrix::MakeFromMatrix(const Matrix& mat) {
   if (mat.Cols() != mat.Rows())
      throw std::runtime_error("Bad matrix sizes (matrix should be squared)");

   diag.resize(mat.Rows());
   ia.resize(mat.Rows() + 1);
   std::size_t asize = 0;

   for (std::size_t i = 1; i < mat.Rows(); i++)
   {
      bool flag = false;
      for (std::size_t j = 0; j < i; j++)
      {
         if (mat(i, j) != 0 || mat(j, i) != 0)
         {
            flag = true;
         }
         if (flag)
         {
            asize++;
         }
      }
   }

   al.resize(asize); au.resize(asize);

   int s = 0;
   for (std::size_t i = 0; i < mat.Rows(); i++)
   {
      diag[i] = mat(i, i);
      ia[i] = s;
      bool flag = false;
      for (std::size_t j = 0; j < i; j++)
      {
         if (!flag && (mat(i, j) != 0 || mat(j, i) != 0))
         {
            flag = true;
         }
         if (flag)
         {
            au[s] = mat(i, j);
            al[s] = mat(j, i);
            s++;
         }
      }
   }
   ia[mat.Rows()] = s;

   type = ProfileMatrixType::ProfileOnly;
}

inline void ProfileMatrix::LUdecompose() {
   for (size_t i = 0; i < Size(); i++)
   {
      size_t i0 = ia[i];
      size_t i1 = ia[i + 1];
      size_t j = i - (ia[i + 1] - ia[i]);
      double bdi = 0;
      for (size_t k = ia[i]; k < ia[i + 1]; k++, j++)
      {
         size_t ki = ia[i];
         size_t kj = ia[j];
         size_t dif = k - ia[i] - ia[j + 1] + ia[j];
         if (dif < 0)
            kj -= dif;
         else
            ki += dif;
         double bal = 0;
         double bau = 0;
         for (ki; ki < k; ki++, kj++)
         {
            bal += al[ki] * au[kj];
            bau += au[ki] * al[kj];
         }
         al[k] = al[k] - bal;
         au[k] = au[k] - bau;
         au[k] = au[k] / diag[j];
         bdi += al[k] * au[k];
      }
      diag[i] -= bdi;
   }

   type = ProfileMatrixType::LUdecomposed;
}
