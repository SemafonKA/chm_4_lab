#pragma once

#include "ProfileMatrix.h"

namespace LU {
   class ProfileSolver {
   private:

      ProfileSolver() {}

      static void _Direct(const ProfileMatrix& mat, std::vector<double>& x, const std::vector<double>& F);

      static void _Reverse(const ProfileMatrix& mat, std::vector<double>& x, const std::vector<double>& F);


   public:
      
      static void Solve(const ProfileMatrix& mat, std::vector<double>& x, const std::vector<double>& F) {
         if (!mat.isLU()){
            throw std::runtime_error("Profile matrix is not LU decomposed, that ProfileSolver needs.");
         }
         _Direct(mat, x, F);
         _Reverse(mat, x, F);
      }
   };
}