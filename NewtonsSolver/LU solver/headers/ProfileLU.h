#pragma once

#include "ProfileMatrix.h"

namespace LU {
   class ProfileSolver {
   private:

      ProfileSolver() {}

      static void _Direct(const ProfileMatrix& mat, std::vector<double>& x);

      static void _Reverse(const ProfileMatrix& mat, std::vector<double>& x);


   public:
      
      static void Solve(const ProfileMatrix& mat, std::vector<double>& x, const std::vector<double>& F) {
         if (!mat.isLU()){
            throw std::runtime_error("Profile matrix is not LU decomposed, that ProfileSolver needs.");
         }
         x = F;
         _Direct(mat, x);
         _Reverse(mat, x);
      }
   };
}