#pragma once

#include "Matrix.h"
#include <stdexcept>

class ProfileMatrix {
public:

   enum class ProfileMatrixType
   {
      Empty,
      ProfileOnly,
      LUdecomposed
   };


public:

   std::vector<double> diag;
   std::vector<std::size_t> ia;
   std::vector<double> al;
   std::vector<double> au;

   ProfileMatrixType type = ProfileMatrixType::Empty;


public:

   ProfileMatrix() {}

   ProfileMatrix(std::size_t diagSize, std::size_t aSize)
      : diag(diagSize), ia(diagSize+1), al(aSize), au(aSize) { }


public:

   // Return size of matrix (number of dimention of matrix, size of diag)
   inline std::size_t Size(void) const {
      return diag.size();
   }
   // Return number of elements in lower or upper triangles of matrix (size of al)
   inline std::size_t Asize(void) const {
      return al.size();
   }

   bool isEmpty() const { return type == ProfileMatrixType::Empty; }
   bool isLU() const { return type == ProfileMatrixType::LUdecomposed; }

   void MakeFromMatrix(const Matrix& mat);

   void LUdecompose();
};