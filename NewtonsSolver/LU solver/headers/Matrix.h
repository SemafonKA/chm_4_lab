#pragma once

#include <vector>

class Matrix {
public:

   std::vector<std::vector<double>> elems;


public:
   Matrix() {}

   Matrix(size_t cols, size_t rows) {
      elems.resize(cols);
      for (size_t i = 0; i < elems.size(); i++)
      {
         elems[i].resize(rows);
      }
   }

   Matrix(std::vector<std::vector<double>>&& initMat) {
      elems = std::move(initMat);
   }

   Matrix(const std::vector<std::vector<double>>& initMat) {
      elems.resize(initMat.size());
      for (size_t i = 0; i < elems.size(); i++)
      {
         elems[i].resize(initMat.size());
         for (size_t k = 0; k < elems[i].size(); k++)
         {
            elems[i][k] = initMat[i][k];
         }
      }
   }

   Matrix(const Matrix& initMat) : Matrix(initMat.elems) {}


public:

   inline std::size_t Cols(void) const {
      if (Rows() == 0) return 0;
      return elems[0].size();
   }
   inline std::size_t Rows(void) const {
      return elems.size();
   }

   double& operator() (std::size_t x, std::size_t y) {
      return elems[x][y];
   }

   double operator() (std::size_t x, std::size_t y) const {
      return elems[x][y];
   }

};
