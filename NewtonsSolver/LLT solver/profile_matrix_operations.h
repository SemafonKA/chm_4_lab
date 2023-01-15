#pragma once

#include "predefined_types.h"
#include <iostream>

// Функция преобразует вектор y в вектор x (перезаписывает элементы!)
void GetVectorX(size_t _size, real_t* _matrixDiag, real_t* _matrixL, size_t* _matrixIA, real_t* _vecB) {

   // Сначала отдельно считаем последний элемент
   _vecB[_size - 1] /= _matrixDiag[_size - 1];

   for (size_t j = _size - 1; j > 0; j--) {
      // Вычитаем все элементы, содержащие x_j
      size_t jColFirstInd = j - _matrixIA[j + 1] + _matrixIA[j];
      real_t* jColPtr = _matrixL + (_matrixIA[j] - jColFirstInd);
      for (size_t i = jColFirstInd; i < j; i++) {
         _vecB[i] -= _vecB[j] * jColPtr[i];
      }

      // Высчитываем следующий x_j (x_(j-1))
      _vecB[j - 1] /= _matrixDiag[j - 1];
   }
}

// Функция преобразует вектор b в вектор y (перезаписывает элементы!)
void GetVectorY(size_t _size, real_t* _matrixDiag, real_t* _matrixL, size_t* _matrixIA, real_t* _vecB) {
   for (size_t i = 0; i < _size; i++) {
      accum_t accum = 0.0;
      size_t iStrFirstInd = i - _matrixIA[i + 1] + _matrixIA[i];
      real_t* iStrPtr = _matrixL + (_matrixIA[i] - iStrFirstInd);
      for (size_t j = iStrFirstInd; j < i; j++) {
         accum += iStrPtr[j] * _vecB[j];
      }
      _vecB[i] = (_vecB[i] - accum) / _matrixDiag[i];
   }
}

// Функция преобразует переданную матрицу в матрицу L (перезаписывает элементы!)
void GetMatrixL(size_t _size, real_t* _matrixDiag, real_t* _matrixAL, size_t* _matrixIA) {
   _matrixDiag[0] = std::sqrt(_matrixDiag[0]);     // Находим отдельно первый элемент
   accum_t accum = 0.0;
   for (size_t i = 1; i < _size; i++) {
      // Находим все элементы i-ой строки до диагонального элемента
      size_t iStrLen = _matrixIA[i + 1] - _matrixIA[i];
      size_t iStrFirstInd = i - iStrLen;
      real_t* iStrPtr = _matrixAL + (_matrixIA[i] - iStrFirstInd);   // Указатель на строку со смещением до "нулевой" позиции
      // Идём от первого элемента в i-ой строке
      for (size_t j = iStrFirstInd, jElemInd = _matrixIA[i]; j < i; j++, jElemInd++) {
         size_t jStrFirstInd = j - _matrixIA[j + 1] + _matrixIA[j];
         real_t* jStrPtr = _matrixAL + (_matrixIA[j] - jStrFirstInd);

         accum = 0.0;
         // k выбираем как максимальный начальный индекс для двух строк
         for (size_t k = std::max(iStrFirstInd, jStrFirstInd); k < j; k++) {
            accum += iStrPtr[k] * jStrPtr[k];
         }
         _matrixAL[jElemInd] = (_matrixAL[jElemInd] - accum) / _matrixDiag[j];
      }

      // Находим диагональный элемент на этой строке
      accum = 0.0;
      for (size_t k = _matrixIA[i]; k < _matrixIA[i + 1]; k++) {
         accum += _matrixAL[k] * _matrixAL[k];
      }
      // Если пошли мнимые корни
      if (accum > _matrixDiag[i]) {
         std::cout << "Во время вычисления матрицы возникли мнимые корни. Проверьте матрицу и повторите попытку." << std::endl;
         throw std::exception("Illegal matrix.");
      }

      if (IsAlmostEq(_matrixDiag[i], accum)) {
         std::cout << "Определитель матрицы равен нулю. Проверьте матрицу и повторите попытку" << std::endl;
         throw std::exception("Illegal matrix.");
      }
      else {
         _matrixDiag[i] = std::sqrt(_matrixDiag[i] - accum);  // Находим i-ый диагональный элемент
      }
   }
}
