#pragma once
#include "LU solver/headers/ProfileLU.h"
#include "json.hpp"
#include <cmath>
#include <functional>
#include <algorithm>
#include <iostream>
#include <format>

using json = nlohmann::json;

namespace Newtons {

   
   namespace Vec {
      inline double Scalar(const std::vector<double>& l, const std::vector<double>& r);

      inline double Norm(const std::vector<double>& vec);

      // ans = left + coef * right
      inline void AddVec(
         const std::vector<double>& left,
         double coef,
         const std::vector<double>& right,
         std::vector<double>& ans) 
      {
         for (size_t i = 0; i < right.size(); i++)
         {
            ans[i] = left[i] + coef * right[i];
         }
      }
   }


   class NewtonsSolver {
   public:

      struct TraceElement {
         int iterationNum{};
         std::vector<double> prevX;
         std::vector<double> X;
         std::vector<double> dX;
         double eps{};
         double prevEps{};

         TraceElement() noexcept {}

         TraceElement(const TraceElement& elem) noexcept {
            iterationNum = elem.iterationNum;
            prevX = (elem.prevX);
            X = (elem.X);
            dX = (elem.dX);
            eps = elem.eps;
            prevEps = elem.prevEps;
         }

         TraceElement(TraceElement&& elem) noexcept {
            iterationNum = elem.iterationNum;
            prevX = std::move(elem.prevX);
            X = std::move(elem.X);
            dX = std::move(elem.dX);
            eps = elem.eps;
            prevEps = elem.prevEps;
         }
      };

      class TraceVector {
      private:
         std::vector<TraceElement> _traceVec;

      public:
         TraceElement& operator[](size_t ind) {
            return _traceVec[ind];
         }

         void Push(TraceElement&& elem) {
            _traceVec.push_back(elem);
         }
         void Push(TraceElement& elem) {
            _traceVec.push_back(elem);
         }

         void Clear() {
            _traceVec.clear();
         }

         json ToJson() {
            json js;
            std::vector<json> semi_jsons;
            semi_jsons.resize(_traceVec.size());
            for (size_t i = 0; i < _traceVec.size(); i++)
            {
               semi_jsons[i]["iterNum"] = _traceVec[i].iterationNum;
               semi_jsons[i]["prevX"] = _traceVec[i].prevX;
               semi_jsons[i]["X"] = _traceVec[i].X;
               semi_jsons[i]["dX"] = _traceVec[i].dX;
               semi_jsons[i]["eps"] = _traceVec[i].eps;
               semi_jsons[i]["prevEps"] = _traceVec[i].prevEps;
            }
            js = semi_jsons;

            return js;
         }
      };

   private:

      // Переменные для матриц и векторов, используемых в солвере
      ProfileMatrix _profMat;
      Matrix _mat;
      std::vector<double> _F;
      std::vector<double> _x;
      std::vector<double> _dx;
      std::vector<double> _dx_trim;

      // Переменные для хранения фунцкий и их дифференциалов
      std::function<double(std::size_t, const std::vector<double>&)> _functions;
      std::function<double(std::size_t, const std::size_t, std::vector<double>&)> _differentials;

      // Переменные для хранения количества функций и переменных
      size_t _funcCount;
      size_t _varCount;

      // Вектор для формирования масок
      std::vector<std::pair<size_t, double>> _pairVec;

      // Для обозначения статуса маски
      enum class MaskType {
         None,
         MoreVars,
         MoreFuncs
      };
      MaskType _maskType = MaskType::None;

      // Маска для исключения лишних переменных или функций
      std::vector<bool> _mask;

      // Указатель на массив для трассировки метода (получение результата вычислений на каждом шагу)
      TraceVector* _traceVector;

   public:

      /// <summary>
      /// Инициализатор солвера методом Ньютона
      /// </summary>
      /// <param name="variableCount"> - количество переменных в системе</param>
      /// <param name="funcCount"> - количество функций в системе</param>
      /// <param name="functions"> - функция, в которой заданы все функции F_j системы</param>
      /// <param name="differentials"> - функция, в которой заданы все дифференциалы системы. </param>
      NewtonsSolver(
         size_t variableCount,
         size_t funcCount,
         std::function<double(size_t, const std::vector<double>&)> functions,
         std::function<double(size_t, size_t, const std::vector<double>&)> differentials)
      {
         _varCount = variableCount;
         _funcCount = funcCount;

         _x.resize(variableCount);
         _dx.resize(variableCount);

         _functions = functions;
         _differentials = differentials;

         size_t minSize = std::min(variableCount, funcCount);
         _mat.resize(minSize, minSize);
         _F.resize(minSize);

         if (variableCount != funcCount)
         {
            if (variableCount > funcCount)
            {
               _dx_trim.resize(funcCount);
            }

            _mask.resize(std::max(variableCount, funcCount));
            _pairVec.resize(std::max(variableCount, funcCount));
         }
      }


   private:

      // Определяет тип маски и получает маску, меняет _mask и _maskType
      inline void _GetMask();

      // Находит матрицу Якоби с учётом маски и записывает её в _mat
      void _GetJacobi();

      // Находит вектор функций F для решения системы
      void _GetF();

      // Находит норму вектора из значений фунций F в точке x
      double _GetNormF(const std::vector<double>& x) {
         double res = 0;
         double t;
         for (size_t i = 0; i < _funcCount; i++)
         {
            t = _functions(i, x);
            res += t * t;
         }
         return std::sqrt(res);
      }

   public:

      // Минимальное значение невязки вектора решения
      double minEps = 1e-5;

      // Максимальное число итераций
      int maxIter = 100;

      // Критический коэффициент, после которого метод завершается с ошибкой сходимости
      // По умолчанию равен числу, соответствующему 6 дроблениям коэффициента на 2
      double criticalCoef = 1 / (1 << 6);

      // Метод для решения системы нелинейных уравнений
      // - init_x - начальное приближение, в том числе итоговое решение
      // - eps - полученная невязка решения
      // Возможный возврат:
      // -1 - ошибка сходимости (при любом допустимо возможном шаге невязка возрастает)
      // -2 - выход по превышению числа итераций
      // Положительное число - число итераций сходимости метода
      int Solve(std::vector<double>& init_x, double& eps, const bool debugOutput = false);

      class TraceVector;

      void EnableTracing(TraceVector& traceVector) {
         _traceVector = &traceVector;
         traceVector.Clear();
      }
   };


}