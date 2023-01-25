#include "NewtonsSolver.h"

namespace Newtons {


   void NewtonsSolver::_GetMask() {
      if (_varCount == _funcCount)
      {
         _maskType = MaskType::None;
      }
       else if (_varCount < _funcCount)
       {
           size_t delta = _funcCount - _varCount;

           // Найдём для каждой функции их абсолютное значение в точке х
           for (size_t i = 0; i < _funcCount; i++)
           {
               _pairVec[i].first = i;
               _pairVec[i].second = std::abs(_functions(i, _x));
           }

           // Сортируем эту последовательность по возрастанию значений
           std::sort(_pairVec.begin(), _pairVec.end(), [](const std::pair<size_t, double>& l, const std::pair<size_t, double>& r)
               {
                   return l.second < r.second;
               }
           );

           // Исключаем индексы первых delta элементов из расчётов (помечаем их нулями в маске)
           size_t i;
           for (i = 0; i < delta; i++)
           {
               _mask[_pairVec[i].first] = false;
           }

           // Остальные помечаем единицами
           for (i; i < _funcCount; i++)
           {
               _mask[_pairVec[i].first] = true;
           }

           _maskType = MaskType::MoreFuncs;
       }

       else // _varCount > _funcCount
       {
           size_t delta = _varCount - _funcCount;

           // Найдём для каждой переменной максимальное абсолютное значение среди прозводных функций
           for (size_t i = 0; i < _varCount; i++)
           {
               _pairVec[i].first = i;

               // Перебираем производные всех функций, находим максимум для i переменной
               double a = 0;
               for (size_t k = 0; k < _funcCount; k++)
               {
                   a = std::max(a, std::abs(_differentials(k, i, _x)));
               }
               _pairVec[i].second = a;
           }

           // Сортируем эту последовательность по возрастанию значений
           std::sort(_pairVec.begin(), _pairVec.end(), [](const std::pair<size_t, double>& l, const std::pair<size_t, double>& r)
               {
                   return l.second < r.second;
               }
           );

           // Исключаем индексы первых delta элементов из расчётов (помечаем их нулями в маске)
           size_t i;
           for (i = 0; i < delta; i++)
           {
               _mask[_pairVec[i].first] = false;
           }

           // Остальные помечаем единицами
           for (i; i < _varCount; i++)
           {
               _mask[_pairVec[i].first] = true;
           }

           _maskType = MaskType::MoreVars;
       }
   }

   void NewtonsSolver::_GetJacobi() {
       switch (_maskType)
       {
           case MaskType::None:
           {
               for (size_t func = 0; func < _funcCount; func++)
               {
                   for (size_t var = 0; var < _varCount; var++)
                   {
                       _mat(func, var) = _differentials(func, var, _x);
                   }
               }
           }
           break;

           case MaskType::MoreVars:
           {
               for (size_t func = 0; func < _funcCount; func++)
               {
                   for (size_t var = 0, dvar = 0; var < _varCount; var++)
                   {
                       // Считаем только для тех переменных, которые не исключены маской из расчётов
                       if (_mask[var])
                       {
                           _mat(func, dvar) = _differentials(func, var, _x);
                           dvar++;
                       }
                   }
               }
           }
           break;

           case MaskType::MoreFuncs:
           {
               for (size_t func = 0, dfunc = 0; func < _funcCount; func++)
               {
                   // Считаем только те функции, которые не исключены маской из расчётов
                   if (_mask[func])
                   {
                       for (size_t var = 0; var < _varCount; var++)
                       {
                           _mat(dfunc, var) = _differentials(func, var, _x);
                       }
                       dfunc++;
                   }
               }
           }
           break;
       }
   }

   void NewtonsSolver::_GetF() {
      if (_maskType == MaskType::MoreFuncs)
      {
         for (size_t i = 0, k = 0; i < _funcCount; i++)
         {
            // Считаем только те функции, которые не исключены из расчётов маской
            if (_mask[i])
            {
               _F[k] = -_functions(i, _x);
               k++;
            }
         }
      }
      else
      {
         for (size_t i = 0; i < _funcCount; i++)
         {
            _F[i] = -_functions(i, _x);
         }
      }
   }

   // Метод для решения системы нелинейных уравнений
   // - init_x - начальное приближение, в том числе итоговое решение
   // - eps - полученная невязка решения
   // Возможный возврат:
   // - [-1] - ошибка сходимости (при любом допустимо возможном шаге невязка возрастает)
   // - [-2] - выход по превышению числа итераций
   // - [-3] - ошибка сходимости (метод не может иметь направления движения)
   // - [Положительное число] - число итераций сходимости метода

   int NewtonsSolver::Solve(std::vector<double>& init_x, double& eps, const bool debugOutput) {
      std::swap(init_x, _x);
      eps = _GetNormF(_x);


      int it;
      for (it = 1; it <= maxIter && eps > minEps; it++)
      {
         _GetMask();
         _GetJacobi();
         _profMat.MakeFromMatrix(_mat);
         _profMat.LUdecompose();
         _GetF();

         if (_dx_trim.size() != 0)
         {
            LU::ProfileSolver::Solve(_profMat, _dx_trim, _F);

            // Переписываем обрезанный вектор _dx_trim в полноценный _dx
            for (size_t i = 0, k = 0; i < _varCount; i++)
            {
               if (_mask[i])
               {
                  _dx[i] = _dx_trim[k];
                  k++;
               }
               else
               {
                  _dx[i] = 0;
               }
            }
         }
         else
         {
            LU::ProfileSolver::Solve(_profMat, _dx, _F);
         }

         for (auto& el : _dx)
         {
            if (std::abs(el) == std::numeric_limits<double>::infinity())
            {
               if (debugOutput)
               {
                  std::cout << "Выход по ошибке сходимости: методу некуда идти.\nПопробуйте сместить начальную точку в сторону.\n\n";
               }
               if (_traceVector)
               {
                  TraceElement elem;
                  elem.iterationNum = it;
                  elem.prevX = _x;
                  elem.X = _x;
                  elem.dX = _dx;
                  elem.prevEps = eps;
                  elem.eps = eps;

                  _traceVector->Push(std::move(elem));
               }
               std::swap(_x, init_x);

               return -3;
            }
         }

         double coef = 2;
         double newEps = eps;
         while (eps <= newEps && coef > criticalCoef)
         {
            coef /= 2;
            Vec::AddVec(_x, coef, _dx, init_x);
            newEps = _GetNormF(init_x);
         }

         if (coef <= criticalCoef)
         {
            init_x = _x;
            if (debugOutput)
            {
               std::cout << "Выход по ошибке сходимости - метод пришёл к оптимальной точке\n";
               std::cout << std::format("\tКоэф. \\beta:{:15.5f}\n\n", coef);
            }

            if (_traceVector)
            {
               TraceElement elem;
               elem.iterationNum = it;
               elem.prevX = _x;
               elem.X = init_x;
               elem.dX = _dx;
               elem.prevEps = eps;
               elem.eps = newEps;

               _traceVector->Push(std::move(elem));
            }
            return -1;
         }

         if (_traceVector)
         {
            TraceElement elem;
            elem.iterationNum = it;
            elem.prevX = _x;
            elem.X = init_x;
            elem.dX = _dx;
            elem.prevEps = eps;
            elem.eps = newEps;

            _traceVector->Push(std::move(elem));
         }

         _x = init_x;
         eps = newEps;

         if (debugOutput)
         {
            std::cout << std::format("Текущая итерация: {:>4}, текущая невязка: {:>10.2e}\n\tТекущая точка:", it, eps);
            for (size_t i = 0; i < _x.size(); i++)
            {
               std::cout << std::format("{0:15.5f}", _x[i]);
            }
            std::cout << "\n\tВектор сдвига:";
            for (size_t i = 0; i < _x.size(); i++)
            {
               std::cout << std::format("{0:15.5f}", _dx[i] * coef);
            }
            std::cout << std::format("\n\tКоэф. \\beta:{:15.5f}", coef);
            std::cout << "\n\n";
         }
      }

      if (eps > minEps)
      {
         if (debugOutput)
            std::cout << "Выход по превышению числа итераций\n\n";
         return -2;
      }

      if (debugOutput)
         std::cout << "Выход по достижению требуемой невязки решения\n\n";
      return it - 1;


   }
   


}

double Newtons::Vec::Scalar(const std::vector<double>& l, const std::vector<double>& r) {
   if (l.size() != r.size()) throw std::runtime_error("Size of vectors not same.");

   double res = 0.0;
   for (size_t i = 0; i < l.size(); i++)
   {
      res += l[i] * r[i];
   }

   return res;
}

double Newtons::Vec::Norm(const std::vector<double>& vec) {
   return std::sqrt(Scalar(vec, vec));
}
