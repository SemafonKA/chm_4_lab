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

           // ����� ��� ������ ������� �� ���������� �������� � ����� �
           for (size_t i = 0; i < _funcCount; i++)
           {
               _pairVec[i].first = i;
               _pairVec[i].second = std::abs(_functions(i, _x));
           }

           // ��������� ��� ������������������ �� ����������� ��������
           std::sort(_pairVec.begin(), _pairVec.end(), [](const std::pair<size_t, double>& l, const std::pair<size_t, double>& r)
               {
                   return l.second < r.second;
               }
           );

           // ��������� ������� ������ delta ��������� �� �������� (�������� �� ������ � �����)
           size_t i;
           for (i = 0; i < delta; i++)
           {
               _mask[_pairVec[i].first] = false;
           }

           // ��������� �������� ���������
           for (i; i < _funcCount; i++)
           {
               _mask[_pairVec[i].first] = true;
           }

           _maskType = MaskType::MoreFuncs;
       }

       else // _varCount > _funcCount
       {
           size_t delta = _varCount - _funcCount;

           // ����� ��� ������ ���������� ������������ ���������� �������� ����� ���������� �������
           for (size_t i = 0; i < _varCount; i++)
           {
               _pairVec[i].first = i;

               // ���������� ����������� ���� �������, ������� �������� ��� i ����������
               double a = 0;
               for (size_t k = 0; k < _funcCount; k++)
               {
                   a = std::max(a, std::abs(_differentials(k, i, _x)));
               }
               _pairVec[i].second = a;
           }

           // ��������� ��� ������������������ �� ����������� ��������
           std::sort(_pairVec.begin(), _pairVec.end(), [](const std::pair<size_t, double>& l, const std::pair<size_t, double>& r)
               {
                   return l.second < r.second;
               }
           );

           // ��������� ������� ������ delta ��������� �� �������� (�������� �� ������ � �����)
           size_t i;
           for (i = 0; i < delta; i++)
           {
               _mask[_pairVec[i].first] = false;
           }

           // ��������� �������� ���������
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
                       // ������� ������ ��� ��� ����������, ������� �� ��������� ������ �� ��������
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
                   // ������� ������ �� �������, ������� �� ��������� ������ �� ��������
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
            // ������� ������ �� �������, ������� �� ��������� �� �������� ������
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

   // ����� ��� ������� ������� ���������� ���������
   // - init_x - ��������� �����������, � ��� ����� �������� �������
   // - eps - ���������� ������� �������
   // ��������� �������:
   // - [-1] - ������ ���������� (��� ����� ��������� ��������� ���� ������� ����������)
   // - [-2] - ����� �� ���������� ����� ��������
   // - [������������� �����] - ����� �������� ���������� ������

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

            // ������������ ���������� ������ _dx_trim � ����������� _dx
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

         double coef = 2;
         double newEps = eps;
         while (eps <= newEps && coef > criticalCoef)
         {
            coef /= 2;
            Vec::AddVec(_x, coef, _dx, init_x);
            newEps = _GetNormF(init_x);
         }
         if (coef < criticalCoef)
         {
            init_x = _x;
            if (debugOutput)
            {
               std::cout << "����� �� ������ ����������\n\n";
            }
            return -1;     // ������ - ����� �� ���������� ������
         }

         if (_traceVector)
         {
            TraceElement elem;
            elem.iterationNum = it;
            elem.prevX = _x;
            elem.X = init_x;
            elem.dX = _dx;
            elem.eps = newEps;
            elem.prevEps = eps;

            _traceVector->Push(std::move(elem));
         }

         _x = init_x;
         eps = newEps;

         if (debugOutput)
         {
            std::cout << std::format("������� ��������: {:>4}, ������� �������: {:>10.2e}\n\t������� �����:", it, eps);
            for (size_t i = 0; i < _x.size(); i++)
            {
               std::cout << std::format("{0:15.5f}", _x[i]);
            }
            std::cout << "\n\t������ ������:";
            for (size_t i = 0; i < _x.size(); i++)
            {
               std::cout << std::format("{0:15.5f}", _dx[i] * coef);
            }
            std::cout << "\n\n";
         }
      }

      if (eps > minEps)
      {
         if (debugOutput)
            std::cout << "����� �� ���������� ����� ��������\n\n";
         return -2;
      }

      if (debugOutput)
         std::cout << "����� �� ���������� ��������� ������� �������\n\n";
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
