#include <iostream>
#include <fstream>
#include <vector>


using namespace std;

/// <summary>
/// ��������� ������� F_i
/// </summary>
/// <param name="funcNum"> - ����� ������� ��� ������;</param>
/// <param name="vars"> - ��������� �������</param>
/// <returns>�������� ������� � ����� [vars]</returns>
double F(int funcNum, const vector<double>& vars) {
   switch (funcNum)
   {
      //case 0: return vars[0] + vars[1];
      //case 1: return vars[1];

      default: throw runtime_error("������������ ����� ������� F");
   }
}

/// <summary>
/// ��������� ����������� ������� F_i
/// </summary>
/// <param name="funcNum"> - ����� ������� ��� ������ �����������;</param>
/// <param name="varNum"> - ����� ���������, �� �������� ������ �����������;</param>
/// <param name="vars"> - ��������� �������</param>
/// <returns>�������� ������� � ����� [vars]</returns>
double div_F(int funcNum, int varNum, const vector<double>& vars) {
   switch (funcNum)
   {
      //case 0:
      //{
      //   switch (varNum)
      //   {
      //      case 0:
      //         return 1;
      //      case 1:
      //         return 1;

      //      default: throw runtime_error("������������ ����� ��������� �������");
      //   }
      //}


      //case 1:
      //{
      //   switch (varNum)
      //   {
      //      case 0:
      //         return 0;
      //      case 1:
      //         return 1;

      //      default: throw runtime_error("������������ ����� ��������� �������");
      //   }
      //}

      default: throw runtime_error("������������ ����� ������� F");
   }
}

int main() {


}