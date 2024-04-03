/*
 * @Author: (xiao_huan) 1156492597@qq.com
 * @Date: 2024-04-03 16:36:36
 * @LastEditors: (xiao_huan) 1156492597@qq.com
 * @LastEditTime: 2024-04-03 16:40:16
 * @FilePath: \CPlusPlus Exercise\CPlusPlus-Exercise\C++Class\TestMatrix.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include<iostream>
#include "Matrix.h"

using namespace std;

int main(int argc, char *argv[]){
    int n = 10;
    Matrix<double> matrix(n, n);
    for(int i = 0; i < n; i++){
        matrix[0][i] = 1;
    }
    cout << matrix[0][0];
}