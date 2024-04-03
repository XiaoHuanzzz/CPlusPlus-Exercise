/*
 * @Author: (xiao_huan) 1156492597@qq.com
 * @Date: 2024-04-03 16:10:03
 * @LastEditors: (xiao_huan) 1156492597@qq.com
 * @LastEditTime: 2024-04-03 16:39:18
 * @FilePath: \CPlusPlus Exercise\CPlusPlus-Exercise\C++Class\Matrix.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef MATRIX_H
#define MATIRX_H

#include <vector>

using namespace std;

template<class T>
class Matrix
{
private:
    vector<vector<T>> array;
public:
    Matrix(int rows, int cols): array(rows) {
        for(int i = 0; i < rows; i++)
            array[i].resize(cols);
    }

    const vector<T> & operator[](int row) const{
        return array[row];
    }

    vector<T> & operator[](int row){
        return array[row];
    }

    int numrows() const{
        return array.size();
    }

    int numcols() const{
        return numrows() ? array[0].size() : 0;
    }
};

template<class T>
void copyMatrix(Matrix<T> & from, Matrix<T> & to){
    for(int i = 0; i < to.numrows(); i++){
        to[i] = from[i];
    }
}

#endif