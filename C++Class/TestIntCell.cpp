/*
 * @Author: (xiao_huan) 1156492597@qq.com
 * @Date: 2024-04-02 21:47:58
 * @LastEditors: (xiao_huan) 1156492597@qq.com
 * @LastEditTime: 2024-04-03 10:14:30
 * @FilePath: \CPlusPlus Exercise\C++Class\TestIntCell.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <iostream>
#include "IntCell.h"
using namespace std;

int main(){
    IntCell *m;
    m = new IntCell(10);
    m->write(5);
    cout << m->read();
    delete m;
    return 0;
}
