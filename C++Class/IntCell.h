/*
 * @Author: (xiao_huan) 1156492597@qq.com
 * @Date: 2024-04-02 21:16:27
 * @LastEditors: (xiao_huan) 1156492597@qq.com
 * @LastEditTime: 2024-04-02 21:45:36
 * @FilePath: \CPlusPlus Exercise\C++Class\IntCell.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <iostream>

// Define a preprocessor to define symbols, usually as file names. Determine if it has been read.
#ifndef IntCell_H
#define IntCell_H

class IntCell
{
private:
    int storedValue;

public:
    // implicit conversion not allowed. every single variance constructor should be explicit.
    explicit IntCell(int initialValue = 0)
        : storedValue(initialValue){ }

    // accessor not mutator
    int read() const
        {return storedValue;}

    void write(int x)
        {storedValue = x;}
};

#endif

// int main(int argc, char* argv[]){
//     // IntCell cell = 10; // implicit conversion not allowed.
//     // IntCell cell(); // function declaration
//     IntCell cell1(10); // explicit constructor allowed.
//     IntCell cell;
//     std::cout << cell.read() << std::endl;
//     cell.write(10);
//     std::cout << cell.read();
// }
