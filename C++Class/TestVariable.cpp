/*
 * @Author: (xiao_huan) 1156492597@qq.com
 * @Date: 2024-04-03 10:28:03
 * @LastEditors: (xiao_huan) 1156492597@qq.com
 * @LastEditTime: 2024-04-03 10:36:47
 * @FilePath: \CPlusPlus Exercise\CPlusPlus-Exercise\C++Class\TestVariable.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <iostream>
#include <string>
#include <vector>

using namespace std;

const string & findMax(const vector<string> & arr){ // call by constant reference
    int maxIndex = 0; // attention: local variable
    for(int i = 0; i < arr.size(); i++){
        if(arr[maxIndex] < arr[i] ){
            maxIndex = i; // record the index, not copy the variable
        }
    }
    return arr[maxIndex];
}

int main(int argc, char *argv[]){
    vector<string> arr = {"abc", "dfe", "gcp"};
    string max = findMax(arr);
    cout << max;
}