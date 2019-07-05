#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

struct Element {
    int row;
    int col;
    double value;
};

struct SparseMatrix {
    vector<Element> element;
    int nrows, ncols, nvalues;
    double rowsum, colsum;
};

SparseMatrix read(string mtx) {
    ifstream fin(mtx.c_str());
    SparseMatrix mat;
    
    int nrows, ncols, nvalues;
    double rowsum, colsum;

    // ignore header
    while (fin.peek() == '%') fin.ignore(2048, '\n');

    // get the number of rows, columns, and values
    fin >> nrows >> ncols >> nvalues;
    mat.nrows = nrows;
    mat.ncols = ncols;
    mat.nvalues = nvalues;
    cout << "# rows: " << mat.ncols << " # cols: " << mat.ncols << " # values: "<< mat.nvalues << endl;

    // populate the matrix
    for (int i=0; i < nvalues; i++){
        Element el;
        fin >> el.row >> el.col >> el.value;
        mat.element.push_back(el);

        // perform row and colsum


    }

    fin.close();


    return mat;
}

int main(int argc, char **argv){ 
    SparseMatrix mtx = read(argv[1]);
//    for (int i=0;i<20;i++){
//        cout << mtx.element[i].row << ", " << mtx.element[i].col << ", " << mtx.element[i].value << "\n";
//    }


    cout<<"You have entered: "<<argc-1<<" arguments."<<endl;
    for (int i=0; i<argc; i++){
        cout<<"Argument "<<i<<" = "<<argv[i]<<endl;
    }
    return 0;
}
