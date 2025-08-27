#ifndef LINALG_H
#define LINALG_H

#include "device.cuh"
#include "interfaces.cuh"

template <typename T> class Vector {
    private:
        int length;
        Device device;
    
    public:
        T *data;

        Vector(constructorOverride()){};
        Vector(int length);
        Vector(Vector<T> &other);
        Vector<T> &operator=(Vector<T> &other);
        ~Vector();
        void setZero();
        int getLength();
        Device getDevice();
        T &operator[](int index);
        Vector<T> &operator+(Vector<T> &other);
        Vector<T> &operator-(Vector<T> &other);
        T operator*(Vector<T> &other);
        Vector<T> &operator*(T scalar);
        Vector<T> &operator/(T scalar);
        Vector<T> &operator+=(Vector<T> &other);
        Vector<T> &operator-=(Vector<T> &other);
};

template <typename T> Vector<T> &operator*(T scalar, Vector<T> &other);

class Matrix {
    protected:
        int rows, columns; 
        Device device;
    
    public:
        int getRows();
        int getCols();
};

template <typename T> class DenseMatrix: public Matrix {
    private:
        T *data;
    public:
        DenseMatrix(Device device, int rows, int cols);
        DenseMatrix(DenseMatrix<T> &other);
        DenseMatrix<T> &operator=(DenseMatrix<T> &other);
        ~DenseMatrix();
        void setZero();
        T *operator[](int index);
        DenseMatrix<T> &operator*(DenseMatrix<T> &other);
        DenseMatrix<T> &operator+(DenseMatrix<T> &other);
        DenseMatrix<T> &operator-(DenseMatrix<T> &other);
        Vector<T> &operator*(Vector<T> &other);
        DenseMatrix<T> &operator*(T scalar);
        typedef T type;
};

template <typename T> DenseMatrix<T> &operator*(T scalar, DenseMatrix<T> &other);

template <typename T> class SparseMatrix: public Matrix {
    private:
        T *values, *col_idx;
        int elementsPerRow;
    public:
        SparseMatrix(Device device, int rows, int cols, int elementsPerRow, T *values, T *col_idx);
        SparseMatrix(SparseMatrix<T> &other);
        SparseMatrix<T> &operator=(SparseMatrix<T> other);
        ~SparseMatrix();
        DenseMatrix<T> &toDense();
        Vector<T> &operator*(Vector<T> &other);
        SparseMatrix<T> &operator*(T scalar);
};

template <typename T> SparseMatrix<T> &operator*(T scalar, SparseMatrix<T> &other);

#endif