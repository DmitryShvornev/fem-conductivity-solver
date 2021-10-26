#pragma once

#include "Vector.h"

class Matrix
{
public:
	Matrix(int dimension1, int dimension2);
	Matrix(const Matrix &matrix);
	double** getElements() const { return elements; };
	int getRows() const { return rows; };
	int getColumns() const { return columns; };
	void setElement(int row, int column, double element);
	void addElement(int row, int column, double element);

	Matrix transpose() const;
	
	Matrix& operator =(const Matrix &matrix);
	Matrix operator *(const Matrix &matrix);
	Vector operator *(const Vector &vector) const;
	Matrix operator *(const double &number);

	~Matrix()
	{
		for (int i = 0; i < rows; i++)
			delete[] elements[i];
		delete[] elements;
	};

private:
	int rows;
	int columns;
	double** elements;
};