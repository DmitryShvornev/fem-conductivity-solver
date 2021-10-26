#include "Matrix.h"

Matrix::Matrix(int dimension1, int dimension2) 
	: rows(dimension1), columns(dimension2)
{
	elements = new double*[rows];
	for (int i = 0; i < rows; i++)
	{
		elements[i] = new double[columns];
		for (int j = 0; j < columns; j++)
			elements[i][j] = 0.;
	}
}

Matrix::Matrix(const Matrix &matrix)
	: rows(matrix.rows), columns(matrix.columns) 
{
	elements = new double*[rows];
	for (int i = 0; i < rows; i++)
	{
		elements[i] = new double[columns];
		for (int j = 0; j < columns; j++)
			elements[i][j] = matrix.getElements()[i][j];
	}
}

void Matrix::setElement(int row, int column, double element)
{
	(*this).elements[row][column] = element;
}

void Matrix::addElement(int row, int column, double element)
{
	(*this).elements[row][column] += element;
}

Matrix Matrix::transpose() const
{
	Matrix result(columns, rows);
	for (int i = 0; i < rows; i++) 
		for (int j = 0; j < columns; j++) 
			result.setElement(j, i, (*this).getElements()[i][j]);
	return result;
}

Matrix& Matrix::operator =(const Matrix &matrix) 
{
	rows = matrix.getRows();
	columns = matrix.getColumns();
	for (int i = 0; i < rows; i++) 
		for (int j = 0; j < columns; j++)
			elements[i][j] = matrix.getElements()[i][j];
	return (*this);
}

Matrix Matrix::operator *(const Matrix &matrix)
{
	Matrix result(rows, matrix.columns);
	for (int i = 0; i < rows; i++) 
		for (int j = 0; j < matrix.columns; j++) 
			for (int k = 0; k < columns; k++) 
				result.addElement(
					i, j, (*this).getElements()[i][k] * matrix.getElements()[k][j]);
	return result;
}

Vector Matrix::operator *(const Vector &vector) const
{
	Vector result(rows);
	for (int i = 0; i < rows; i++) 
		for (int j = 0; j < columns; j++) 
			result[i] += (*this).getElements()[i][j] * vector[j];
	return result;
}

Matrix Matrix::operator *(const double &number) 
{
	Matrix result(*this);
	for (int i = 0; i < rows; i++) 
		for (int j = 0; j < columns; j++)
			result.setElement(i, j, (*this).getElements()[i][j] * number);
	return  result;
}