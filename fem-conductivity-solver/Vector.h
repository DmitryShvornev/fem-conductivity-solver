#pragma once
#include <math.h>

class Vector
{
public:
	Vector(int dimension);
	Vector(const Vector &vector);
	int getSize() const { return size; };
	void addElement(int position, double element);
	double findNorm() const;

	Vector& operator =(const Vector &vector);
	Vector operator *(const double &number) const;
	double operator *(const Vector &vector) const;
	Vector& operator +=(const  Vector &vector);
	Vector& operator -=(const  Vector &vector);
	Vector operator +(const  Vector &vector) const;
	Vector operator -(const  Vector &vector) const;
	double& operator [](int index) const { return elements[index]; };

	~Vector() { delete[] elements; };

private:
	int size;
	double* elements;
};