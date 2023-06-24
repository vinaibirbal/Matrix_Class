#include <iostream>
#include <math.h>
#include "SafeMatrix.h"

using namespace std;

void SafeMatrix::transpose() 
{
	if(_numRows >= 0)
	{
		float* tdata = new (std::nothrow) float [_numRows * _numCols];
		if (tdata == 0) 
		{
			cerr << "Uable to allocate space for 1000 Tasks"<< endl;
		}
		int temp;
		temp = _numRows;
		_numRows = _numCols;
		_numCols = temp;
		for(int i = 0; i < _numRows; i++)
		{
			for(int j = 0; j < _numCols; j ++)
			{
				tdata[i*_numCols + j] = _data[j*_numRows + i];
			}
		}
		delete []_data;
		_data = tdata;
	}
	else
	{
		_numRows = -2;
		cerr << "Unable to allocate data" << endl;
	}
}
bool SafeMatrix::appendRow(const int cols, const float data[]) 
{
	if(cols != _numCols)
	{
		return false;
	}
	float* tdata = new(std::nothrow) float[(_numRows + 1) * _numCols];
	if (tdata == 0) 
	{
		cerr << "Uable to allocate space for 1000 Tasks"<< endl;
		return false;
	}
	for(int i = 0; i < _numRows; i++)
	{
		for(int j = 0; j < _numCols; j ++)
		{
			tdata[i*_numCols + j] = _data[i*_numCols + j];
		}
	}
	for(int i = 0; i < _numCols ; i++)
	{
		tdata[_numRows*_numCols + i ] = data[i];
	}
	_numRows ++;
	delete []_data;
	_data = tdata;
	return true;
}

bool SafeMatrix::appendColumn(const int rows, const float data[]) 
{
	if(rows != _numRows)
	{
		return false;
	}
	transpose();
	appendRow(rows, data);
	transpose();
	return true;
}

Dimensions SafeMatrix::dimensions() const 
{
	Dimensions* Dim = new(std::nothrow) Dimensions;
	if (Dim == 0) 
	{
		cerr << "Uable to allocate space for 1000 Tasks"<< endl;
		Dim -> rows = NOT_A_MATRIX;
		Dim -> cols = -1;
		return *Dim ;
	}
	Dim->rows = _numRows;
	Dim->cols = _numCols;
	return *Dim;
}

  // Operators
float& SafeMatrix::     operator()(int i, int j){

	return MATRIX(*this,i,j);
}

SafeMatrix SafeMatrix::operator+(const SafeMatrix& m) const 
{
	if(m._numCols != _numCols || m._numRows != _numRows)
	{
		return SafeMatrix(-2,-2);
	}
	float* tdata = new(std::nothrow) float[_numRows * _numCols];
	if (tdata == 0) 
	{
		cerr << "Uable to allocate space for 1000 Tasks"<< endl;
		return SafeMatrix(-2,-2);
	}
	for(int i = 0; i < _numRows; i++)
	{
		for(int j = 0; j < _numCols; j ++)
		{
			tdata[i*_numCols + j] = _data[i*_numCols + j] + m._data[i*_numCols + j];
		}
	}
	SafeMatrix mout(_numRows,_numCols);
	mout._data = tdata;
	return mout;
}

SafeMatrix SafeMatrix::operator*(const SafeMatrix& m) const 
{
	if(m._numRows != _numCols)
	{
		return SafeMatrix(-2,-2);
	}
	SafeMatrix mout(_numRows, m._numCols, 0);
	for(int i = 0; i < _numRows; i++)
	{
		for(int j = 0; j < m._numCols; j ++)
		{
			for(int a = 0; a < _numCols; a++)
			{
				mout._data[i*_numCols + j] += _data[i*_numCols + a]* m._data[a*_numCols + j];
			}
		}
	}
	return mout;
}
void SafeMatrix::operator=(const SafeMatrix& m) 
{
	_numRows = m._numRows;
	_numCols = m._numCols;
	_dataSpaceAllocated = m._dataSpaceAllocated;
	_data = new (std::nothrow) float[m._dataSpaceAllocated];
	for(int i = 0; i < _numRows; i++)
	{
		for(int j = 0; j < _numCols; j ++)
		{
			_data[i*_numCols + j] = m._data[i*_numCols + j];
		}
	}
}
  // Constructors/Destructor
SafeMatrix::SafeMatrix()
{

	_numRows = 0;
	_numCols = 0;
	_dataSpaceAllocated = 0;
	
	_data = new (std::nothrow) float[_dataSpaceAllocated];
	
}

SafeMatrix::SafeMatrix(const SafeMatrix& m)
{
	
	_numRows = m._numRows;
	_numCols = m._numCols;
	_dataSpaceAllocated = m._dataSpaceAllocated;
	
	_data = new (std::nothrow) float[m._dataSpaceAllocated];
	
	for(int i = 0; i < _dataSpaceAllocated; i++)
	{
		_data[i] = m._data[i];
	}
}
	
SafeMatrix::SafeMatrix(const int rows, const int cols){                           
	if(rows < 0||cols < 0)
	{
		_numRows = -1;
		_numCols = -1;
		_data = NULL;
	}
	else
	{
		_numRows = rows;
		_numCols = cols;
		_dataSpaceAllocated = rows*cols;
		_data = new (std::nothrow) float[_dataSpaceAllocated];
	}
}

SafeMatrix::SafeMatrix(const int rows, const int cols, const float initVal){      
// Initialized matrix
	if(rows < 0||cols < 0)
	{
		_numRows = -1;
		_numCols = -1;
		_data = NULL;
	}
	else
	{
		_numRows = rows;
		_numCols = cols;
		_dataSpaceAllocated = rows*cols;
		
		_data = new (std::nothrow) float[_dataSpaceAllocated];
		
		for(int i = 0; i < _dataSpaceAllocated; i++)
		{
			_data[i] = initVal;
		}
	}
}

SafeMatrix::~SafeMatrix(){
	_numRows = MATRIX_DELETED;
	delete []_data;
}
  