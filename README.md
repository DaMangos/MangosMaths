# MangoMath
Light weight templated matrix and vector class.

To access the i j element without bounds checking, use the A(i, j)
notation. To access the i j element with bounds checking, use the
A.At(i, j) notation and finaly to access the underling array, use
the A[k] notation. If your matrix is of type vector with size <= 4
i.e. Matrix<Type, M, 1> where M = 1, 2, 3, 4, then you can access
the elements the same as if it were a M by 1 matrix, and by A.x and
A.y if M == 2, A.x, A.y, or A.z if M == 3, and A.x, A.y, A.z, or
A.w if M == 4. Depending on the size of the matrix and the data type,
some methords relevent to that matirx will become visable.  All
Matrix types iterator retriever methords that can be use in the
standard library's algorithms. These methords have the same name as
the the ones the standard library uses for it's containers. To accses
typedefs of common matrice types the define MGO_USE_MATRIX_DEFINES
and to use common math mathematical constant then define MGO_USE_MATHS_DEFINES.
