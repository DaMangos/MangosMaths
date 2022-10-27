# MangoMath
mgo::Matrix is a container that encapsulates fixed M x N size Matrix.

This container acts very similar to a std::array holding M * N asthmatic types and has same semantics as a struct holding a C-style array. However When (M, N) equal (2, 1), (3, 1), or (4, 1) This container has same semantics as a struct holding 2, 3, or 4 asthmatic types. As this container is ment to mimic a M x N matrix, either M or N must be geater than 1 and M and N must be greater than 0 so the smallest matrix one can create is a 2 element matrix i.e. 2 x 1 or 1 x 2.

To access common typedef of the matrices then define MGO_USE_MATRIX_TYPES and to access series of mathatical constent macros then define MGO_USE_MATHS_DEFINES. 

Element access
--------------

  At(i, j)                                      :Rreturns the element at the i, j location with bounds cheacking. O(1)

  operator(i, j)                                :Returns the element at the i, j location without bounds cheacking. O(1)

  operator[k]                                   :Returns the element at the k location in the array without bounds cheacking. O(1)


Iterators
---------

  Begin()                                       :Returns an iterator to the beginning. O(1)
  
  End()                                         :Returns an iterator to the End. O(1)
  
  ReverseBegin()                                :Returns a reverse iterator to the beginning. O(1)
  
  ReverseEnd()                                  :Returns a reverse iterator to the end. O(1)
  
  ConstBegin()                                  :Returns an const iterator to the beginning. O(1)
  
  ConstEnd()                                    :Returns an const iterator to the End. O(1)
    
  ConstReverseBegin()                           :Returns a const reverse iterator to the beginning. O(1)
  
  ConstReverseEnd()                             :Returns a const reverse iterator to the end. O(1)
  
Operations
----------

  Fill(value)                                   :Fill the container with specified value. O(N * M)
  
  Swap(x, y, space)                             :Swaps rows or colums x and y. O(M) or O(N) respectively.
  (requires M > 1 && N > 1)
  
  Transpose()                                   :Transposes the matrix. O(N * M). 
  
  RowEchelon()                                  :Upper triangulates the matrix. O(M * N * (1.5 * M - 0.5))
  (requires M > 1 && N > 1 && f:oating point)
  
  ReducedRowEchelon()                           :Preforms Gaussian elimination on the matrix. O(M * N * (2M - 1))
  (requires M > 1 && N > 1 && f:oating point)

  Det()                                         :Finds the determinant of the matrix. O(M^2 * (1.5 * M - 0.5))
  (requires M == N && M > 1)

  Inverse()                                     :Finds the Inverse of the matrix. O(M^2 * (4 * M - 2))
  (requires M == N && M > 1)

  Trace()                                       :Finds the Trace of the matrix. O(M)
  (requires M == N && M > 1)
  
  Dot()                                         :Finds the dot product on two vectors. O(M)
  (requires N = 1)
  
  LengthSquared()                               :Finds the dot product on two vectors. O(M)
  (requires N = 1)
  
  Length()                                      :Finds the euclidean distances. O(M)
  (requires N = 1 && floating point)
  
  Angle(vector)                                 :Finds the angle between two vectors. O(3M)
  (requires N = 1 && floating point)
  
  Normalize()                                   :Returns the normalize vector. O(2M)
  (requires N = 1 && floating point)

  Theta()                                       :Returns the polar angle of the vector O(M)
  (requires N = 1 && N = 3 && floating point)
    
  Phi()                                         :Returns the azimuth angle of the vector O(M)
  (requires N = 1 && M < 4 && floating point)

  Cross(vector)                                 :Finds the dot product on two vectors. O(M * (M -1))
  (requires N = 1 && N = 3 && floating point)

  



