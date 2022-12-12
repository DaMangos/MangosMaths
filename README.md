# MangoMath
mgo::matrix is a container that encapsulates fixed M x N size matrix.

This container acts very similar to a std::array holding M * N asthmatic types and has same semantics as a struct holding a C-style array. However When (M, N) equal (2, 1), (3, 1), or (4, 1) This container has same semantics as a struct holding 2, 3, or 4 asthmatic types. As this container is ment to mimic a M x N matrix, either M or N must be geater than 1 and M and N must be greater than 0 so the smallest matrix one can create is a 2 element matrix i.e. 2 x 1 or 1 x 2.

Element access
--------------

  at(i, j)                                      :Rreturns the element at the i, j location with bounds cheacking. O(1)

  operator(i, j)                                :Returns the element at the i, j location without bounds cheacking. O(1)

  operator[k]                                   :Returns the element at the k location in the array without bounds cheacking. O(1)


Iterators
---------

  begin()                                       :Returns an iterator to the beginning. O(1)
  
  end()                                         :Returns an iterator to the End. O(1)
  
  rbegin()                                      :Returns a reverse iterator to the beginning. O(1)
  
  rend()                                        :Returns a reverse iterator to the end. O(1)
  
  cbegin()                                      :Returns an const iterator to the beginning. O(1)
  
  cend()                                        :Returns an const iterator to the End. O(1)
    
  crbegin()                                     :Returns a const reverse iterator to the beginning. O(1)
  
  crend()                                       :Returns a const reverse iterator to the end. O(1)
  
Operations
----------

  data()                                        :Returns a pointer to the underling array 

  fill(value)                                   :Fill the container with specified value. O(N * M)
  
  swap(x, y, space)                             :Swaps rows or colums x and y. O(M) or O(N) respectively.
  (requires M > 1 && N > 1)
  
  transpose()                                   :Transposes the matrix. O(N * M). 
  
  row_echelon()                                 :Upper triangulates the matrix. O(M * N * (1.5 * M - 0.5))
  (requires M > 1 && N > 1 && f:oating point)
  
  reduced_row_echelon()                         :Preforms Gaussian elimination on the matrix. O(M * N * (2M - 1))
  (requires M > 1 && N > 1 && f:oating point)

  det()                                         :Finds the determinant of the matrix. O(M^2 * (1.5 * M - 0.5))
  (requires M == N && M > 1)

  inverse()                                     :Finds the Inverse of the matrix. O(M^2 * (4 * M - 2))
  (requires M == N && M > 1)

  trace()                                       :Finds the Trace of the matrix. O(M)
  (requires M == N && M > 1)
  
  dot()                                         :Finds the dot product on two vectors. O(M)
  (requires N = 1)
  
  length_squared()                              :Finds the dot product on two vectors. O(M)
  (requires N = 1)
  
  length()                                      :Finds the euclidean distances. O(M)
  (requires N = 1 && floating point)
  
  angle(vector)                                 :Finds the angle between two vectors. O(3M)
  (requires N = 1 && floating point)
  
  normalize()                                   :Returns the normalize vector. O(2M)
  (requires N = 1 && floating point)

  theta()                                       :Returns the polar angle of the vector O(M)
  (requires N = 1 && N = 3 && floating point)
    
  phi()                                         :Returns the azimuth angle of the vector O(M)
  (requires N = 1 && M < 4 && floating point)

  cross(vector)                                 :Finds the dot product on two vectors. O(M * (M -1))
  (requires N = 1 && N = 3 && floating point)
