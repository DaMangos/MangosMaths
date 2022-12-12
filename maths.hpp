#pragma once
#define _USE_MATHS_DEFINES
#define M_4PI       12.5663706143591729538505735331180115
#define M_3PI       9.42477796076937971538793014983850865
#define M_2PI       6.28318530717958647692528676655900577
#define M_3PI_2     4.71238898038468985769396507491925433
#define M_4PI_3     4.18879020478639098461685784437267051
#define M_3PI_4     2.35619449019234492884698253745962716
#define M_2PI_3     2.09439510239319549230842892218633526
#define M_Ef        2.71828182845904523536028747135266250f
#define M_LOG2Ef    1.44269504088896340735992468100189214f
#define M_LOG10Ef   0.43429448190325182765112891891660508f
#define M_LN2f      0.69314718055994530941723212145817657f
#define M_LN10f     2.30258509299404568401799145468436421f
#define M_4PIf      12.5663706143591729538505735331180115f
#define M_3PIf      9.42477796076937971538793014983850865f
#define M_2PIf      6.28318530717958647692528676655900577f
#define M_3PI_2f    4.71238898038468985769396507491925433f
#define M_4PI_3f    4.18879020478639098461685784437267051f
#define M_PIf       3.14159265358979323846264338327950288f
#define M_3PI_4f    2.35619449019234492884698253745962716f
#define M_2PI_3f    2.09439510239319549230842892218633526f
#define M_PI_2f     1.57079632679489661923132169163975144f
#define M_PI_4f     0.78539816339744830961566084581987572f
#define M_1_PIf     0.31830988618379067153776752674502872f
#define M_2_PIf     0.63661977236758134307553505349005745f
#define M_2_SQRTPIf 1.12837916709551257389615890312154517f
#define M_SQRT2f    1.41421356237309504880168872420969808f
#define M_SQRT1_2f  0.70710678118654752440084436210484904f
#include <type_traits>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <ostream>
namespace mgo
{
    template<typename T>
    concept arithmetic_type = std::is_arithmetic_v<T>;
    
#pragma mark - mgo::matrix<T, M, N>
    template<arithmetic_type T, std::size_t M, std::size_t N>
    requires ((M > 0 && N > 1) || (M > 1 && N > 0))
    class matrix final
    {
    private:
        typedef T                                           type;
        typedef matrix<T, M, N>                             matrixMxN;
        typedef matrix<T, N, M>                             matrixNxM;
        typedef matrix<T, M == 1 ? 2 : M, M == 1 ? 2 : M>   matrixMxM;
        typedef matrix<T, 1, N == 1 ? 2 : N>                matrix1xN;
        typedef matrix<T, M == 1 ? 2 : M, 1>                vectorM;
        typedef matrix<T, N == 1 ? 2 : N, 1>                vectorN;
        
        type elems_[M * N];
        
    public:
        enum init_type {columns, rows};
        
        constexpr matrix() noexcept
        :
        elems_{}
        {}
        
        constexpr explicit matrix(type value) noexcept
        requires (M > 1 && N > 1)
        :
        elems_{}
        {
            for (std::size_t k = 0; k < M && k < N; k++)
                (*this)(k, k) = value;
        }
        
        constexpr explicit matrix(type value) noexcept
        requires (M == 1 || N == 1)
        {
            this->fill(value);
        }
        
        template<arithmetic_type H>
        constexpr explicit matrix(const matrix<H, M, N>& value) noexcept
        {
            for (std::size_t k = 0; k < M * N; k++)
                (*this)[k] = static_cast<type>(value[k]);
        }
        
        constexpr matrix(std::initializer_list<type> list, init_type space = rows) noexcept
        requires (N > 1 && M > 1)
        :
        elems_{}
        {
            auto it = list.begin();
            if (space == rows)
                for (std::size_t i = 0; i < M; i++)
                    for (std::size_t j = 0; j < N && it != list.end(); j++, it++)
                        (*this)(i, j) = *it;
            if (space == columns)
                for (std::size_t j = 0; j < N; j++)
                    for (std::size_t i = 0; i < M && it != list.end(); i++, it++)
                        (*this)(i, j) = *it;
        }
        
        constexpr matrix(std::initializer_list<type> list) noexcept
        requires (M == 1 || N == 1)
        :
        elems_{}
        {
            auto it = list.begin();
            for (std::size_t k = 0; k < M * N && it != list.end(); k++, it++)
                (*this)[k] = *it;
        }
        
        constexpr matrix(const std::initializer_list<vectorM>& list, init_type space = columns) noexcept
        requires (M > 1 && N > 1)
        :
        elems_{}
        {
            auto it = list.begin();
            if (space == rows)
                for (std::size_t i = 0; i < M; i++, it++)
                    for (std::size_t j = 0; j < N && it != list.end(); j++)
                        if (j < M)
                            (*this)(i, j) = (*it)(j, 0);
            if (space == columns)
                for (std::size_t j = 0; j < N; j++, it++)
                    for (std::size_t i = 0; i < M && it != list.end(); i++)
                        (*this)(i, j) = (*it)(i, 0);
        }
        
        constexpr matrix(const std::initializer_list<matrix1xN>& list, init_type space = rows) noexcept
        requires (M > 1 && N > 1)
        :
        elems_{}
        {
            auto it = list.begin();
            if (space == rows)
                for (std::size_t i = 0; i < M; i++, it++)
                    for (std::size_t j = 0; j < N && it != list.end(); j++)
                        (*this)(i, j) = (*it)(0, j);
            if (space == columns)
                for (std::size_t j = 0; j < N; j++, it++)
                    for (std::size_t i = 0; i < M && it != list.end(); i++)
                        if (i < N)
                            (*this)(i, j) = (*it)(0, i);
        }
        
        constexpr type* begin() noexcept
        {
            return this->elems_;
        }
        
        constexpr const type* begin() const noexcept
        {
            return this->elems_;
        }
        
        constexpr type* end() noexcept
        {
            return this->begin() + M * N;
        }
        
        constexpr const type* end() const noexcept
        {
            return this->begin() + M * N;
        }
        
        constexpr std::reverse_iterator<type*>rbegin() noexcept
        {
            return std::reverse_iterator<type*>(this->end());
        }
        
        constexpr std::reverse_iterator<const type*> rbegin() const noexcept
        {
            return std::reverse_iterator<const type*>(this->end());
        }
        
        constexpr std::reverse_iterator<type*>rend() noexcept
        {
            return std::reverse_iterator<type*>(this->begin());
        }
        
        constexpr std::reverse_iterator<const type*> rend() const noexcept
        {
            return std::reverse_iterator<const type*>(this->begin());
        }
        
        constexpr const type* cbegin() const noexcept
        {
            return this->begin();
        }
        
        constexpr const type* cend() const noexcept
        {
            return this->end();
        }
        
        constexpr std::reverse_iterator<const type*> crbegin() const noexcept
        {
            return this->rbegin();
        }
        
        constexpr std::reverse_iterator<const type*> crend() const noexcept
        {
            return this->rend();
        }
        
        constexpr std::size_t size() const noexcept
        {
            return M * N;
        }
        
        constexpr std::size_t max_size() const noexcept
        {
            return M * N;
        }
        
        constexpr bool empty() const noexcept
        {
            return false;
        }
        
        constexpr type* data() noexcept
        {
            this->begin();
        }
        
        constexpr const type* data() const noexcept
        {
            this->begin();
        }
        
        constexpr type& at(std::size_t i, std::size_t j)
        {
            if (j + i * N >= M * N)
                throw std::out_of_range("matrix::at");
            return this->elems_[j + i * N];
        }
        
        constexpr const type& at(std::size_t i, std::size_t j) const
        {
            if (j + i * N >= M * N)
                throw std::out_of_range("matrix::at");
            return this->elems_[j + i * N];
        }
        
        constexpr type& operator()(std::size_t i, std::size_t j) noexcept
        {
            return this->elems_[j + i * N];
        }
        
        constexpr const type& operator()(std::size_t i, std::size_t j) const noexcept
        {
            return this->elems_[j + i * N];
        }
        
        constexpr type& operator[](std::size_t k) noexcept
        {
            return this->elems_[k];
        }
        
        constexpr const type& operator[](std::size_t k) const noexcept
        {
            return this->elems_[k];
        }
        
        template<std::size_t P>
        constexpr matrix<type, M, P> operator*(const matrix<type, N, P>& value) const noexcept
        {
            matrix<type, M, P> out;
            for (std::size_t i = 0; i < M; i++)
                for (std::size_t j = 0; j < N; j++)
                    for (std::size_t k = 0; k < P; k++)
                        out(i, k) += (*this)(i, j) * value(j, k);
            return out;
        }
        
        constexpr type operator*(const vectorN& value) const noexcept
        requires (M == 1)
        {
            vectorN out;
            for (std::size_t k = 0; k < N; k++)
                out += (*this)[k] * value[k];
            return out;
        }
        
        constexpr matrixMxN operator*(type value) const noexcept
        {
            matrixMxN out;
            for (std::size_t k = 0; k < M * N; k++)
                out = (*this)[k] * value;
            return out;
        }
        
        constexpr matrixMxN operator/(type value) const noexcept
        {
            matrixMxN out;
            for (std::size_t k = 0; k < M * N; k++)
                out = (*this)[k] / value;
            return out;
        }
        
        constexpr matrixMxN operator+(const matrixMxN& value) const noexcept
        {
            matrixMxN out;
            for (std::size_t k = 0; k < M * N; k++)
                out[k] = (*this)[k] + value[k];
            return out;
        }
        
        constexpr matrixMxN operator-(const matrixMxN& value) const noexcept
        {
            matrixMxN out;
            for (std::size_t k = 0; k < M * N; k++)
                out[k] = (*this)[k] - value[k];
            return out;
        }
        
        constexpr matrixMxM operator*=(const matrixMxM& value) noexcept
        requires (M == N)
        {
            *this = *this * value;
            return *this;
        }
        
        constexpr matrixMxM operator*=(type value) noexcept
        {
            *this = *this * value;
            return *this;
        }
        
        constexpr matrixMxM operator/=(type value) noexcept
        {
            *this = *this / value;
            return *this;
        }
        
        constexpr matrixMxN operator+=(const matrixMxN& value) noexcept
        {
            *this = *this + value;
            return *this;
        }
        
        constexpr matrixMxN operator-=(const matrixMxN& value) noexcept
        {
            *this = *this - value;
            return *this;
        }
        
        constexpr void fill(const type& value) noexcept
        {
            for (std::size_t k = 0; k < M * N; k++)
                (*this)[k] = value;
        }
        
        constexpr void swap(std::size_t x, std::size_t y, init_type space) noexcept
        requires (M > 1 && N > 1)
        {
            if (space == rows)
                for (std::size_t i = 0; i < N; i++)
                    std::swap(this->at(i, x), this->at(i, y));
            if (space == columns)
                for (std::size_t j = 0; j < M; j++)
                    std::swap(this->at(x, j), this->at(x, j));
        }
        
        constexpr matrixNxM transpose() const noexcept
        {
            matrixNxM out;
            for (std::size_t i = 0; i < M; i++)
                for (std::size_t j = 0; j < N; j++)
                    out(j, i) = (*this)(i, j);
            return out;
        }
        
        constexpr matrixMxN row_echelon() const noexcept
        requires (M > 1 && N > 1 && std::is_floating_point_v<type>)
        {
            matrixMxN out(*this);
            for (std::size_t k = 0; k < M; k++)
            {
                for (std::size_t i = k; out(k, k) == 0; i++)
                {
                    if (i == M)
                    {k++; break;}
                    
                    if (out(i, k) != 0)
                    {out.swap(k, i, rows); break;}
                }
                for (std::size_t i = k + 1; i < M; i++)
                {
                    type x = out(i, k) / out(k, k);
                    for (std::size_t j = 0; j < N; j++)
                        out(i, j) -= x * out(k, j);
                }
            }
            return out;
        }
        
        constexpr matrixMxN reduced_row_echelon() const noexcept
        requires (M > 1 && N > 1 && std::is_floating_point_v<type>)
        {
            matrixMxN out(*this);
            for (std::size_t k = 0; k < M; k++)
            {
                for (std::size_t i = k; out(k, k) == 0; i++)
                {
                    if (i == M)
                    {k++; break;}
                    
                    if (out(i, k) != 0)
                    {out.swap(k, i, rows); break;}
                }
                for (std::size_t i = 0; i < M && k < M; i++)
                    if (i == k)
                    {
                        type x = out(k, k);
                        for (std::size_t j = 0; j < N; j++)
                            out(i, j) /= x;
                    }
                    else
                    {
                        type x = out(i, k) / out(k, k);
                        for (std::size_t j = 0; j < N; j++)
                            out(i, j) -= x * out(k, j);
                    }
            }
            return out;
        }
        
        constexpr type det() const noexcept
        requires (M == N && M > 1 && N > 1)
        {
            matrixMxM in(*this);
            type out = 1;
            for (std::size_t k = 0; k < M; k++)
            {
                for (std::size_t i = k; in(k, k) == 0; i++)
                {
                    if (i == M)
                        return 0;
                    if (in(i, k) != 0)
                    {
                        in.swap(k, i, rows);
                        out *= -1;
                        break;
                    }
                }
                for (std::size_t i = k + 1; i < M; i++)
                {
                    type x = in(i, k) / in(k, k);
                    for (std::size_t j = 0; j < M; j++)
                        in(i, j) -= x * in(k, j);
                }
                out *= in(k, k);
            }
            return out;
        }
        
        constexpr matrixMxM Inverse() const noexcept
        requires (M == N && M > 1 && N > 1)
        {
            matrixMxM in(*this);
            matrixMxM out(1);
            for (std::size_t k = 0; k < M; k++)
            {
                for (std::size_t i = k; in(k, k) == 0; i++)
                {
                    if (i == M)
                        return matrixMxM(0);
                    if (in(i, k) != 0)
                    {
                        in.swap(k, i, rows);
                        out.swap(k, i, rows);
                        break;
                    }
                }
                for (std::size_t i = 0; i < M && k < M; i++)
                {
                    if (i == k)
                    {
                        type x = in(k, k);
                        for (std::size_t j = 0; j < M; j++)
                        {
                            in(i, j) /= x;
                            out(i, j) /= x;
                        }
                    }
                    else
                    {
                        type x = in(i, k) / in(k, k);
                        for (std::size_t j = 0; j < M; j++)
                        {
                            in(i, j) -= x * in(k, j);
                            out(i, j) -= x * out(k, j);
                        }
                    }
                }
            }
            return out;
        }
        
        constexpr type Trace() const noexcept
        requires (M == N && M > 1 && N > 1)
        {
            type out = 1;
            for (std::size_t k = 0; k < M; k++)
                out *= (*this)(k, k);
            return out;
        }
        
        constexpr type dot(const vectorM& value) const
        requires (N == 1)
        {
            type out;
            for (std::size_t k = 0; k < M; k++)
                out += (*this)[k] * value[k];
            return out;
        }
        
        constexpr type length_squared() const
        requires (N == 1)
        {
            return this->dot(*this);
        }
        
        constexpr type length() const
        requires (N == 1 && std::is_floating_point_v<type>)
        {
            return sqrt(this->length_squared());
        }
        
        constexpr type angle(const vectorM& value) const
        requires (N == 1 && std::is_floating_point_v<type>)
        {
            return acos(this->dot(value) / (this->length() * value.length));
        }
        
        constexpr vectorM normalize() const
        requires (N == 1 && std::is_floating_point_v<type>)
        {
            return *this / this->length();
        }
    };
    
#pragma mark - mgo::matrix<T, 2, 1>
    template<arithmetic_type T>
    class matrix<T, 2, 1> final
    {
    private:
        typedef T                               type;
        typedef matrix<T, 1, 2>                 matrix1x2;
        typedef matrix<T, 2, 1>                 vector2;
        
    public:
        enum init_type {spherical, cartesian};
        
        type x, y;
        
        constexpr matrix() noexcept
        :
        x(0), y(0)
        {}
        
        constexpr explicit matrix(type value) noexcept
        :
        x(value), y(value)
        {}
        
        template<arithmetic_type H>
        constexpr explicit matrix(const matrix<H, 2, 1>& value) noexcept
        :
        x(static_cast<type>(value.x)),
        y(static_cast<type>(value.y))
        {}
        
        constexpr matrix(std::initializer_list<type> list) noexcept
        :
        x(list.size() > 0 ? *(list.begin())     : 0),
        y(list.size() > 1 ? *(list.begin() + 1) : 0)
        {}
        
        constexpr matrix(type x, type y, init_type coords = cartesian) noexcept
        requires (std::is_floating_point_v<type>)
        :
        x(coords ? x : fabs(x) * cos(y)),
        y(coords ? y : fabs(x) * sin(y))
        {}
        
        template<std::size_t P>
        constexpr explicit matrix(const matrix<type, P, 1>& value) noexcept
        :
        x(P > 0 ? value[0] : 0),
        y(P > 1 ? value[1] : 0)
        {}
        
        constexpr type* begin() noexcept
        {
            return &this->x;
        }
        
        constexpr const type* begin() const noexcept
        {
            return &this->x;
        }
        
        constexpr type* end() noexcept
        {
            return this->begin() + 2;
        }
        
        constexpr const type* end() const noexcept
        {
            return this->begin() + 2;
        }
        
        constexpr std::reverse_iterator<type*>rbegin() noexcept
        {
            return std::reverse_iterator<type*>(this->end());
        }
        
        constexpr std::reverse_iterator<const type*> rbegin() const noexcept
        {
            return std::reverse_iterator<const type*>(this->end());
        }
        
        constexpr std::reverse_iterator<type*>rend() noexcept
        {
            return std::reverse_iterator<type*>(this->begin());
        }
        
        constexpr std::reverse_iterator<const type*> rend() const noexcept
        {
            return std::reverse_iterator<const type*>(this->begin());
        }
        
        constexpr const type* cbegin() const noexcept
        {
            return this->begin();
        }
        
        constexpr const type* cend() const noexcept
        {
            return this->end();
        }
        
        constexpr std::reverse_iterator<const type*> crbegin() const noexcept
        {
            return this->rbegin();
        }
        
        constexpr std::reverse_iterator<const type*> crend() const noexcept
        {
            return this->rend();
        }
        
        constexpr std::size_t size() const noexcept
        {
            return 2;
        }
        
        constexpr std::size_t max_size() const noexcept
        {
            return 2;
        }
        
        constexpr bool empty() const noexcept
        {
            return false;
        }
        
        constexpr type* data() noexcept
        {
            return this->begin();
        }
        
        constexpr const type* data() const noexcept
        {
            this->begin();
        }
        
        constexpr type& at(std::size_t i, std::size_t j)
        {
            if (j > 1 || i != 1)
                throw std::out_of_range("matrix::at");
            return i == 0 ? this->x : this->y;
        }
        
        constexpr const type& at(std::size_t i, std::size_t j) const
        {
            if (j > 1 || i != 1)
                throw std::out_of_range("matrix::at");
            return i == 0 ? this->x : this->y;
        }
        
        constexpr type& operator()(std::size_t i, std::size_t j) noexcept
        {
            return i == 0 ? this->x : this->y;
        }
        
        constexpr const type& operator()(std::size_t i, std::size_t j) const noexcept
        {
            return i == 0 ? this->x : this->y;
        }
        
        constexpr type& operator[](std::size_t k) noexcept
        {
            return k == 0 ? this->x : this->y;
        }
        
        constexpr const type& operator[](std::size_t k) const noexcept
        {
            return k == 0 ? this->x : this->y;
        }
        
        template<std::size_t P>
        constexpr matrix<type, 2, P> operator*(const matrix<type, 1, P>& value) const noexcept
        {
            matrix<type, 2, P> out;
            for (std::size_t i = 0; i < 2; i++)
                for (std::size_t k = 0; k < P; k++)
                    out(i, k) += (*this)(i, 0) * value(0, k);
            return out;
        }
        
        constexpr vector2 operator*(type value) const noexcept
        {
            return vector2({this->x * value, this->y * value});
        }
        
        constexpr vector2 operator/(type value) const noexcept
        {
            return vector2({this->x / value, this->y / value});
        }
        
        constexpr vector2 operator+(const vector2& value) const noexcept
        {
            return vector2({this->x + value.x, this->y + value.y});
        }
        
        constexpr vector2 operator-(const vector2& value) const noexcept
        {
            return vector2({this->x - value.x, this->y - value.y});
        }
        
        constexpr vector2 operator*=(type value) noexcept
        {
            *this = *this * value;
            return *this;
        }
        
        constexpr vector2 operator/=(type value) noexcept
        {
            *this = *this / value;
            return *this;
        }
        
        constexpr vector2 operator+=(const vector2& value) noexcept
        {
            *this = *this + value;
            return *this;
        }
        
        constexpr vector2 operator-=(const vector2& value) noexcept
        {
            *this = *this - value;
            return *this;
        }
        
        constexpr void fill(const type& value) noexcept
        {
            for (std::size_t k = 0; k < 2; k++)
                (*this)[k] = value;
        }
        
        constexpr matrix1x2 transpose() const noexcept
        {
            return matrix1x2({this->x, this->y});
        }
        
        constexpr type dot(const vector2& value) const
        {
            return this->x * value.x + this->y * value.y;
        }
        
        constexpr type length_squared() const
        {
            return this->dot(*this);
        }
        
        constexpr type length() const
        requires (std::is_floating_point_v<type>)
        {
            return sqrt(this->length_squared());
        }
        
        constexpr type angle(const vector2& value) const
        requires (std::is_floating_point_v<type>)
        {
            return acos(this->dot(value) / (this->length() * value.length));
        }
        
        constexpr vector2 normalize() const
        requires (std::is_floating_point_v<type>)
        {
            return *this / this->length();}
        
        constexpr auto phi() const
        requires (std::is_floating_point_v<type>)
        {
            return fmod(atan2(this->y, this->x), static_cast<type>(M_2PI));
        }
    };
    
#pragma mark - mgo::matrix<T, 3, 1>
    template<arithmetic_type T>
    class matrix<T, 3, 1> final
    {
    private:
        typedef T                               type;
        typedef matrix<T, 1, 3>                 matrix1x3;
        typedef matrix<T, 3, 1>                 vector3;
        
    public:
        enum init_type {spherical, cartesian};
        
        type x, y, z;
        
        constexpr matrix() noexcept
        :
        x(0), y(0), z(0)
        {}
        
        constexpr explicit matrix(type value) noexcept
        :
        x(value), y(value), z(value)
        {}
        
        template<arithmetic_type H>
        constexpr explicit matrix(const matrix<H, 3, 1>& value) noexcept
        :
        x(static_cast<type>(value.x)),
        y(static_cast<type>(value.y)),
        z(static_cast<type>(value.z))
        {}
        
        constexpr matrix(std::initializer_list<type> list) noexcept
        :
        x(list.size() > 0 ? *(list.begin())     : 0),
        y(list.size() > 1 ? *(list.begin() + 1) : 0),
        z(list.size() > 2 ? *(list.begin() + 2) : 0)
        {}
        
        constexpr matrix(type x, type y, type z, init_type coords = cartesian) noexcept
        requires (std::is_floating_point_v<type>)
        :
        x(coords ? x : fabs(x) * sin(y) * cos(z)),
        y(coords ? y : fabs(x) * sin(y) * sin(z)),
        z(coords ? z : fabs(x) * cos(y))
        {}
        
        template<std::size_t P>
        constexpr explicit matrix(const matrix<type, P, 1>& value) noexcept
        :
        x(P > 0 ? value[0] : 0),
        y(P > 1 ? value[1] : 0),
        z(P > 2 ? value[2] : 0)
        {}
        
        constexpr type* begin() noexcept
        {
            return &this->x;
        }
        
        constexpr const type* begin() const noexcept
        {
            return &this->x;
        }
        
        constexpr type* end() noexcept
        {
            return this->begin() + 3;
        }
        
        constexpr const type* end() const noexcept
        {
            return this->begin() + 3;
        }
        
        constexpr std::reverse_iterator<type*>rbegin() noexcept
        {
            return std::reverse_iterator<type*>(this->end());
        }
        
        constexpr std::reverse_iterator<const type*> rbegin() const noexcept
        {
            return std::reverse_iterator<const type*>(this->end());
        }
        
        constexpr std::reverse_iterator<type*>rend() noexcept
        {
            return std::reverse_iterator<type*>(this->begin());
        }
        
        constexpr std::reverse_iterator<const type*> rend() const noexcept
        {
            return std::reverse_iterator<const type*>(this->begin());
        }
        
        constexpr const type* cbegin() const noexcept
        {
            return this->begin();
        }
        
        constexpr const type* cend() const noexcept
        {
            return this->end();
        }
        
        constexpr std::reverse_iterator<const type*> crbegin() const noexcept
        {
            return this->rbegin();
        }
        
        constexpr std::reverse_iterator<const type*> crend() const noexcept
        {
            return this->rend();
        }
        
        constexpr std::size_t size() const noexcept
        {
            return 3;
        }
        
        constexpr std::size_t max_size() const noexcept
        {
            return 3;
        }
        
        constexpr bool empty() const noexcept
        {
            return false;
        }
        
        constexpr type* data() noexcept
        {
            return this->begin();
        }
        
        constexpr const type* data() const noexcept
        {
            this->begin();
        }
        
        constexpr type& at(std::size_t i, std::size_t j)
        {
            if (j > 2 || i != 1)
                throw std::out_of_range("matrix::at");
            return i == 0 ? this->x : (i == 1 ? this->y : this->z);
        }
        
        constexpr const type& at(std::size_t i, std::size_t j) const
        {
            if (j > 2 || i != 1)
                throw std::out_of_range("matrix::at");
            return i == 0 ? this->x : (i == 1 ? this->y : this->z);
        }
        
        constexpr type& operator()(std::size_t i, std::size_t j) noexcept
        {
            return i == 0 ? this->x : (i == 1 ? this->y : this->z);
        }
        
        constexpr const type& operator()(std::size_t i, std::size_t j) const noexcept
        {
            return i == 0 ? this->x : (i == 1 ? this->y : this->z);
        }
        
        constexpr type& operator[](std::size_t k) noexcept
        {
            return k == 0 ? this->x : (k == 1 ? this->y : this->z);
        }
        
        constexpr const type& operator[](std::size_t k) const noexcept
        {
            return k == 0 ? this->x : (k == 1 ? this->y : this->z);
        }
        
        template<std::size_t P>
        constexpr matrix<type, 3, P> operator*(const matrix<type, 1, P>& value) const noexcept
        {
            matrix<type, 3, P> out;
            for (std::size_t i = 0; i < 3; i++)
                for (std::size_t k = 0; k < P; k++)
                    out(i, k) += (*this)(i, 0) * value(0, k);
            return out;
        }
        
        constexpr vector3 operator*(type value) const noexcept
        {
            return vector3({this->x * value, this->y * value, this->z * value});
        }
        
        constexpr vector3 operator/(type value) const noexcept
        {
            return vector3({this->x / value, this->y / value, this->z / value});
        }
        
        constexpr vector3 operator+(const vector3& value) const noexcept
        {
            return vector3({this->x + value.x, this->y + value.y, this->z + value.z});
        }
        
        constexpr vector3 operator-(const vector3& value) const noexcept
        {
            return vector3({this->x - value.x, this->y - value.y, this->z - value.z});
        }
        
        constexpr vector3 operator*=(type value) noexcept
        {
            *this = *this * value;
            return *this;
        }
        
        constexpr vector3 operator/=(type value) noexcept
        {
            *this = *this / value;
            return *this;
        }
        
        constexpr vector3 operator+=(const vector3& value) noexcept
        {
            *this = *this + value;
            return *this;
        }
        
        constexpr vector3 operator-=(const vector3& value) noexcept
        {
            *this = *this - value;
            return *this;
        }
        
        constexpr void fill(const type& value) noexcept
        {
            for (std::size_t k = 0; k < 3; k++)
                (*this)[k] = value;
        }
        
        constexpr matrix1x3 transpose() const noexcept
        {
            return matrix1x4({this->x, this->y, this->z});
        }
        
        constexpr type dot(const vector3& value) const
        {
            return this->x * value.x + this->y * value.y + this->z * value.z;
        }
        
        constexpr type length_squared() const
        {
            return this->dot(*this);
        }
        
        constexpr type length() const
        requires (std::is_floating_point_v<type>)
        {
            return sqrt(this->length_squared());
        }
        
        constexpr type angle(const vector3& value) const
        requires (std::is_floating_point_v<type>)
        {
            return acos(this->dot(value) / (this->length() * value.length));
        }
        
        constexpr vector3 normalize() const
        requires (std::is_floating_point_v<type>)
        {
            return *this / this->length();
        }
        
        constexpr auto theta() const
        requires (std::is_floating_point_v<type>)
        {
            return acos(this->z / this->length());
        }
        
        constexpr auto phi() const
        requires (std::is_floating_point_v<type>)
        {
            return fmod(atan2(this->y, this->x), static_cast<type>(M_2PI));
        }
        
        constexpr vector3 cross(const vector3& value) const
        {
            return vector3((this->y * value.z) - (this->z * value.y),
                           (this->z * value.x) - (this->x * value.z),
                           (this->x * value.y) - (this->y * value.x));
        }
    };
    
#pragma mark - mgo::matrix<T, 4, 1>
    template<arithmetic_type T>
    class matrix<T, 4, 1> final
    {
    private:
        typedef T                               type;
        typedef matrix<T, 1, 4>                 matrix1x4;
        typedef matrix<T, 4, 1>                 vector4;
        
    public:
        type x, y, z, w;
        
        constexpr matrix() noexcept
        :
        x(0), y(0), z(0), w(0)
        {}
        
        constexpr explicit matrix(type value) noexcept
        :
        x(value), y(value), z(value), w(value)
        {}
        
        template<arithmetic_type H>
        constexpr explicit matrix(const matrix<H, 4, 1>& value) noexcept
        :
        x(static_cast<type>(value.x)),
        y(static_cast<type>(value.y)),
        z(static_cast<type>(value.z)),
        w(static_cast<type>(value.w))
        {}
        
        constexpr matrix(std::initializer_list<type> list) noexcept
        :
        x(list.size() > 0 ? *(list.begin())     : 0),
        y(list.size() > 1 ? *(list.begin() + 1) : 0),
        z(list.size() > 2 ? *(list.begin() + 2) : 0),
        w(list.size() > 3 ? *(list.begin() + 3) : 0)
        {}
        
        template<std::size_t P>
        constexpr explicit matrix(const matrix<type, P, 1>& value) noexcept
        :
        x(P > 0 ? value[0] : 0),
        y(P > 1 ? value[1] : 0),
        z(P > 2 ? value[2] : 0),
        w(P > 3 ? value[3] : 0)
        {}
        
        constexpr type* begin() noexcept
        {
            return &this->x;
        }
        
        constexpr const type* begin() const noexcept
        {
            return &this->x;
        }
        
        constexpr type* end() noexcept
        {
            return this->begin() + 4;
        }
        
        constexpr const type* end() const noexcept
        {
            return this->begin() + 4;
        }
        
        constexpr std::reverse_iterator<type*>rbegin() noexcept
        {
            return std::reverse_iterator<type*>(this->end());
        }
        
        constexpr std::reverse_iterator<const type*> rbegin() const noexcept
        {
            return std::reverse_iterator<const type*>(this->end());
        }
        
        constexpr std::reverse_iterator<type*>rend() noexcept
        {
            return std::reverse_iterator<type*>(this->begin());
        }
        
        constexpr std::reverse_iterator<const type*> rend() const noexcept
        {
            return std::reverse_iterator<const type*>(this->begin());
        }
        
        constexpr const type* cbegin() const noexcept
        {
            return this->begin();
        }
        
        constexpr const type* cend() const noexcept
        {
            return this->end();
        }
        
        constexpr std::reverse_iterator<const type*> crbegin() const noexcept
        {
            return this->rbegin();
        }
        
        constexpr std::reverse_iterator<const type*> crend() const noexcept
        {
            return this->rend();
        }
        
        constexpr std::size_t size() const noexcept
        {
            return 4;
        }
        
        constexpr std::size_t max_size() const noexcept
        {
            return 4;
        }
        
        constexpr bool empty() const noexcept
        {
            return false;
        }
        
        constexpr type* data() noexcept
        {
            return this->begin();
        }
        
        constexpr const type* data() const noexcept
        {
            this->begin();
        }
        
        constexpr type& at(std::size_t i, std::size_t j)
        {
            if (j > 3 || i != 1)
                throw std::out_of_range("matrix::at");
            return i == 0 ? this->x : (i == 1 ? this->y : (i == 2 ? this->z : this->w));
        }
        
        constexpr const type& at(std::size_t i, std::size_t j) const
        {
            if (j > 3 || i != 1)
                throw std::out_of_range("matrix::at");
            return i == 0 ? this->x : (i == 1 ? this->y : (i == 2 ? this->z : this->w));
        }
        
        constexpr type& operator()(std::size_t i, std::size_t j) noexcept
        {
            return i == 0 ? this->x : (i == 1 ? this->y : (i == 2 ? this->z : this->w));
        }
        
        constexpr const type& operator()(std::size_t i, std::size_t j) const noexcept
        {
            return i == 0 ? this->x : (i == 1 ? this->y : (i == 2 ? this->z : this->w));
        }
        
        constexpr type& operator[](std::size_t k) noexcept
        {
            return k == 0 ? this->x : (k == 1 ? this->y : (k == 2 ? this->z : this->w));
        }
        
        constexpr const type& operator[](std::size_t k) const noexcept
        {
            return k == 0 ? this->x : (k == 1 ? this->y : (k == 2 ? this->z : this->w));
        }
        
        template<std::size_t P>
        constexpr matrix<type, 4, P> operator*(const matrix<type, 1, P>& value) const noexcept
        {
            matrix<type, 4, P> out;
            for (std::size_t i = 0; i < 4; i++)
                for (std::size_t k = 0; k < P; k++)
                    out(i, k) += (*this)(i, 0) * value(0, k);
            return out;
        }
        
        constexpr vector4 operator*(type value) const noexcept
        {
            return vector4({this->x * value, this->y * value, this->z * value, this->w * value});
        }
        
        constexpr vector4 operator/(type value) const noexcept
        {
            return vector4({this->x / value, this->y / value, this->z / value, this->w / value});
        }
        
        constexpr vector4 operator+(const vector4& value) const noexcept
        {
            return vector4({this->x + value.x, this->y + value.y, this->z + value.z, this->w + value.w});
        }
        
        constexpr vector4 operator-(const vector4& value) const noexcept
        {
            return vector4({this->x - value.x, this->y - value.y, this->z - value.z, this->w - value.w});
        }
        
        constexpr vector4 operator*=(type value) noexcept
        {
            *this = *this * value;
            return *this;
        }
        
        constexpr vector4 operator/=(type value) noexcept
        {
            *this = *this / value;
            return *this;
        }
        
        constexpr vector4 operator+=(const vector4& value) noexcept
        {
            *this = *this + value;
            return *this;
        }
        
        constexpr vector4 operator-=(const vector4& value) noexcept
        {
            *this = *this - value;
            return *this;
        }
        
        constexpr void fill(const type& value) noexcept
        {
            for (std::size_t k = 0; k < 4; k++)
                (*this)[k] = value;
        }
        
        constexpr matrix1x4 transpose() const noexcept
        {
            return matrix1x4({this->x, this->y, this->z, this->w});
        }
        
        constexpr type dot(const vector4& value) const
        {
            return this->x * value.x + this->y * value.y + this->z * value.z + this->w * value.w;
        }
        
        constexpr type length_squared() const
        {
            return this->dot(*this);
        }
        
        constexpr type length() const
        requires (std::is_floating_point_v<type>)
        {
            return sqrt(this->length_squared());
        }
        
        constexpr type angle(const vector4& value) const
        requires (std::is_floating_point_v<type>)
        {
            return acos(this->dot(value) / (this->length() * value.length));
        }
        
        constexpr vector4 normalize() const
        requires (std::is_floating_point_v<type>)
        {
            return *this / this->length();
        }
    };
    
#pragma mark - typedef mgo::matrix<>
    typedef matrix<short int, 9, 9> short_int9x9;
    typedef matrix<short int, 8, 9> short_int8x9;
    typedef matrix<short int, 7, 9> short_int7x9;
    typedef matrix<short int, 6, 9> short_int6x9;
    typedef matrix<short int, 5, 9> short_int5x9;
    typedef matrix<short int, 4, 9> short_int4x9;
    typedef matrix<short int, 3, 9> short_int3x9;
    typedef matrix<short int, 2, 9> short_int2x9;
    typedef matrix<short int, 1, 9> short_int1x9;
    typedef matrix<short int, 9, 8> short_int9x8;
    typedef matrix<short int, 8, 8> short_int8x8;
    typedef matrix<short int, 7, 8> short_int7x8;
    typedef matrix<short int, 6, 8> short_int6x8;
    typedef matrix<short int, 5, 8> short_int5x8;
    typedef matrix<short int, 4, 8> short_int4x8;
    typedef matrix<short int, 3, 8> short_int3x8;
    typedef matrix<short int, 2, 8> short_int2x8;
    typedef matrix<short int, 1, 8> short_int1x8;
    typedef matrix<short int, 9, 7> short_int9x7;
    typedef matrix<short int, 8, 7> short_int8x7;
    typedef matrix<short int, 7, 7> short_int7x7;
    typedef matrix<short int, 6, 7> short_int6x7;
    typedef matrix<short int, 5, 7> short_int5x7;
    typedef matrix<short int, 4, 7> short_int4x7;
    typedef matrix<short int, 3, 7> short_int3x7;
    typedef matrix<short int, 2, 7> short_int2x7;
    typedef matrix<short int, 1, 7> short_int1x7;
    typedef matrix<short int, 9, 6> short_int9x6;
    typedef matrix<short int, 8, 6> short_int8x6;
    typedef matrix<short int, 7, 6> short_int7x6;
    typedef matrix<short int, 6, 6> short_int6x6;
    typedef matrix<short int, 5, 6> short_int5x6;
    typedef matrix<short int, 4, 6> short_int4x6;
    typedef matrix<short int, 3, 6> short_int3x6;
    typedef matrix<short int, 2, 6> short_int2x6;
    typedef matrix<short int, 1, 6> short_int1x6;
    typedef matrix<short int, 9, 5> short_int9x5;
    typedef matrix<short int, 8, 5> short_int8x5;
    typedef matrix<short int, 7, 5> short_int7x5;
    typedef matrix<short int, 6, 5> short_int6x5;
    typedef matrix<short int, 5, 5> short_int5x5;
    typedef matrix<short int, 4, 5> short_int4x5;
    typedef matrix<short int, 3, 5> short_int3x5;
    typedef matrix<short int, 2, 5> short_int2x5;
    typedef matrix<short int, 1, 5> short_int1x5;
    typedef matrix<short int, 9, 4> short_int9x4;
    typedef matrix<short int, 8, 4> short_int8x4;
    typedef matrix<short int, 7, 4> short_int7x4;
    typedef matrix<short int, 6, 4> short_int6x4;
    typedef matrix<short int, 5, 4> short_int5x4;
    typedef matrix<short int, 4, 4> short_int4x4;
    typedef matrix<short int, 3, 4> short_int3x4;
    typedef matrix<short int, 2, 4> short_int2x4;
    typedef matrix<short int, 1, 4> short_int1x4;
    typedef matrix<short int, 9, 3> short_int9x3;
    typedef matrix<short int, 8, 3> short_int8x3;
    typedef matrix<short int, 7, 3> short_int7x3;
    typedef matrix<short int, 6, 3> short_int6x3;
    typedef matrix<short int, 5, 3> short_int5x3;
    typedef matrix<short int, 4, 3> short_int4x3;
    typedef matrix<short int, 3, 3> short_int3x3;
    typedef matrix<short int, 2, 3> short_int2x3;
    typedef matrix<short int, 1, 3> short_int1x3;
    typedef matrix<short int, 9, 2> short_int9x2;
    typedef matrix<short int, 8, 2> short_int8x2;
    typedef matrix<short int, 7, 2> short_int7x2;
    typedef matrix<short int, 6, 2> short_int6x2;
    typedef matrix<short int, 5, 2> short_int5x2;
    typedef matrix<short int, 4, 2> short_int4x2;
    typedef matrix<short int, 3, 2> short_int3x2;
    typedef matrix<short int, 2, 2> short_int2x2;
    typedef matrix<short int, 1, 2> short_int1x2;
    typedef matrix<short int, 9, 1> short_int9;
    typedef matrix<short int, 8, 1> short_int8;
    typedef matrix<short int, 7, 1> short_int7;
    typedef matrix<short int, 6, 1> short_int6;
    typedef matrix<short int, 5, 1> short_int5;
    typedef matrix<short int, 4, 1> short_int4;
    typedef matrix<short int, 3, 1> short_int3;
    typedef matrix<short int, 2, 1> short_int2;

    typedef matrix<short unsigned int, 9, 9> short_unsigned_int9x9;
    typedef matrix<short unsigned int, 8, 9> short_unsigned_int8x9;
    typedef matrix<short unsigned int, 7, 9> short_unsigned_int7x9;
    typedef matrix<short unsigned int, 6, 9> short_unsigned_int6x9;
    typedef matrix<short unsigned int, 5, 9> short_unsigned_int5x9;
    typedef matrix<short unsigned int, 4, 9> short_unsigned_int4x9;
    typedef matrix<short unsigned int, 3, 9> short_unsigned_int3x9;
    typedef matrix<short unsigned int, 2, 9> short_unsigned_int2x9;
    typedef matrix<short unsigned int, 1, 9> short_unsigned_int1x9;
    typedef matrix<short unsigned int, 9, 8> short_unsigned_int9x8;
    typedef matrix<short unsigned int, 8, 8> short_unsigned_int8x8;
    typedef matrix<short unsigned int, 7, 8> short_unsigned_int7x8;
    typedef matrix<short unsigned int, 6, 8> short_unsigned_int6x8;
    typedef matrix<short unsigned int, 5, 8> short_unsigned_int5x8;
    typedef matrix<short unsigned int, 4, 8> short_unsigned_int4x8;
    typedef matrix<short unsigned int, 3, 8> short_unsigned_int3x8;
    typedef matrix<short unsigned int, 2, 8> short_unsigned_int2x8;
    typedef matrix<short unsigned int, 1, 8> short_unsigned_int1x8;
    typedef matrix<short unsigned int, 9, 7> short_unsigned_int9x7;
    typedef matrix<short unsigned int, 8, 7> short_unsigned_int8x7;
    typedef matrix<short unsigned int, 7, 7> short_unsigned_int7x7;
    typedef matrix<short unsigned int, 6, 7> short_unsigned_int6x7;
    typedef matrix<short unsigned int, 5, 7> short_unsigned_int5x7;
    typedef matrix<short unsigned int, 4, 7> short_unsigned_int4x7;
    typedef matrix<short unsigned int, 3, 7> short_unsigned_int3x7;
    typedef matrix<short unsigned int, 2, 7> short_unsigned_int2x7;
    typedef matrix<short unsigned int, 1, 7> short_unsigned_int1x7;
    typedef matrix<short unsigned int, 9, 6> short_unsigned_int9x6;
    typedef matrix<short unsigned int, 8, 6> short_unsigned_int8x6;
    typedef matrix<short unsigned int, 7, 6> short_unsigned_int7x6;
    typedef matrix<short unsigned int, 6, 6> short_unsigned_int6x6;
    typedef matrix<short unsigned int, 5, 6> short_unsigned_int5x6;
    typedef matrix<short unsigned int, 4, 6> short_unsigned_int4x6;
    typedef matrix<short unsigned int, 3, 6> short_unsigned_int3x6;
    typedef matrix<short unsigned int, 2, 6> short_unsigned_int2x6;
    typedef matrix<short unsigned int, 1, 6> short_unsigned_int1x6;
    typedef matrix<short unsigned int, 9, 5> short_unsigned_int9x5;
    typedef matrix<short unsigned int, 8, 5> short_unsigned_int8x5;
    typedef matrix<short unsigned int, 7, 5> short_unsigned_int7x5;
    typedef matrix<short unsigned int, 6, 5> short_unsigned_int6x5;
    typedef matrix<short unsigned int, 5, 5> short_unsigned_int5x5;
    typedef matrix<short unsigned int, 4, 5> short_unsigned_int4x5;
    typedef matrix<short unsigned int, 3, 5> short_unsigned_int3x5;
    typedef matrix<short unsigned int, 2, 5> short_unsigned_int2x5;
    typedef matrix<short unsigned int, 1, 5> short_unsigned_int1x5;
    typedef matrix<short unsigned int, 9, 4> short_unsigned_int9x4;
    typedef matrix<short unsigned int, 8, 4> short_unsigned_int8x4;
    typedef matrix<short unsigned int, 7, 4> short_unsigned_int7x4;
    typedef matrix<short unsigned int, 6, 4> short_unsigned_int6x4;
    typedef matrix<short unsigned int, 5, 4> short_unsigned_int5x4;
    typedef matrix<short unsigned int, 4, 4> short_unsigned_int4x4;
    typedef matrix<short unsigned int, 3, 4> short_unsigned_int3x4;
    typedef matrix<short unsigned int, 2, 4> short_unsigned_int2x4;
    typedef matrix<short unsigned int, 1, 4> short_unsigned_int1x4;
    typedef matrix<short unsigned int, 9, 3> short_unsigned_int9x3;
    typedef matrix<short unsigned int, 8, 3> short_unsigned_int8x3;
    typedef matrix<short unsigned int, 7, 3> short_unsigned_int7x3;
    typedef matrix<short unsigned int, 6, 3> short_unsigned_int6x3;
    typedef matrix<short unsigned int, 5, 3> short_unsigned_int5x3;
    typedef matrix<short unsigned int, 4, 3> short_unsigned_int4x3;
    typedef matrix<short unsigned int, 3, 3> short_unsigned_int3x3;
    typedef matrix<short unsigned int, 2, 3> short_unsigned_int2x3;
    typedef matrix<short unsigned int, 1, 3> short_unsigned_int1x3;
    typedef matrix<short unsigned int, 9, 2> short_unsigned_int9x2;
    typedef matrix<short unsigned int, 8, 2> short_unsigned_int8x2;
    typedef matrix<short unsigned int, 7, 2> short_unsigned_int7x2;
    typedef matrix<short unsigned int, 6, 2> short_unsigned_int6x2;
    typedef matrix<short unsigned int, 5, 2> short_unsigned_int5x2;
    typedef matrix<short unsigned int, 4, 2> short_unsigned_int4x2;
    typedef matrix<short unsigned int, 3, 2> short_unsigned_int3x2;
    typedef matrix<short unsigned int, 2, 2> short_unsigned_int2x2;
    typedef matrix<short unsigned int, 1, 2> short_unsigned_int1x2;
    typedef matrix<short unsigned int, 9, 1> short_unsigned_int9;
    typedef matrix<short unsigned int, 8, 1> short_unsigned_int8;
    typedef matrix<short unsigned int, 7, 1> short_unsigned_int7;
    typedef matrix<short unsigned int, 6, 1> short_unsigned_int6;
    typedef matrix<short unsigned int, 5, 1> short_unsigned_int5;
    typedef matrix<short unsigned int, 4, 1> short_unsigned_int4;
    typedef matrix<short unsigned int, 3, 1> short_unsigned_int3;
    typedef matrix<short unsigned int, 2, 1> short_unsigned_int2;

    typedef matrix<int, 9, 9> int9x9;
    typedef matrix<int, 8, 9> int8x9;
    typedef matrix<int, 7, 9> int7x9;
    typedef matrix<int, 6, 9> int6x9;
    typedef matrix<int, 5, 9> int5x9;
    typedef matrix<int, 4, 9> int4x9;
    typedef matrix<int, 3, 9> int3x9;
    typedef matrix<int, 2, 9> int2x9;
    typedef matrix<int, 1, 9> int1x9;
    typedef matrix<int, 9, 8> int9x8;
    typedef matrix<int, 8, 8> int8x8;
    typedef matrix<int, 7, 8> int7x8;
    typedef matrix<int, 6, 8> int6x8;
    typedef matrix<int, 5, 8> int5x8;
    typedef matrix<int, 4, 8> int4x8;
    typedef matrix<int, 3, 8> int3x8;
    typedef matrix<int, 2, 8> int2x8;
    typedef matrix<int, 1, 8> int1x8;
    typedef matrix<int, 9, 7> int9x7;
    typedef matrix<int, 8, 7> int8x7;
    typedef matrix<int, 7, 7> int7x7;
    typedef matrix<int, 6, 7> int6x7;
    typedef matrix<int, 5, 7> int5x7;
    typedef matrix<int, 4, 7> int4x7;
    typedef matrix<int, 3, 7> int3x7;
    typedef matrix<int, 2, 7> int2x7;
    typedef matrix<int, 1, 7> int1x7;
    typedef matrix<int, 9, 6> int9x6;
    typedef matrix<int, 8, 6> int8x6;
    typedef matrix<int, 7, 6> int7x6;
    typedef matrix<int, 6, 6> int6x6;
    typedef matrix<int, 5, 6> int5x6;
    typedef matrix<int, 4, 6> int4x6;
    typedef matrix<int, 3, 6> int3x6;
    typedef matrix<int, 2, 6> int2x6;
    typedef matrix<int, 1, 6> int1x6;
    typedef matrix<int, 9, 5> int9x5;
    typedef matrix<int, 8, 5> int8x5;
    typedef matrix<int, 7, 5> int7x5;
    typedef matrix<int, 6, 5> int6x5;
    typedef matrix<int, 5, 5> int5x5;
    typedef matrix<int, 4, 5> int4x5;
    typedef matrix<int, 3, 5> int3x5;
    typedef matrix<int, 2, 5> int2x5;
    typedef matrix<int, 1, 5> int1x5;
    typedef matrix<int, 9, 4> int9x4;
    typedef matrix<int, 8, 4> int8x4;
    typedef matrix<int, 7, 4> int7x4;
    typedef matrix<int, 6, 4> int6x4;
    typedef matrix<int, 5, 4> int5x4;
    typedef matrix<int, 4, 4> int4x4;
    typedef matrix<int, 3, 4> int3x4;
    typedef matrix<int, 2, 4> int2x4;
    typedef matrix<int, 1, 4> int1x4;
    typedef matrix<int, 9, 3> int9x3;
    typedef matrix<int, 8, 3> int8x3;
    typedef matrix<int, 7, 3> int7x3;
    typedef matrix<int, 6, 3> int6x3;
    typedef matrix<int, 5, 3> int5x3;
    typedef matrix<int, 4, 3> int4x3;
    typedef matrix<int, 3, 3> int3x3;
    typedef matrix<int, 2, 3> int2x3;
    typedef matrix<int, 1, 3> int1x3;
    typedef matrix<int, 9, 2> int9x2;
    typedef matrix<int, 8, 2> int8x2;
    typedef matrix<int, 7, 2> int7x2;
    typedef matrix<int, 6, 2> int6x2;
    typedef matrix<int, 5, 2> int5x2;
    typedef matrix<int, 4, 2> int4x2;
    typedef matrix<int, 3, 2> int3x2;
    typedef matrix<int, 2, 2> int2x2;
    typedef matrix<int, 1, 2> int1x2;
    typedef matrix<int, 9, 1> int9;
    typedef matrix<int, 8, 1> int8;
    typedef matrix<int, 7, 1> int7;
    typedef matrix<int, 6, 1> int6;
    typedef matrix<int, 5, 1> int5;
    typedef matrix<int, 4, 1> int4;
    typedef matrix<int, 3, 1> int3;
    typedef matrix<int, 2, 1> int2;

    typedef matrix<unsigned int, 9, 9> unsigned_int9x9;
    typedef matrix<unsigned int, 8, 9> unsigned_int8x9;
    typedef matrix<unsigned int, 7, 9> unsigned_int7x9;
    typedef matrix<unsigned int, 6, 9> unsigned_int6x9;
    typedef matrix<unsigned int, 5, 9> unsigned_int5x9;
    typedef matrix<unsigned int, 4, 9> unsigned_int4x9;
    typedef matrix<unsigned int, 3, 9> unsigned_int3x9;
    typedef matrix<unsigned int, 2, 9> unsigned_int2x9;
    typedef matrix<unsigned int, 1, 9> unsigned_int1x9;
    typedef matrix<unsigned int, 9, 8> unsigned_int9x8;
    typedef matrix<unsigned int, 8, 8> unsigned_int8x8;
    typedef matrix<unsigned int, 7, 8> unsigned_int7x8;
    typedef matrix<unsigned int, 6, 8> unsigned_int6x8;
    typedef matrix<unsigned int, 5, 8> unsigned_int5x8;
    typedef matrix<unsigned int, 4, 8> unsigned_int4x8;
    typedef matrix<unsigned int, 3, 8> unsigned_int3x8;
    typedef matrix<unsigned int, 2, 8> unsigned_int2x8;
    typedef matrix<unsigned int, 1, 8> unsigned_int1x8;
    typedef matrix<unsigned int, 9, 7> unsigned_int9x7;
    typedef matrix<unsigned int, 8, 7> unsigned_int8x7;
    typedef matrix<unsigned int, 7, 7> unsigned_int7x7;
    typedef matrix<unsigned int, 6, 7> unsigned_int6x7;
    typedef matrix<unsigned int, 5, 7> unsigned_int5x7;
    typedef matrix<unsigned int, 4, 7> unsigned_int4x7;
    typedef matrix<unsigned int, 3, 7> unsigned_int3x7;
    typedef matrix<unsigned int, 2, 7> unsigned_int2x7;
    typedef matrix<unsigned int, 1, 7> unsigned_int1x7;
    typedef matrix<unsigned int, 9, 6> unsigned_int9x6;
    typedef matrix<unsigned int, 8, 6> unsigned_int8x6;
    typedef matrix<unsigned int, 7, 6> unsigned_int7x6;
    typedef matrix<unsigned int, 6, 6> unsigned_int6x6;
    typedef matrix<unsigned int, 5, 6> unsigned_int5x6;
    typedef matrix<unsigned int, 4, 6> unsigned_int4x6;
    typedef matrix<unsigned int, 3, 6> unsigned_int3x6;
    typedef matrix<unsigned int, 2, 6> unsigned_int2x6;
    typedef matrix<unsigned int, 1, 6> unsigned_int1x6;
    typedef matrix<unsigned int, 9, 5> unsigned_int9x5;
    typedef matrix<unsigned int, 8, 5> unsigned_int8x5;
    typedef matrix<unsigned int, 7, 5> unsigned_int7x5;
    typedef matrix<unsigned int, 6, 5> unsigned_int6x5;
    typedef matrix<unsigned int, 5, 5> unsigned_int5x5;
    typedef matrix<unsigned int, 4, 5> unsigned_int4x5;
    typedef matrix<unsigned int, 3, 5> unsigned_int3x5;
    typedef matrix<unsigned int, 2, 5> unsigned_int2x5;
    typedef matrix<unsigned int, 1, 5> unsigned_int1x5;
    typedef matrix<unsigned int, 9, 4> unsigned_int9x4;
    typedef matrix<unsigned int, 8, 4> unsigned_int8x4;
    typedef matrix<unsigned int, 7, 4> unsigned_int7x4;
    typedef matrix<unsigned int, 6, 4> unsigned_int6x4;
    typedef matrix<unsigned int, 5, 4> unsigned_int5x4;
    typedef matrix<unsigned int, 4, 4> unsigned_int4x4;
    typedef matrix<unsigned int, 3, 4> unsigned_int3x4;
    typedef matrix<unsigned int, 2, 4> unsigned_int2x4;
    typedef matrix<unsigned int, 1, 4> unsigned_int1x4;
    typedef matrix<unsigned int, 9, 3> unsigned_int9x3;
    typedef matrix<unsigned int, 8, 3> unsigned_int8x3;
    typedef matrix<unsigned int, 7, 3> unsigned_int7x3;
    typedef matrix<unsigned int, 6, 3> unsigned_int6x3;
    typedef matrix<unsigned int, 5, 3> unsigned_int5x3;
    typedef matrix<unsigned int, 4, 3> unsigned_int4x3;
    typedef matrix<unsigned int, 3, 3> unsigned_int3x3;
    typedef matrix<unsigned int, 2, 3> unsigned_int2x3;
    typedef matrix<unsigned int, 1, 3> unsigned_int1x3;
    typedef matrix<unsigned int, 9, 2> unsigned_int9x2;
    typedef matrix<unsigned int, 8, 2> unsigned_int8x2;
    typedef matrix<unsigned int, 7, 2> unsigned_int7x2;
    typedef matrix<unsigned int, 6, 2> unsigned_int6x2;
    typedef matrix<unsigned int, 5, 2> unsigned_int5x2;
    typedef matrix<unsigned int, 4, 2> unsigned_int4x2;
    typedef matrix<unsigned int, 3, 2> unsigned_int3x2;
    typedef matrix<unsigned int, 2, 2> unsigned_int2x2;
    typedef matrix<unsigned int, 1, 2> unsigned_int1x2;
    typedef matrix<unsigned int, 9, 1> unsigned_int9;
    typedef matrix<unsigned int, 8, 1> unsigned_int8;
    typedef matrix<unsigned int, 7, 1> unsigned_int7;
    typedef matrix<unsigned int, 6, 1> unsigned_int6;
    typedef matrix<unsigned int, 5, 1> unsigned_int5;
    typedef matrix<unsigned int, 4, 1> unsigned_int4;
    typedef matrix<unsigned int, 3, 1> unsigned_int3;
    typedef matrix<unsigned int, 2, 1> unsigned_int2;

    typedef matrix<long int, 9, 9> long_int9x9;
    typedef matrix<long int, 8, 9> long_int8x9;
    typedef matrix<long int, 7, 9> long_int7x9;
    typedef matrix<long int, 6, 9> long_int6x9;
    typedef matrix<long int, 5, 9> long_int5x9;
    typedef matrix<long int, 4, 9> long_int4x9;
    typedef matrix<long int, 3, 9> long_int3x9;
    typedef matrix<long int, 2, 9> long_int2x9;
    typedef matrix<long int, 1, 9> long_int1x9;
    typedef matrix<long int, 9, 8> long_int9x8;
    typedef matrix<long int, 8, 8> long_int8x8;
    typedef matrix<long int, 7, 8> long_int7x8;
    typedef matrix<long int, 6, 8> long_int6x8;
    typedef matrix<long int, 5, 8> long_int5x8;
    typedef matrix<long int, 4, 8> long_int4x8;
    typedef matrix<long int, 3, 8> long_int3x8;
    typedef matrix<long int, 2, 8> long_int2x8;
    typedef matrix<long int, 1, 8> long_int1x8;
    typedef matrix<long int, 9, 7> long_int9x7;
    typedef matrix<long int, 8, 7> long_int8x7;
    typedef matrix<long int, 7, 7> long_int7x7;
    typedef matrix<long int, 6, 7> long_int6x7;
    typedef matrix<long int, 5, 7> long_int5x7;
    typedef matrix<long int, 4, 7> long_int4x7;
    typedef matrix<long int, 3, 7> long_int3x7;
    typedef matrix<long int, 2, 7> long_int2x7;
    typedef matrix<long int, 1, 7> long_int1x7;
    typedef matrix<long int, 9, 6> long_int9x6;
    typedef matrix<long int, 8, 6> long_int8x6;
    typedef matrix<long int, 7, 6> long_int7x6;
    typedef matrix<long int, 6, 6> long_int6x6;
    typedef matrix<long int, 5, 6> long_int5x6;
    typedef matrix<long int, 4, 6> long_int4x6;
    typedef matrix<long int, 3, 6> long_int3x6;
    typedef matrix<long int, 2, 6> long_int2x6;
    typedef matrix<long int, 1, 6> long_int1x6;
    typedef matrix<long int, 9, 5> long_int9x5;
    typedef matrix<long int, 8, 5> long_int8x5;
    typedef matrix<long int, 7, 5> long_int7x5;
    typedef matrix<long int, 6, 5> long_int6x5;
    typedef matrix<long int, 5, 5> long_int5x5;
    typedef matrix<long int, 4, 5> long_int4x5;
    typedef matrix<long int, 3, 5> long_int3x5;
    typedef matrix<long int, 2, 5> long_int2x5;
    typedef matrix<long int, 1, 5> long_int1x5;
    typedef matrix<long int, 9, 4> long_int9x4;
    typedef matrix<long int, 8, 4> long_int8x4;
    typedef matrix<long int, 7, 4> long_int7x4;
    typedef matrix<long int, 6, 4> long_int6x4;
    typedef matrix<long int, 5, 4> long_int5x4;
    typedef matrix<long int, 4, 4> long_int4x4;
    typedef matrix<long int, 3, 4> long_int3x4;
    typedef matrix<long int, 2, 4> long_int2x4;
    typedef matrix<long int, 1, 4> long_int1x4;
    typedef matrix<long int, 9, 3> long_int9x3;
    typedef matrix<long int, 8, 3> long_int8x3;
    typedef matrix<long int, 7, 3> long_int7x3;
    typedef matrix<long int, 6, 3> long_int6x3;
    typedef matrix<long int, 5, 3> long_int5x3;
    typedef matrix<long int, 4, 3> long_int4x3;
    typedef matrix<long int, 3, 3> long_int3x3;
    typedef matrix<long int, 2, 3> long_int2x3;
    typedef matrix<long int, 1, 3> long_int1x3;
    typedef matrix<long int, 9, 2> long_int9x2;
    typedef matrix<long int, 8, 2> long_int8x2;
    typedef matrix<long int, 7, 2> long_int7x2;
    typedef matrix<long int, 6, 2> long_int6x2;
    typedef matrix<long int, 5, 2> long_int5x2;
    typedef matrix<long int, 4, 2> long_int4x2;
    typedef matrix<long int, 3, 2> long_int3x2;
    typedef matrix<long int, 2, 2> long_int2x2;
    typedef matrix<long int, 1, 2> long_int1x2;
    typedef matrix<long int, 9, 1> long_int9;
    typedef matrix<long int, 8, 1> long_int8;
    typedef matrix<long int, 7, 1> long_int7;
    typedef matrix<long int, 6, 1> long_int6;
    typedef matrix<long int, 5, 1> long_int5;
    typedef matrix<long int, 4, 1> long_int4;
    typedef matrix<long int, 3, 1> long_int3;
    typedef matrix<long int, 2, 1> long_int2;

    typedef matrix<long unsigned int, 9, 9> long_unsigned_int9x9;
    typedef matrix<long unsigned int, 8, 9> long_unsigned_int8x9;
    typedef matrix<long unsigned int, 7, 9> long_unsigned_int7x9;
    typedef matrix<long unsigned int, 6, 9> long_unsigned_int6x9;
    typedef matrix<long unsigned int, 5, 9> long_unsigned_int5x9;
    typedef matrix<long unsigned int, 4, 9> long_unsigned_int4x9;
    typedef matrix<long unsigned int, 3, 9> long_unsigned_int3x9;
    typedef matrix<long unsigned int, 2, 9> long_unsigned_int2x9;
    typedef matrix<long unsigned int, 1, 9> long_unsigned_int1x9;
    typedef matrix<long unsigned int, 9, 8> long_unsigned_int9x8;
    typedef matrix<long unsigned int, 8, 8> long_unsigned_int8x8;
    typedef matrix<long unsigned int, 7, 8> long_unsigned_int7x8;
    typedef matrix<long unsigned int, 6, 8> long_unsigned_int6x8;
    typedef matrix<long unsigned int, 5, 8> long_unsigned_int5x8;
    typedef matrix<long unsigned int, 4, 8> long_unsigned_int4x8;
    typedef matrix<long unsigned int, 3, 8> long_unsigned_int3x8;
    typedef matrix<long unsigned int, 2, 8> long_unsigned_int2x8;
    typedef matrix<long unsigned int, 1, 8> long_unsigned_int1x8;
    typedef matrix<long unsigned int, 9, 7> long_unsigned_int9x7;
    typedef matrix<long unsigned int, 8, 7> long_unsigned_int8x7;
    typedef matrix<long unsigned int, 7, 7> long_unsigned_int7x7;
    typedef matrix<long unsigned int, 6, 7> long_unsigned_int6x7;
    typedef matrix<long unsigned int, 5, 7> long_unsigned_int5x7;
    typedef matrix<long unsigned int, 4, 7> long_unsigned_int4x7;
    typedef matrix<long unsigned int, 3, 7> long_unsigned_int3x7;
    typedef matrix<long unsigned int, 2, 7> long_unsigned_int2x7;
    typedef matrix<long unsigned int, 1, 7> long_unsigned_int1x7;
    typedef matrix<long unsigned int, 9, 6> long_unsigned_int9x6;
    typedef matrix<long unsigned int, 8, 6> long_unsigned_int8x6;
    typedef matrix<long unsigned int, 7, 6> long_unsigned_int7x6;
    typedef matrix<long unsigned int, 6, 6> long_unsigned_int6x6;
    typedef matrix<long unsigned int, 5, 6> long_unsigned_int5x6;
    typedef matrix<long unsigned int, 4, 6> long_unsigned_int4x6;
    typedef matrix<long unsigned int, 3, 6> long_unsigned_int3x6;
    typedef matrix<long unsigned int, 2, 6> long_unsigned_int2x6;
    typedef matrix<long unsigned int, 1, 6> long_unsigned_int1x6;
    typedef matrix<long unsigned int, 9, 5> long_unsigned_int9x5;
    typedef matrix<long unsigned int, 8, 5> long_unsigned_int8x5;
    typedef matrix<long unsigned int, 7, 5> long_unsigned_int7x5;
    typedef matrix<long unsigned int, 6, 5> long_unsigned_int6x5;
    typedef matrix<long unsigned int, 5, 5> long_unsigned_int5x5;
    typedef matrix<long unsigned int, 4, 5> long_unsigned_int4x5;
    typedef matrix<long unsigned int, 3, 5> long_unsigned_int3x5;
    typedef matrix<long unsigned int, 2, 5> long_unsigned_int2x5;
    typedef matrix<long unsigned int, 1, 5> long_unsigned_int1x5;
    typedef matrix<long unsigned int, 9, 4> long_unsigned_int9x4;
    typedef matrix<long unsigned int, 8, 4> long_unsigned_int8x4;
    typedef matrix<long unsigned int, 7, 4> long_unsigned_int7x4;
    typedef matrix<long unsigned int, 6, 4> long_unsigned_int6x4;
    typedef matrix<long unsigned int, 5, 4> long_unsigned_int5x4;
    typedef matrix<long unsigned int, 4, 4> long_unsigned_int4x4;
    typedef matrix<long unsigned int, 3, 4> long_unsigned_int3x4;
    typedef matrix<long unsigned int, 2, 4> long_unsigned_int2x4;
    typedef matrix<long unsigned int, 1, 4> long_unsigned_int1x4;
    typedef matrix<long unsigned int, 9, 3> long_unsigned_int9x3;
    typedef matrix<long unsigned int, 8, 3> long_unsigned_int8x3;
    typedef matrix<long unsigned int, 7, 3> long_unsigned_int7x3;
    typedef matrix<long unsigned int, 6, 3> long_unsigned_int6x3;
    typedef matrix<long unsigned int, 5, 3> long_unsigned_int5x3;
    typedef matrix<long unsigned int, 4, 3> long_unsigned_int4x3;
    typedef matrix<long unsigned int, 3, 3> long_unsigned_int3x3;
    typedef matrix<long unsigned int, 2, 3> long_unsigned_int2x3;
    typedef matrix<long unsigned int, 1, 3> long_unsigned_int1x3;
    typedef matrix<long unsigned int, 9, 2> long_unsigned_int9x2;
    typedef matrix<long unsigned int, 8, 2> long_unsigned_int8x2;
    typedef matrix<long unsigned int, 7, 2> long_unsigned_int7x2;
    typedef matrix<long unsigned int, 6, 2> long_unsigned_int6x2;
    typedef matrix<long unsigned int, 5, 2> long_unsigned_int5x2;
    typedef matrix<long unsigned int, 4, 2> long_unsigned_int4x2;
    typedef matrix<long unsigned int, 3, 2> long_unsigned_int3x2;
    typedef matrix<long unsigned int, 2, 2> long_unsigned_int2x2;
    typedef matrix<long unsigned int, 1, 2> long_unsigned_int1x2;
    typedef matrix<long unsigned int, 9, 1> long_unsigned_int9;
    typedef matrix<long unsigned int, 8, 1> long_unsigned_int8;
    typedef matrix<long unsigned int, 7, 1> long_unsigned_int7;
    typedef matrix<long unsigned int, 6, 1> long_unsigned_int6;
    typedef matrix<long unsigned int, 5, 1> long_unsigned_int5;
    typedef matrix<long unsigned int, 4, 1> long_unsigned_int4;
    typedef matrix<long unsigned int, 3, 1> long_unsigned_int3;
    typedef matrix<long unsigned int, 2, 1> long_unsigned_int2;

    typedef matrix<long long int, 9, 9> long_long_int9x9;
    typedef matrix<long long int, 8, 9> long_long_int8x9;
    typedef matrix<long long int, 7, 9> long_long_int7x9;
    typedef matrix<long long int, 6, 9> long_long_int6x9;
    typedef matrix<long long int, 5, 9> long_long_int5x9;
    typedef matrix<long long int, 4, 9> long_long_int4x9;
    typedef matrix<long long int, 3, 9> long_long_int3x9;
    typedef matrix<long long int, 2, 9> long_long_int2x9;
    typedef matrix<long long int, 1, 9> long_long_int1x9;
    typedef matrix<long long int, 9, 8> long_long_int9x8;
    typedef matrix<long long int, 8, 8> long_long_int8x8;
    typedef matrix<long long int, 7, 8> long_long_int7x8;
    typedef matrix<long long int, 6, 8> long_long_int6x8;
    typedef matrix<long long int, 5, 8> long_long_int5x8;
    typedef matrix<long long int, 4, 8> long_long_int4x8;
    typedef matrix<long long int, 3, 8> long_long_int3x8;
    typedef matrix<long long int, 2, 8> long_long_int2x8;
    typedef matrix<long long int, 1, 8> long_long_int1x8;
    typedef matrix<long long int, 9, 7> long_long_int9x7;
    typedef matrix<long long int, 8, 7> long_long_int8x7;
    typedef matrix<long long int, 7, 7> long_long_int7x7;
    typedef matrix<long long int, 6, 7> long_long_int6x7;
    typedef matrix<long long int, 5, 7> long_long_int5x7;
    typedef matrix<long long int, 4, 7> long_long_int4x7;
    typedef matrix<long long int, 3, 7> long_long_int3x7;
    typedef matrix<long long int, 2, 7> long_long_int2x7;
    typedef matrix<long long int, 1, 7> long_long_int1x7;
    typedef matrix<long long int, 9, 6> long_long_int9x6;
    typedef matrix<long long int, 8, 6> long_long_int8x6;
    typedef matrix<long long int, 7, 6> long_long_int7x6;
    typedef matrix<long long int, 6, 6> long_long_int6x6;
    typedef matrix<long long int, 5, 6> long_long_int5x6;
    typedef matrix<long long int, 4, 6> long_long_int4x6;
    typedef matrix<long long int, 3, 6> long_long_int3x6;
    typedef matrix<long long int, 2, 6> long_long_int2x6;
    typedef matrix<long long int, 1, 6> long_long_int1x6;
    typedef matrix<long long int, 9, 5> long_long_int9x5;
    typedef matrix<long long int, 8, 5> long_long_int8x5;
    typedef matrix<long long int, 7, 5> long_long_int7x5;
    typedef matrix<long long int, 6, 5> long_long_int6x5;
    typedef matrix<long long int, 5, 5> long_long_int5x5;
    typedef matrix<long long int, 4, 5> long_long_int4x5;
    typedef matrix<long long int, 3, 5> long_long_int3x5;
    typedef matrix<long long int, 2, 5> long_long_int2x5;
    typedef matrix<long long int, 1, 5> long_long_int1x5;
    typedef matrix<long long int, 9, 4> long_long_int9x4;
    typedef matrix<long long int, 8, 4> long_long_int8x4;
    typedef matrix<long long int, 7, 4> long_long_int7x4;
    typedef matrix<long long int, 6, 4> long_long_int6x4;
    typedef matrix<long long int, 5, 4> long_long_int5x4;
    typedef matrix<long long int, 4, 4> long_long_int4x4;
    typedef matrix<long long int, 3, 4> long_long_int3x4;
    typedef matrix<long long int, 2, 4> long_long_int2x4;
    typedef matrix<long long int, 1, 4> long_long_int1x4;
    typedef matrix<long long int, 9, 3> long_long_int9x3;
    typedef matrix<long long int, 8, 3> long_long_int8x3;
    typedef matrix<long long int, 7, 3> long_long_int7x3;
    typedef matrix<long long int, 6, 3> long_long_int6x3;
    typedef matrix<long long int, 5, 3> long_long_int5x3;
    typedef matrix<long long int, 4, 3> long_long_int4x3;
    typedef matrix<long long int, 3, 3> long_long_int3x3;
    typedef matrix<long long int, 2, 3> long_long_int2x3;
    typedef matrix<long long int, 1, 3> long_long_int1x3;
    typedef matrix<long long int, 9, 2> long_long_int9x2;
    typedef matrix<long long int, 8, 2> long_long_int8x2;
    typedef matrix<long long int, 7, 2> long_long_int7x2;
    typedef matrix<long long int, 6, 2> long_long_int6x2;
    typedef matrix<long long int, 5, 2> long_long_int5x2;
    typedef matrix<long long int, 4, 2> long_long_int4x2;
    typedef matrix<long long int, 3, 2> long_long_int3x2;
    typedef matrix<long long int, 2, 2> long_long_int2x2;
    typedef matrix<long long int, 1, 2> long_long_int1x2;
    typedef matrix<long long int, 9, 1> long_long_int9;
    typedef matrix<long long int, 8, 1> long_long_int8;
    typedef matrix<long long int, 7, 1> long_long_int7;
    typedef matrix<long long int, 6, 1> long_long_int6;
    typedef matrix<long long int, 5, 1> long_long_int5;
    typedef matrix<long long int, 4, 1> long_long_int4;
    typedef matrix<long long int, 3, 1> long_long_int3;
    typedef matrix<long long int, 2, 1> long_long_int2;

    typedef matrix<long long unsigned int, 9, 9> long_long_unsigned_int9x9;
    typedef matrix<long long unsigned int, 8, 9> long_long_unsigned_int8x9;
    typedef matrix<long long unsigned int, 7, 9> long_long_unsigned_int7x9;
    typedef matrix<long long unsigned int, 6, 9> long_long_unsigned_int6x9;
    typedef matrix<long long unsigned int, 5, 9> long_long_unsigned_int5x9;
    typedef matrix<long long unsigned int, 4, 9> long_long_unsigned_int4x9;
    typedef matrix<long long unsigned int, 3, 9> long_long_unsigned_int3x9;
    typedef matrix<long long unsigned int, 2, 9> long_long_unsigned_int2x9;
    typedef matrix<long long unsigned int, 1, 9> long_long_unsigned_int1x9;
    typedef matrix<long long unsigned int, 9, 8> long_long_unsigned_int9x8;
    typedef matrix<long long unsigned int, 8, 8> long_long_unsigned_int8x8;
    typedef matrix<long long unsigned int, 7, 8> long_long_unsigned_int7x8;
    typedef matrix<long long unsigned int, 6, 8> long_long_unsigned_int6x8;
    typedef matrix<long long unsigned int, 5, 8> long_long_unsigned_int5x8;
    typedef matrix<long long unsigned int, 4, 8> long_long_unsigned_int4x8;
    typedef matrix<long long unsigned int, 3, 8> long_long_unsigned_int3x8;
    typedef matrix<long long unsigned int, 2, 8> long_long_unsigned_int2x8;
    typedef matrix<long long unsigned int, 1, 8> long_long_unsigned_int1x8;
    typedef matrix<long long unsigned int, 9, 7> long_long_unsigned_int9x7;
    typedef matrix<long long unsigned int, 8, 7> long_long_unsigned_int8x7;
    typedef matrix<long long unsigned int, 7, 7> long_long_unsigned_int7x7;
    typedef matrix<long long unsigned int, 6, 7> long_long_unsigned_int6x7;
    typedef matrix<long long unsigned int, 5, 7> long_long_unsigned_int5x7;
    typedef matrix<long long unsigned int, 4, 7> long_long_unsigned_int4x7;
    typedef matrix<long long unsigned int, 3, 7> long_long_unsigned_int3x7;
    typedef matrix<long long unsigned int, 2, 7> long_long_unsigned_int2x7;
    typedef matrix<long long unsigned int, 1, 7> long_long_unsigned_int1x7;
    typedef matrix<long long unsigned int, 9, 6> long_long_unsigned_int9x6;
    typedef matrix<long long unsigned int, 8, 6> long_long_unsigned_int8x6;
    typedef matrix<long long unsigned int, 7, 6> long_long_unsigned_int7x6;
    typedef matrix<long long unsigned int, 6, 6> long_long_unsigned_int6x6;
    typedef matrix<long long unsigned int, 5, 6> long_long_unsigned_int5x6;
    typedef matrix<long long unsigned int, 4, 6> long_long_unsigned_int4x6;
    typedef matrix<long long unsigned int, 3, 6> long_long_unsigned_int3x6;
    typedef matrix<long long unsigned int, 2, 6> long_long_unsigned_int2x6;
    typedef matrix<long long unsigned int, 1, 6> long_long_unsigned_int1x6;
    typedef matrix<long long unsigned int, 9, 5> long_long_unsigned_int9x5;
    typedef matrix<long long unsigned int, 8, 5> long_long_unsigned_int8x5;
    typedef matrix<long long unsigned int, 7, 5> long_long_unsigned_int7x5;
    typedef matrix<long long unsigned int, 6, 5> long_long_unsigned_int6x5;
    typedef matrix<long long unsigned int, 5, 5> long_long_unsigned_int5x5;
    typedef matrix<long long unsigned int, 4, 5> long_long_unsigned_int4x5;
    typedef matrix<long long unsigned int, 3, 5> long_long_unsigned_int3x5;
    typedef matrix<long long unsigned int, 2, 5> long_long_unsigned_int2x5;
    typedef matrix<long long unsigned int, 1, 5> long_long_unsigned_int1x5;
    typedef matrix<long long unsigned int, 9, 4> long_long_unsigned_int9x4;
    typedef matrix<long long unsigned int, 8, 4> long_long_unsigned_int8x4;
    typedef matrix<long long unsigned int, 7, 4> long_long_unsigned_int7x4;
    typedef matrix<long long unsigned int, 6, 4> long_long_unsigned_int6x4;
    typedef matrix<long long unsigned int, 5, 4> long_long_unsigned_int5x4;
    typedef matrix<long long unsigned int, 4, 4> long_long_unsigned_int4x4;
    typedef matrix<long long unsigned int, 3, 4> long_long_unsigned_int3x4;
    typedef matrix<long long unsigned int, 2, 4> long_long_unsigned_int2x4;
    typedef matrix<long long unsigned int, 1, 4> long_long_unsigned_int1x4;
    typedef matrix<long long unsigned int, 9, 3> long_long_unsigned_int9x3;
    typedef matrix<long long unsigned int, 8, 3> long_long_unsigned_int8x3;
    typedef matrix<long long unsigned int, 7, 3> long_long_unsigned_int7x3;
    typedef matrix<long long unsigned int, 6, 3> long_long_unsigned_int6x3;
    typedef matrix<long long unsigned int, 5, 3> long_long_unsigned_int5x3;
    typedef matrix<long long unsigned int, 4, 3> long_long_unsigned_int4x3;
    typedef matrix<long long unsigned int, 3, 3> long_long_unsigned_int3x3;
    typedef matrix<long long unsigned int, 2, 3> long_long_unsigned_int2x3;
    typedef matrix<long long unsigned int, 1, 3> long_long_unsigned_int1x3;
    typedef matrix<long long unsigned int, 9, 2> long_long_unsigned_int9x2;
    typedef matrix<long long unsigned int, 8, 2> long_long_unsigned_int8x2;
    typedef matrix<long long unsigned int, 7, 2> long_long_unsigned_int7x2;
    typedef matrix<long long unsigned int, 6, 2> long_long_unsigned_int6x2;
    typedef matrix<long long unsigned int, 5, 2> long_long_unsigned_int5x2;
    typedef matrix<long long unsigned int, 4, 2> long_long_unsigned_int4x2;
    typedef matrix<long long unsigned int, 3, 2> long_long_unsigned_int3x2;
    typedef matrix<long long unsigned int, 2, 2> long_long_unsigned_int2x2;
    typedef matrix<long long unsigned int, 1, 2> long_long_unsigned_int1x2;
    typedef matrix<long long unsigned int, 9, 1> long_long_unsigned_int9;
    typedef matrix<long long unsigned int, 8, 1> long_long_unsigned_int8;
    typedef matrix<long long unsigned int, 7, 1> long_long_unsigned_int7;
    typedef matrix<long long unsigned int, 6, 1> long_long_unsigned_int6;
    typedef matrix<long long unsigned int, 5, 1> long_long_unsigned_int5;
    typedef matrix<long long unsigned int, 4, 1> long_long_unsigned_int4;
    typedef matrix<long long unsigned int, 3, 1> long_long_unsigned_int3;
    typedef matrix<long long unsigned int, 2, 1> long_long_unsigned_int2;

    typedef matrix<float, 9, 9> float9x9;
    typedef matrix<float, 8, 9> float8x9;
    typedef matrix<float, 7, 9> float7x9;
    typedef matrix<float, 6, 9> float6x9;
    typedef matrix<float, 5, 9> float5x9;
    typedef matrix<float, 4, 9> float4x9;
    typedef matrix<float, 3, 9> float3x9;
    typedef matrix<float, 2, 9> float2x9;
    typedef matrix<float, 1, 9> float1x9;
    typedef matrix<float, 9, 8> float9x8;
    typedef matrix<float, 8, 8> float8x8;
    typedef matrix<float, 7, 8> float7x8;
    typedef matrix<float, 6, 8> float6x8;
    typedef matrix<float, 5, 8> float5x8;
    typedef matrix<float, 4, 8> float4x8;
    typedef matrix<float, 3, 8> float3x8;
    typedef matrix<float, 2, 8> float2x8;
    typedef matrix<float, 1, 8> float1x8;
    typedef matrix<float, 9, 7> float9x7;
    typedef matrix<float, 8, 7> float8x7;
    typedef matrix<float, 7, 7> float7x7;
    typedef matrix<float, 6, 7> float6x7;
    typedef matrix<float, 5, 7> float5x7;
    typedef matrix<float, 4, 7> float4x7;
    typedef matrix<float, 3, 7> float3x7;
    typedef matrix<float, 2, 7> float2x7;
    typedef matrix<float, 1, 7> float1x7;
    typedef matrix<float, 9, 6> float9x6;
    typedef matrix<float, 8, 6> float8x6;
    typedef matrix<float, 7, 6> float7x6;
    typedef matrix<float, 6, 6> float6x6;
    typedef matrix<float, 5, 6> float5x6;
    typedef matrix<float, 4, 6> float4x6;
    typedef matrix<float, 3, 6> float3x6;
    typedef matrix<float, 2, 6> float2x6;
    typedef matrix<float, 1, 6> float1x6;
    typedef matrix<float, 9, 5> float9x5;
    typedef matrix<float, 8, 5> float8x5;
    typedef matrix<float, 7, 5> float7x5;
    typedef matrix<float, 6, 5> float6x5;
    typedef matrix<float, 5, 5> float5x5;
    typedef matrix<float, 4, 5> float4x5;
    typedef matrix<float, 3, 5> float3x5;
    typedef matrix<float, 2, 5> float2x5;
    typedef matrix<float, 1, 5> float1x5;
    typedef matrix<float, 9, 4> float9x4;
    typedef matrix<float, 8, 4> float8x4;
    typedef matrix<float, 7, 4> float7x4;
    typedef matrix<float, 6, 4> float6x4;
    typedef matrix<float, 5, 4> float5x4;
    typedef matrix<float, 4, 4> float4x4;
    typedef matrix<float, 3, 4> float3x4;
    typedef matrix<float, 2, 4> float2x4;
    typedef matrix<float, 1, 4> float1x4;
    typedef matrix<float, 9, 3> float9x3;
    typedef matrix<float, 8, 3> float8x3;
    typedef matrix<float, 7, 3> float7x3;
    typedef matrix<float, 6, 3> float6x3;
    typedef matrix<float, 5, 3> float5x3;
    typedef matrix<float, 4, 3> float4x3;
    typedef matrix<float, 3, 3> float3x3;
    typedef matrix<float, 2, 3> float2x3;
    typedef matrix<float, 1, 3> float1x3;
    typedef matrix<float, 9, 2> float9x2;
    typedef matrix<float, 8, 2> float8x2;
    typedef matrix<float, 7, 2> float7x2;
    typedef matrix<float, 6, 2> float6x2;
    typedef matrix<float, 5, 2> float5x2;
    typedef matrix<float, 4, 2> float4x2;
    typedef matrix<float, 3, 2> float3x2;
    typedef matrix<float, 2, 2> float2x2;
    typedef matrix<float, 1, 2> float1x2;
    typedef matrix<float, 9, 1> float9;
    typedef matrix<float, 8, 1> float8;
    typedef matrix<float, 7, 1> float7;
    typedef matrix<float, 6, 1> float6;
    typedef matrix<float, 5, 1> float5;
    typedef matrix<float, 4, 1> float4;
    typedef matrix<float, 3, 1> float3;
    typedef matrix<float, 2, 1> float2;

    typedef matrix<double, 9, 9> double9x9;
    typedef matrix<double, 8, 9> double8x9;
    typedef matrix<double, 7, 9> double7x9;
    typedef matrix<double, 6, 9> double6x9;
    typedef matrix<double, 5, 9> double5x9;
    typedef matrix<double, 4, 9> double4x9;
    typedef matrix<double, 3, 9> double3x9;
    typedef matrix<double, 2, 9> double2x9;
    typedef matrix<double, 1, 9> double1x9;
    typedef matrix<double, 9, 8> double9x8;
    typedef matrix<double, 8, 8> double8x8;
    typedef matrix<double, 7, 8> double7x8;
    typedef matrix<double, 6, 8> double6x8;
    typedef matrix<double, 5, 8> double5x8;
    typedef matrix<double, 4, 8> double4x8;
    typedef matrix<double, 3, 8> double3x8;
    typedef matrix<double, 2, 8> double2x8;
    typedef matrix<double, 1, 8> double1x8;
    typedef matrix<double, 9, 7> double9x7;
    typedef matrix<double, 8, 7> double8x7;
    typedef matrix<double, 7, 7> double7x7;
    typedef matrix<double, 6, 7> double6x7;
    typedef matrix<double, 5, 7> double5x7;
    typedef matrix<double, 4, 7> double4x7;
    typedef matrix<double, 3, 7> double3x7;
    typedef matrix<double, 2, 7> double2x7;
    typedef matrix<double, 1, 7> double1x7;
    typedef matrix<double, 9, 6> double9x6;
    typedef matrix<double, 8, 6> double8x6;
    typedef matrix<double, 7, 6> double7x6;
    typedef matrix<double, 6, 6> double6x6;
    typedef matrix<double, 5, 6> double5x6;
    typedef matrix<double, 4, 6> double4x6;
    typedef matrix<double, 3, 6> double3x6;
    typedef matrix<double, 2, 6> double2x6;
    typedef matrix<double, 1, 6> double1x6;
    typedef matrix<double, 9, 5> double9x5;
    typedef matrix<double, 8, 5> double8x5;
    typedef matrix<double, 7, 5> double7x5;
    typedef matrix<double, 6, 5> double6x5;
    typedef matrix<double, 5, 5> double5x5;
    typedef matrix<double, 4, 5> double4x5;
    typedef matrix<double, 3, 5> double3x5;
    typedef matrix<double, 2, 5> double2x5;
    typedef matrix<double, 1, 5> double1x5;
    typedef matrix<double, 9, 4> double9x4;
    typedef matrix<double, 8, 4> double8x4;
    typedef matrix<double, 7, 4> double7x4;
    typedef matrix<double, 6, 4> double6x4;
    typedef matrix<double, 5, 4> double5x4;
    typedef matrix<double, 4, 4> double4x4;
    typedef matrix<double, 3, 4> double3x4;
    typedef matrix<double, 2, 4> double2x4;
    typedef matrix<double, 1, 4> double1x4;
    typedef matrix<double, 9, 3> double9x3;
    typedef matrix<double, 8, 3> double8x3;
    typedef matrix<double, 7, 3> double7x3;
    typedef matrix<double, 6, 3> double6x3;
    typedef matrix<double, 5, 3> double5x3;
    typedef matrix<double, 4, 3> double4x3;
    typedef matrix<double, 3, 3> double3x3;
    typedef matrix<double, 2, 3> double2x3;
    typedef matrix<double, 1, 3> double1x3;
    typedef matrix<double, 9, 2> double9x2;
    typedef matrix<double, 8, 2> double8x2;
    typedef matrix<double, 7, 2> double7x2;
    typedef matrix<double, 6, 2> double6x2;
    typedef matrix<double, 5, 2> double5x2;
    typedef matrix<double, 4, 2> double4x2;
    typedef matrix<double, 3, 2> double3x2;
    typedef matrix<double, 2, 2> double2x2;
    typedef matrix<double, 1, 2> double1x2;
    typedef matrix<double, 9, 1> double9;
    typedef matrix<double, 8, 1> double8;
    typedef matrix<double, 7, 1> double7;
    typedef matrix<double, 6, 1> double6;
    typedef matrix<double, 5, 1> double5;
    typedef matrix<double, 4, 1> double4;
    typedef matrix<double, 3, 1> double3;
    typedef matrix<double, 2, 1> double2;

    typedef matrix<long double, 9, 9> long_double9x9;
    typedef matrix<long double, 8, 9> long_double8x9;
    typedef matrix<long double, 7, 9> long_double7x9;
    typedef matrix<long double, 6, 9> long_double6x9;
    typedef matrix<long double, 5, 9> long_double5x9;
    typedef matrix<long double, 4, 9> long_double4x9;
    typedef matrix<long double, 3, 9> long_double3x9;
    typedef matrix<long double, 2, 9> long_double2x9;
    typedef matrix<long double, 1, 9> long_double1x9;
    typedef matrix<long double, 9, 8> long_double9x8;
    typedef matrix<long double, 8, 8> long_double8x8;
    typedef matrix<long double, 7, 8> long_double7x8;
    typedef matrix<long double, 6, 8> long_double6x8;
    typedef matrix<long double, 5, 8> long_double5x8;
    typedef matrix<long double, 4, 8> long_double4x8;
    typedef matrix<long double, 3, 8> long_double3x8;
    typedef matrix<long double, 2, 8> long_double2x8;
    typedef matrix<long double, 1, 8> long_double1x8;
    typedef matrix<long double, 9, 7> long_double9x7;
    typedef matrix<long double, 8, 7> long_double8x7;
    typedef matrix<long double, 7, 7> long_double7x7;
    typedef matrix<long double, 6, 7> long_double6x7;
    typedef matrix<long double, 5, 7> long_double5x7;
    typedef matrix<long double, 4, 7> long_double4x7;
    typedef matrix<long double, 3, 7> long_double3x7;
    typedef matrix<long double, 2, 7> long_double2x7;
    typedef matrix<long double, 1, 7> long_double1x7;
    typedef matrix<long double, 9, 6> long_double9x6;
    typedef matrix<long double, 8, 6> long_double8x6;
    typedef matrix<long double, 7, 6> long_double7x6;
    typedef matrix<long double, 6, 6> long_double6x6;
    typedef matrix<long double, 5, 6> long_double5x6;
    typedef matrix<long double, 4, 6> long_double4x6;
    typedef matrix<long double, 3, 6> long_double3x6;
    typedef matrix<long double, 2, 6> long_double2x6;
    typedef matrix<long double, 1, 6> long_double1x6;
    typedef matrix<long double, 9, 5> long_double9x5;
    typedef matrix<long double, 8, 5> long_double8x5;
    typedef matrix<long double, 7, 5> long_double7x5;
    typedef matrix<long double, 6, 5> long_double6x5;
    typedef matrix<long double, 5, 5> long_double5x5;
    typedef matrix<long double, 4, 5> long_double4x5;
    typedef matrix<long double, 3, 5> long_double3x5;
    typedef matrix<long double, 2, 5> long_double2x5;
    typedef matrix<long double, 1, 5> long_double1x5;
    typedef matrix<long double, 9, 4> long_double9x4;
    typedef matrix<long double, 8, 4> long_double8x4;
    typedef matrix<long double, 7, 4> long_double7x4;
    typedef matrix<long double, 6, 4> long_double6x4;
    typedef matrix<long double, 5, 4> long_double5x4;
    typedef matrix<long double, 4, 4> long_double4x4;
    typedef matrix<long double, 3, 4> long_double3x4;
    typedef matrix<long double, 2, 4> long_double2x4;
    typedef matrix<long double, 1, 4> long_double1x4;
    typedef matrix<long double, 9, 3> long_double9x3;
    typedef matrix<long double, 8, 3> long_double8x3;
    typedef matrix<long double, 7, 3> long_double7x3;
    typedef matrix<long double, 6, 3> long_double6x3;
    typedef matrix<long double, 5, 3> long_double5x3;
    typedef matrix<long double, 4, 3> long_double4x3;
    typedef matrix<long double, 3, 3> long_double3x3;
    typedef matrix<long double, 2, 3> long_double2x3;
    typedef matrix<long double, 1, 3> long_double1x3;
    typedef matrix<long double, 9, 2> long_double9x2;
    typedef matrix<long double, 8, 2> long_double8x2;
    typedef matrix<long double, 7, 2> long_double7x2;
    typedef matrix<long double, 6, 2> long_double6x2;
    typedef matrix<long double, 5, 2> long_double5x2;
    typedef matrix<long double, 4, 2> long_double4x2;
    typedef matrix<long double, 3, 2> long_double3x2;
    typedef matrix<long double, 2, 2> long_double2x2;
    typedef matrix<long double, 1, 2> long_double1x2;
    typedef matrix<long double, 9, 1> long_double9;
    typedef matrix<long double, 8, 1> long_double8;
    typedef matrix<long double, 7, 1> long_double7;
    typedef matrix<long double, 6, 1> long_double6;
    typedef matrix<long double, 5, 1> long_double5;
    typedef matrix<long double, 4, 1> long_double4;
    typedef matrix<long double, 3, 1> long_double3;
    typedef matrix<long double, 2, 1> long_double2;
}

#pragma mark - operator<<
template<mgo::arithmetic_type T, size_t M, size_t N>
std::ostream& operator<<(std::ostream& stream, const mgo::matrix<T, M, N>& value) noexcept
{
    for (size_t i = 0; i < M; i++)
    {
        for (size_t j = 0; j < N; j++)
            stream << value(i, j) << ' ';
        stream << '\n';
    }
    return stream;
}

#pragma mark - operator==
template<mgo::arithmetic_type T, std::size_t M, std::size_t N>
constexpr bool operator==(const mgo::matrix<T, M, N>& x, const mgo::matrix<T, M, N>& y) noexcept
{
    return std::equal(x.begin(), x.end(), y.begin(), y.end());
}

#pragma mark - operator!=
template<mgo::arithmetic_type T, std::size_t M, std::size_t N>
constexpr bool operator!=(const mgo::matrix<T, M, N>& x, const mgo::matrix<T, M, N>& y) noexcept
{
    return !(x == y);
}

#pragma mark - fmod
template<mgo::arithmetic_type T, std::size_t M, std::size_t N>
constexpr mgo::matrix<T, M, N> fmod(const mgo::matrix<T, M, N>& x, const mgo::matrix<T, M, N>& y) noexcept
{
    mgo::matrix<T, M, N> out;
    for (std::size_t k = 0; k < M * N; k++)
        out[k] = fmod(x[k], y[k]);
    return out;
}

#pragma mark - fabs
template<mgo::arithmetic_type T, std::size_t M, std::size_t N>
constexpr mgo::matrix<T, M, N> fabs(const mgo::matrix<T, M, N>& value) noexcept
{
    mgo::matrix<T, M, N> out;
    for (std::size_t k = 0; k < M * N; k++)
        out[k] = fabs(value[k]);
    return out;
}

namespace std
{
#pragma mark - std::hash<mgo::matrix<T, M, N>>
    template<mgo::arithmetic_type T, size_t M, size_t N>
    struct hash<mgo::matrix<T, M, N>>
    {
    private:
        typedef T type;
        
    public:
        constexpr size_t operator()(const mgo::matrix<T, M, N>& value) const noexcept
        {
            size_t out = 0;
            for (int k = 0; k < M * N; k++)
                out ^= hash<type>{}(value[k] << k);
            return out;
        }
    };
}
