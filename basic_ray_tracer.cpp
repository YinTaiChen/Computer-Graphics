// how to compile:
// c++ -o basic_ray_tracer basic_ray_tracer.cpp -std=c++11
// Author: YinTaiChen 

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <memory>
#include <utility>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <iomanip>

using namespace std;

#if defined __linux__ || defined __APPLE__
// "Compiled for Linux
#else
// Windows doesn't define these values by default, Linux does
#define M_PI 3.141592653589793
#endif


template<typename T>
class Vec3
{
public:
    T x, y, z;
    Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(T xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
    Vec3& normalize()
    {
        T nor2 = length2();
        if (nor2 > 0) {
            T invNor = 1 / sqrt(nor2);
            x *= invNor, y *= invNor, z *= invNor;
        }
        return *this;
    }
    Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
    Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
    T dot(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
    T dotProduct(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3<T> crossProduct(const Vec3<T> &v) const { return Vec3<T>(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);}
    Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
    Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
    Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
    Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
    Vec3<T>& operator /= (const Vec3<T> &v) { x /= v.x, y /= v.y, z /= v.z; return *this; }
    Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
    T length2() const { return x * x + y * y + z * z; }
    T length() const { return sqrt(length2()); }
    const T& operator [] (uint8_t i) const { return (&x)[i]; }
    T& operator [] (uint8_t i) { return (&x)[i]; }
    friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v)
    {
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
};

typedef Vec3<float> Vec3f;

template<typename T>
class Matrix44
{
	public:
		
		T x[4][4] = {{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}};
		
		Matrix44() {};
		
		Matrix44 (T a, T b, T c, T d, T e, T f, T g, T h,
					T i, T j, T k, T l, T m, T n, T o, T p)
		{
			x[0][0] = a;
			x[0][1] = b;
			x[0][2] = c;
			x[0][3] = d;
			x[1][0] = e;
			x[1][1] = f;
			x[1][2] = g;
			x[1][3] = h;
			x[2][0] = i;
			x[2][1] = j;
			x[2][2] = k;
			x[2][3] = l;
			x[3][0] = m;
			x[3][1] = n;
			x[3][2] = o;
			x[3][3] = p;
		}
		
		const T* operator [] (uint8_t i) const { return x[i]; }
		T* operator [] (uint8_t i) { return x[i]; }
		
		// Multiply the current matrix with another matrix (rhs)
		Matrix44 operator * (const Matrix44& v) const
		{
			Matrix44 tmp;
			multiply (*this, v, tmp);
			
			return tmp;
		}
		
		static void multiply(const Matrix44<T> &a, const Matrix44& b, Matrix44 &c)
		{
			for (uint8_t i = 0; i < 4; ++i){
				for (uint8_t j = 0; j < 4; ++j){
					c[i][j] = a[i][0] * b[0][j] + a[i][1] * b [1][j] + a[i][2] * b[2][j] + a[i][3] * b[3][j];
				}
			}
		}
		
		// return a transposed copy of the current matrix
		Matrix44 transposed() const
		{
			Matrix44 t;
			for (uint8_t i = 0; i < 4; ++i) {
				for (uint8_t j = 0; j < 4; ++j) {
					t[i][j] = x[j][i];
				}
			}
		}
		
		// transpose a matrix itself
		Matrix44& transpose()
		{
			Matrix44 tmp (x[0][0], x[1][0], x[2][0], x[3][0],
							x[1][0], x[1][1], x[1][2], x[1][3],
							x[2][0], x[2][1], x[2][2], x[2][3],
							x[3][0], x[3][1], x[3][2], x[3][3]);
							
			*this = tmp;
			
			return *this;
		}
		
		// point-matrix multiplication
		template<typename S>
		void multVecMatrix(const Vec3<S> &src, Vec3<S> &dst) const
		{
			S a, b, c, w;
			
			a = src[0] * x[0][0] + src[1] * x[1][0] + src[2] * x[2][0] + x[3][0];
			b = src[0] * x[0][1] + src[1] * x[1][1] + src[2] * x[2][1] + x[3][1];
			c = src[0] * x[0][2] + src[1] * x[1][2] + src[2] * x[2][2] + x[3][2];
			w = src[0] * x[0][3] + src[1] * x[1][3] + src[2] * x[2][3] + x[3][3];
			
			dst.x = a / w;
			dst.y = b / w;
			dst.z = c / w;
		}
		
		// vector-matrix multiplication
		template<typename S>
		void multDirMatrix(const Vec3<S> &src, Vec3<S> &dst) const
		{
			S a, b, c;
			
			a = src[0] * x[0][0] + src[1] * x[1][0] + src[2] * x[2][0];
			b = src[0] * x[0][1] + src[1] * x[1][1] + src[2] * x[2][1];
			c = src[0] * x[0][2] + src[1] * x[1][2] + src[2] * x[2][2];
			
			dst.x = a;
			dst.y = b;
			dst.z = c;
		}
		
		// compute the inverse of the matirx using the Gauss-Jordan elimination method
		Matrix44 inverse() const
		{
			int i, j, k;
			Matrix44 s;
			Matrix44 t (*this);
			
			// forward elimination
			for (i = 0; i < 3; i++){
				int pivot = i;
				
				T pivotsize = t[i][i];
				
				if (pivotsize < 0)
					pivotsize = -pivotsize;
					
					for (j = i + 1; j < 4; j++) {
						T tmp = t[j][i];
						
						if (tmp < 0)
							tmp = -tmp;
							
							if (tmp > pivotsize) {
								pivot = j;
								pivotsize = tmp;
							}
					}
					
				if (pivotsize == 0){
					// cannot invert singular matrix
					return Matrix44();
				}
				
				if (pivot != i){
					for (j = 0; j < 4; j++){
						T tmp;
						
						tmp = t[i][j];
						t[i][j] = t[pivot][j];
						t[pivot][j] = tmp;
						
						tmp = s[i][j];
						s[i][j] = s[pivot][j];
						s[pivot][j] = tmp;
					}
				}
				
				for (j = i + 1; j < 4; j++){
					T f = t[j][i] / t[i][i];
					
					for (k = 0; k < 4; k++){
						t[j][k] -= f * t[i][k];
						s[j][k] -= f * s[i][k];
					}
				}
			}
			
			// backward substitution
			for (i = 3; i >= 0; --i){
				T f;
				
				if ((f == t[i][i]) == 0){
					// cannot invert singular matrix
					return Matrix44();
				}
				
				for (j = 0; j < 4; j++){
					t[i][j] /= f;
					s[i][j] /= f;
				}
				
				for (j = 0; j < i; j++){
					f = t[j][i];
					
					for (k = 0; k < 4; k++){
						t[j][k] -= f * t[i][k];
						s[j][k] -= f * s[i][k];
					}
				}
			}
			
			return s;
		}
		
		// set current matrix to its inverse
		const Matrix44<T>& invert()
		{
			*this = inverse();
			return *this;
		}
		
		friend std::ostream& operator << (std::ostream &s, const Matrix44 &m)
		{
			std::ios_base::fmtflags oldFlags = s.flags();
			int width = 12; // total with of the displayed number
			s.precision(5); // control the number of displayed decimals
			s.setf (std::ios_base::fixed);
			
			s << "[" << std::setw (width) << m[0][0] <<
				 " " << std::setw (width) << m[0][1] <<
				 " " << std::setw (width) << m[0][2] <<
				 " " << std::setw (width) << m[0][3] << "\n" <<
				 
				 " " << std::setw (width) << m[1][0] <<
				 " " << std::setw (width) << m[1][1] <<
				 " " << std::setw (width) << m[1][2] <<
				 " " << std::setw (width) << m[1][3] << "\n" <<
				 
				 " " << std::setw (width) << m[2][0] <<
				 " " << std::setw (width) << m[2][1] <<
				 " " << std::setw (width) << m[2][2] <<
				 " " << std::setw (width) << m[2][3] << "\n" <<
				 
				 " " << std::setw (width) << m[3][0] <<
				 " " << std::setw (width) << m[3][1] <<
				 " " << std::setw (width) << m[3][2] <<
				 " " << std::setw (width) << m[3][3] << "]";
			
			s.flags (oldFlags);
			return s;
		}
};

typedef Matrix44<float> Matrix44f;

class Sphere
{
	public:
	    Vec3f center;                           /// position of the sphere
	    float radius, radius2;                  /// sphere radius and radius^2
	    Vec3f surfaceColor, emissionColor;      /// surface color and emission (light)
	    float transparency, reflection;         /// surface transparency and reflectivity
	    Sphere(
	        const Vec3f &c,
	        const float &r,
	        const Vec3f &sc,
	        const float &refl = 0,
	        const float &transp = 0,
	        const Vec3f &ec = 0) :
	        center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
	        transparency(transp), reflection(refl)
	    { /* empty */ }
	
	    // Compute a ray-sphere intersection using the geometric solution
	    bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
	    {
	        Vec3f l = center - rayorig;
	        float tca = l.dot(raydir);
	        if (tca < 0) return false;
	        float d2 = l.dot(l) - tca * tca;
	        if (d2 > radius2) return false;
	        float thc = sqrt(radius2 - d2);
	        t0 = tca - thc;
	        t1 = tca + thc;
	        
	        return true;
	    }
};

class Triangle
{
	public:
		Vec3f v0, v1, v2;
		Triangle(
			const Vec3f &vertice0,
			const Vec3f &vertice1,
			const Vec3f &vertice2) :
			v0(vertice0), v1(vertice1), v2(vertice2)
			{ /* empty */ }
		
		// Compute a ray-triangle intersection using the geometric solution
		bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t, float &u, float &v) const
		{
			float kEpsilon = 1e-8;
	
			Vec3f v0v1 = v1 - v0;
			Vec3f v0v2 = v2 - v0;
			Vec3f N = v0v1.crossProduct(v0v2);
			float denom = N.dotProduct(N);
			
			// finding P
			
			// check if ray and plane are parallel
			float NdotRayDirection = N.dotProduct(raydir);
			if (fabs(NdotRayDirection) < kEpsilon) // almost 0
				return false; // they are parallel so they don't intersect!
			
			float d = N.dotProduct(v0);
			
			t = (N.dotProduct(rayorig) + d) / NdotRayDirection;
			// check if the triangle is in behind the ray
			if (t < 0) return false; // the triangle is behind
			
			// compute the intersection point
			Vec3f P = rayorig + raydir * t;
			
			Vec3f C; // a vector perpendicular to triangle's plane
			
			// edge 0
			Vec3f edge0 = v1 - v0;
			Vec3f vp0  = P - v0;
			C = edge0.crossProduct(vp0);
			if (N.dotProduct(C) < 0) return false; // P is on the right side
			
			// edge 1
			Vec3f edge1 = v2 - v1;
			Vec3f vp1 = P - v1;
			C = edge1.crossProduct(vp1);
			if ((u = N.dotProduct(C)) < 0) return false; // P is on the right side
			
			// edge 2
			Vec3f edge2 = v0 - v2;
			Vec3f vp2 = P - v2;
			C = edge2.crossProduct(vp2);
			if ((v = N.dotProduct(C)) < 0) return false; // P is on the right side
			
			u /= denom;
			v /= denom;
			
			return true; // this ray hits the triangle
		}
};


// This variable controls the maximum recursion depth
#define MAX_RAY_DEPTH 5

float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}

Vec3f traceSphere(
    const Vec3f &rayorig,
    const Vec3f &raydir,
    const std::vector<Sphere> &spheres,
    const int &depth)
{
    float tnear = INFINITY;
    const Sphere* sphere = NULL;
    // find intersection of this ray with the sphere in the scene
    for (unsigned i = 0; i < spheres.size(); ++i) {
        float t0 = INFINITY, t1 = INFINITY;
        if (spheres[i].intersect(rayorig, raydir, t0, t1)) {
            if (t0 < 0) t0 = t1;
            if (t0 < tnear) {
                tnear = t0;
                sphere = &spheres[i];
            }
        }
    }
    // if there's no intersection return black or background color
    if (!sphere) return Vec3f(0.0);
    else return Vec3f(1.0);
}

void render(
	const std::vector<Sphere> &spheres,
	const std::vector<Triangle> &triangles,
	const Vec3f eye,
	const Vec3f view,
	const float fov,
	const unsigned width,
	const unsigned height)
{
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.);
    
    Vec3f forward = (-view).normalize();
    Vec3f right = Vec3f(0, 1, 0).normalize().crossProduct(forward);
    Vec3f up = forward.crossProduct(right);
    
    Matrix44f camToWorld;
    
    camToWorld[0][0] = right.x; 
    camToWorld[0][1] = right.y; 
    camToWorld[0][2] = right.z; 
    camToWorld[1][0] = up.x; 
    camToWorld[1][1] = up.y; 
    camToWorld[1][2] = up.z; 
    camToWorld[2][0] = forward.x; 
    camToWorld[2][1] = forward.y; 
    camToWorld[2][2] = forward.z; 
 
    camToWorld[3][0] = eye.x; 
    camToWorld[3][1] = eye.y; 
    camToWorld[3][2] = eye.z;
    
    Vec3f orig = eye;
    
    // Ray Tracing
    for (unsigned y = 0; y < height; ++y) {
        for (unsigned x = 0; x < width; ++x, ++pixel) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            // direction of the ray
            Vec3f raydir;
            camToWorld.multDirMatrix(Vec3f(xx, yy, -1), raydir);
            raydir.normalize();
            
            // for spheres
            *pixel = traceSphere(orig, raydir, spheres, 0);
            
            // for triangles
            float t, u, v;
    		for (unsigned k = 0; k < triangles.size(); k++){
    			if (triangles[k].intersect(orig, raydir, t, u, v)){
    				*pixel = Vec3f(1.0, 1.0, 1.0);
				}
			}
        }
    }
    
    // Save result to a PPM image
    std::ofstream ofs("./out_both.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (unsigned i = 0; i < width * height; ++i) {
        ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
               (unsigned char)(std::min(float(1), image[i].y) * 255) <<
               (unsigned char)(std::min(float(1), image[i].z) * 255);
    }
    ofs.close();
    delete [] image;  
}

int main(int argc, char **argv)
{
    std::vector<Sphere> spheres;
    std::vector<Triangle> triangles;
    
    Vec3f eye;
	Vec3f view;
	float fov;
	unsigned width;
	unsigned height;
	
	string line;
	
	ifstream myfile("scene_setting_01.txt");
	if (myfile.is_open())
	{
		while ( getline(myfile, line))
		{
			std::istringstream in(line);
			
			std::string type;
			in >> type;
			
			if (type == "E")
			{
				float x, y, z;
				in >> x >> y >> z;
				eye = Vec3f(x, y, z);
			}
			else if (type == "V")
			{
				float x, y, z;
				in >> x >> y >> z;
				view = Vec3f(x, y, z);
			}
			else if (type == "F")
			{
				float temp;
				in >> temp;
				fov = temp;
			}
			else if (type == "R")
			{
				unsigned w, h;
				in >> w >> h;
				width = w;
				height = h;
			}
			else if (type == "S")
			{
				float x, y, z, r;
				Vec3f position;
				in >> x >> y >> z >> r;
				position = Vec3f(x, y, z);
				
				// position, radius, surface color, reflectivity, transparency, emission color
				spheres.push_back(Sphere(position, r, Vec3f(1.0, 1.0, 1.0), 0.0, 0.0));
			}
			else if (type == "T")
			{
				float v00, v10, v20, v01, v11, v21, v02, v12, v22;
				Vec3f v0, v1, v2;
				in >> v00 >> v10 >> v20 >> v01 >> v11 >> v21 >> v02 >> v12 >> v22;
				v0 = Vec3f(v00, v10, v20);
				v1 = Vec3f(v01, v11, v21);
				v2 = Vec3f(v02, v12, v22);
				
				// vertice0, vertice1, vertice2
				triangles.push_back(Triangle(v0, v1, v2));
			}
		}
		myfile.close();
	}
    
    render(spheres, triangles, eye, view, fov, width, height);
	
    return 0;
}
