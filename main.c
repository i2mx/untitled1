#include <float.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <windows.h>

typedef struct Vector3D {
    double i;
    double j;
    double k;
} Vector3D;

Vector3D Vector3DPrint(const Vector3D v) {
    printf("[%f, %f, %f]", v.i, v.j, v.k);
}

Vector3D Vector3DAddition(const Vector3D a, const Vector3D b) {
    return (Vector3D) { a.i + b.i, a.j + b.j, a.k + b.k };
}

Vector3D Vector3DSubtraction(const Vector3D a, const Vector3D b) {
    return (Vector3D) { a.i - b.i, a.j - b.j, a.k - b.k };
}

double Vector3DDotProduct(const Vector3D a, const Vector3D b) {
    return a.i*b.i + a.j*b.j + a.k+b.k;
}

Vector3D Vector3DCrossProduct(const Vector3D a, const Vector3D b) {
    return (Vector3D) { a.j*b.k - a.k*b.j,
                        a.k*b.i - a.i*b.k,
                        a.i*b.j - a.j*b.i
    };
}

Vector3D Vector3DScalarProduct(double a, const Vector3D v) {
    return (Vector3D) { a * v.i, a * v.j, a * v.k };
}

double Vector3DMagnitude(const Vector3D v) {
    return sqrt(v.i*v.i + v.j*v.j + v.k*v.k);
}

Vector3D Vector3DNormalized(const Vector3D v) {
    if(Vector3DMagnitude(v) == 0) return v;

    return Vector3DScalarProduct(1/Vector3DMagnitude(v), v);
}

typedef struct Matrix3D {
    double a, b, c;
    double d, e, f;
    double g, h, i;
} Matrix3D;

#define Identity3D (Matrix3D) { 1,0,0, 0,1,0, 0,0,1 }
#define Zero3D (Matrix3D) { 0,0,0, 0,0,0, 0,0,0 }
#define Rotation3D(a,b,c) (Matrix3D) { cos(a)*cos(b), cos(a)*sin(b)*sin(c) - sin(a)*cos(c), cos(a)*sin(b)*cos(c) + sin(a)*sin(c),  sin(a)*cos(b), sin(a)*sin(b)*sin(c) + cos(a)*cos(c), sin(a)*sin(b)*cos(c) - cos(a)*sin(c), -sin(b), cos(b)*sin(c), cos(b)*cos(c)}

void Matrix3DPrint(const Matrix3D M) {
    printf("\n%f \t %f \t %f \n%f \t %f \t %f \n%f \t %f \t %f \n\n", M.a, M.b, M.c, M.d, M.e, M.f, M.g, M.h, M.i);
}

Matrix3D Matrix3DMatrixProduct (const Matrix3D A, const Matrix3D B) {
    return (Matrix3D) {
        A.a*B.a + A.b*B.d + A.c*B.g, A.a*B.b + A.b*B.e + A.c*B.h, A.a*B.c + A.b*B.f + A.c*B.i,
        A.d*B.a + A.e*B.d + A.f*B.g, A.d*B.b + A.e*B.e + A.f*B.h, A.d*B.c + A.e*B.f + A.f*B.i,
        A.g*B.a + A.h*B.d + A.i*B.g, A.g*B.b + A.h*B.e + A.i*B.h, A.g*B.c + A.h*B.f + A.i*B.i,
    };
}

Vector3D Matrix3DVectorProduct (const Matrix3D M, const Vector3D v) {
    return (Vector3D) {
        M.a*v.i + M.b*v.j + M.c*v.k,
        M.d*v.i + M.e*v.j + M.f*v.k,
        M.g*v.i + M.h*v.j + M.i*v.k,
    };
}

Matrix3D Matrix3DInverse (const Matrix3D M) {
    const double det = M.a*(M.e*M.i - M.f*M.h) - M.b*(M.d*M.i - M.f*M.g) + M.c*(M.d*M.h - M.e*M.g);

    return (Matrix3D) {
        +(M.e*M.i - M.f*M.h)/det,   -(M.b*M.i - M.c*M.h)/det,  +(M.b*M.f - M.c*M.e)/det,
        -(M.d*M.i - M.f*M.g)/det,   +(M.a*M.i - M.c*M.g)/det,  -(M.a*M.f - M.d*M.c)/det,
        +(M.d*M.h - M.e*M.g)/det,   -(M.a*M.h - M.b*M.g)/det,  +(M.a*M.e - M.b*M.d)/det,
    };
}

typedef struct Triangle {
    Vector3D a,b,c;
} Triangle;

// returns 1 is ray intersects triangle 0 otherwise
int TriangleRayIntersection(const Triangle triangle, const Vector3D start, const Vector3D direction, Vector3D *poi, double *z) {
    const Vector3D edge1 = Vector3DSubtraction(triangle.b, triangle.a);
    const Vector3D edge2 = Vector3DSubtraction(triangle.c, triangle.a);
    const Vector3D target = Vector3DSubtraction(start, triangle.a);

    if(Vector3DDotProduct(Vector3DCrossProduct(edge1, edge2), direction) < DBL_EPSILON) return 0;


    const Matrix3D mat = {
        -direction.i, edge1.i, edge2.i,
        -direction.j, edge1.j, edge2.j,
        -direction.k, edge1.k, edge2.k,
    };

    const Matrix3D inv = Matrix3DInverse(mat);
    const Vector3D v = Matrix3DVectorProduct(inv, target);

    const double t = v.i;
    const double x = v.j;
    const double y = v.k;

    if(poi != NULL) {
        *poi = Vector3DAddition(triangle.a, Vector3DAddition(
            Vector3DScalarProduct(x, edge1),
            Vector3DScalarProduct(y,edge2)
        ));
    }

    if(z != NULL) {
        *z = t;
    }

    return (t > 0 && x > 0 && y > 0 && (1 - x - y)>0);
}

Vector3D TriangleNormal(const Triangle triangle) {
    return Vector3DNormalized(Vector3DCrossProduct(
        Vector3DSubtraction(triangle.b, triangle.a),
        Vector3DSubtraction(triangle.c, triangle.a)
    ));
}

// takes a double from -1 to 1 and outputs a char with that light level
char LightnessToChar(const double l) {
    const int x = (l/2 + 1) * 11.9;
    if(x < 0) return '.';
    if(x > 12) return '@';
    return ".,-~:;=!*#$@"[x];
}

char TrianglesRasterize(const Triangle triangles[], size_t n, const Vector3D start, const Vector3D direction) {
    const Vector3D light = Vector3DNormalized((Vector3D){1,1,-1});

    double minZ = 0;
    double z;
    int mini = 0;

    for (int i = 0; i < n; ++i) {
        const Triangle triangle = triangles[i];

        const int col = TriangleRayIntersection(triangle, start, direction, NULL, &z);

        if(col==0) continue;
        if(z>minZ && minZ) continue;

        minZ = z;
        mini = i;
    }


    Triangle t = triangles[mini];
    Vector3D normal = TriangleNormal(t);

    if (minZ == 0) return ' ';
    return LightnessToChar(Vector3DDotProduct(normal, light));
}

// wackey windows api stuff
void ClearScreen() {
    COORD topLeft  = { 0, 0 };
    HANDLE console = GetStdHandle(STD_OUTPUT_HANDLE);
    CONSOLE_SCREEN_BUFFER_INFO screen;
    DWORD written;

    GetConsoleScreenBufferInfo(console, &screen);
    FillConsoleOutputCharacterA(
        console, ' ', screen.dwSize.X * screen.dwSize.Y, topLeft, &written
    );
    FillConsoleOutputAttribute(
        console, FOREGROUND_GREEN | FOREGROUND_RED | FOREGROUND_BLUE,
        screen.dwSize.X * screen.dwSize.Y, topLeft, &written
    );
    SetConsoleCursorPosition(console, topLeft);
}


#define HEIGHT 30
#define WIDTH 100

#define CAR 2.285714285714286


int main(void) {
    static double time = 0;

    const Vector3D A = {
        0.0,
        0.5,
        0
    };

    const Vector3D B = {
        0.58,
        -0.5,
        0
    };

    const Vector3D C = {
        -0.58,
        -0.50,
        0
    };

    const Vector3D D = {
        0,
        -0.167,
        1.054
    };

    const Triangle triangle1 = {A, B, C};
    const Triangle triangle2 = {B, C, D};
    const Triangle triangle3 = {C, D, A};
    const Triangle triangle4 = {D, A, B};

    const Triangle triangles[] = {triangle1, triangle2, triangle3, triangle4};

    while(1) {
        char buffer[2 * HEIGHT * WIDTH];
        int i = 0;
        for (int y = HEIGHT/2; y >= -HEIGHT/2; --y) {
            for (int x = -WIDTH/2; x <= WIDTH/2; ++x) {
                const Matrix3D R = Rotation3D(time, time*1.2, time*0.7);

                Triangle rotated_triangles[sizeof(triangles)/sizeof(Triangle)];
                for (int j = 0; j < sizeof(triangles)/sizeof(Triangle); ++j) {
                    rotated_triangles[j] = (Triangle) {
                        Matrix3DVectorProduct(R, triangles[j].a),
                        Matrix3DVectorProduct(R, triangles[j].b),
                        Matrix3DVectorProduct(R, triangles[j].c),
                    };
                }

                const char c = TrianglesRasterize(rotated_triangles, 3, (Vector3D){0,0,-10}, (Vector3D){x/CAR,y,120});
                buffer[i++] = c;
            }
            buffer[i++] = '\n';
        }

        // ClearScreen();
        fwrite(buffer,1,i,stdout);
        sleep(0.1);
        time += 0.005;
    }

    return 0;
}
