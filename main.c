#include <float.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <windows.h>

#define EPSILON 1e-8

typedef struct Vector3D {
    double i;
    double j;
    double k;
} Vector3D;

void Vector3DPrint(const Vector3D v) {
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
    const double invdet = 1/det;
    return (Matrix3D) {
        +(M.e*M.i - M.f*M.h)*invdet,   -(M.b*M.i - M.c*M.h)*invdet,  +(M.b*M.f - M.c*M.e)*invdet,
        -(M.d*M.i - M.f*M.g)*invdet,   +(M.a*M.i - M.c*M.g)*invdet,  -(M.a*M.f - M.d*M.c)*invdet,
        +(M.d*M.h - M.e*M.g)*invdet,   -(M.a*M.h - M.b*M.g)*invdet,  +(M.a*M.e - M.b*M.d)*invdet,
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

    if(fabs(Vector3DDotProduct(Vector3DCrossProduct(edge1, edge2), direction)) < EPSILON) return 0;


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

    // const Vector3D p = Vector3DCrossProduct(direction, edge2);
    // const Vector3D q = Vector3DCrossProduct(target, edge1);
    //
    // const double det = Vector3DDotProduct(p, edge1);
    //
    // // culling is done here if we just check the raw value
    // if(fabs(det) < DBL_EPSILON) return 0;
    //
    //
    // const double invdet = 1/det;
    //
    // const double t = Vector3DDotProduct(q, edge2) * invdet;
    // const double x = Vector3DDotProduct(p, target) * invdet;
    // const double y = Vector3DDotProduct(q, direction) * invdet;

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
    const double brightness = (l/2 + 0.5);

    const int x = brightness * 12;


    if(x <= 0) return '.';
    if(x >= 12) return '@';
    // return '.';

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
        if(z<minZ && minZ) continue;

        minZ = z;
        mini = i;
    }


    const Triangle t = triangles[mini];
    const Vector3D normal = TriangleNormal(t);

    if (minZ == 0) return ' ';
    // return '0'+mini;
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
        char buffer[10 * HEIGHT * WIDTH];
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
                const char c = TrianglesRasterize(rotated_triangles, 4, (Vector3D){0,0,-10}, (Vector3D){x/CAR,y,120});
                // here incase you want to color code
                //
                // switch (c) {
                //     case '0':
                //         strcpy(&buffer[i], "\x1b[31m");
                //         // printf("\x1b[31m");
                //     break;
                //     case '1':
                //         strcpy(&buffer[i], "\x1b[32m");
                //         // printf("\x1b[32m");
                //     break;
                //     case '2':
                //         strcpy(&buffer[i], "\x1b[34m");
                //         // printf("\x1b[34m");
                //     break;
                //     case '3':
                //         strcpy(&buffer[i], "\x1b[0m");
                //     // printf("\x1b[34m");
                //     break;
                //     default:
                //         i -= 5;
                //         // strcat("\x1b[0m", &buffer[i]);
                //         // printf("\x1b[0m");
                //     break;


                // }

                // printf("%c",c);

                // i += 5;
                buffer[i++] = c;
            }

            // printf("\n");
            buffer[i++] = '\n';
        }

        // ClearScreen();
        fwrite(buffer,1,i,stdout);
        sleep(0.01);
        time += 0.01;
    }

    return 0;
}
