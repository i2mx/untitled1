#include <stdio.h>
#include <math.h>
#include <unistd.h>

typedef struct Vector3D {
    double i;
    double j;
    double k;
} Vector3D;

void Vector3DPrint(const Vector3D v) {
    printf("[%f, %f, %f]", v.i, v.j, v.k);
}

Vector3D Vector3DAddition(const Vector3D a, const Vector3D b) {
    return (Vector3D) {a.i + b.i, a.j + b.j, a.k + b.k};
}

Vector3D Vector3DSubtraction(const Vector3D a, const Vector3D b) {
    return (Vector3D) {a.i - b.i, a.j - b.j, a.k - b.k};
}

double Vector3DDotProduct(const Vector3D a, const Vector3D b) {
    return a.i * b.i + a.j * b.j + a.k + b.k;
}

Vector3D Vector3DCrossProduct(const Vector3D a, const Vector3D b) {
    return (Vector3D) {a.j * b.k - a.k * b.j,
                       a.k * b.i - a.i * b.k,
                       a.i * b.j - a.j * b.i
    };
}

Vector3D Vector3DScalarProduct(double a, const Vector3D v) {
    return (Vector3D) {a * v.i, a * v.j, a * v.k};
}

double Vector3DMagnitude(const Vector3D v) {
    return sqrt(v.i * v.i + v.j * v.j + v.k * v.k);
}

Vector3D Vector3DNormalized(const Vector3D v) {
    if (Vector3DMagnitude(v) == 0) return v;
    return Vector3DScalarProduct(1 / Vector3DMagnitude(v), v);
}

typedef struct Matrix3D {
    double a, b, c;
    double d, e, f;
    double g, h, i;
} Matrix3D;

#define Matrix3DIdentity (Matrix3D) { 1,0,0, 0,1,0, 0,0,1 }
#define Matrix3DZero (Matrix3D) { 0,0,0, 0,0,0, 0,0,0 }
#define Matrix3DRotation(a, b, c) (Matrix3D) {                                                      \
    cos(a)*cos(b),  cos(a)*sin(b)*sin(c) - sin(a)*cos(c), cos(a)*sin(b)*cos(c) + sin(a)*sin(c),     \
    sin(a)*cos(b),  sin(a)*sin(b)*sin(c) + cos(a)*cos(c),   sin(a)*sin(b)*cos(c) - cos(a)*sin(c),   \
    -sin(b),        cos(b)*sin(c),                          cos(b)*cos(c)                           \
}

void Matrix3DPrint(const Matrix3D M) {
    printf("\n%f \t %f \t %f \n%f \t %f \t %f \n%f \t %f \t %f \n\n", M.a, M.b, M.c, M.d, M.e, M.f, M.g, M.h, M.i);
}

Matrix3D Matrix3DMatrixProduct(const Matrix3D A, const Matrix3D B) {
    return (Matrix3D) {
            A.a * B.a + A.b * B.d + A.c * B.g, A.a * B.b + A.b * B.e + A.c * B.h, A.a * B.c + A.b * B.f + A.c * B.i,
            A.d * B.a + A.e * B.d + A.f * B.g, A.d * B.b + A.e * B.e + A.f * B.h, A.d * B.c + A.e * B.f + A.f * B.i,
            A.g * B.a + A.h * B.d + A.i * B.g, A.g * B.b + A.h * B.e + A.i * B.h, A.g * B.c + A.h * B.f + A.i * B.i,
    };
}

Vector3D Matrix3DVectorProduct(const Matrix3D M, const Vector3D v) {
    return (Vector3D) {
            M.a * v.i + M.b * v.j + M.c * v.k,
            M.d * v.i + M.e * v.j + M.f * v.k,
            M.g * v.i + M.h * v.j + M.i * v.k,
    };
}

Matrix3D Matrix3DInverse(const Matrix3D M) {
    const double det = M.a * (M.e * M.i - M.f * M.h) - M.b * (M.d * M.i - M.f * M.g) + M.c * (M.d * M.h - M.e * M.g);
    const double id = 1 / det;
    return (Matrix3D) {
            +(M.e * M.i - M.f * M.h) * id, -(M.b * M.i - M.c * M.h) * id, +(M.b * M.f - M.c * M.e) * id,
            -(M.d * M.i - M.f * M.g) * id, +(M.a * M.i - M.c * M.g) * id, -(M.a * M.f - M.d * M.c) * id,
            +(M.d * M.h - M.e * M.g) * id, -(M.a * M.h - M.b * M.g) * id, +(M.a * M.e - M.b * M.d) * id,
    };
}

typedef struct Triangle {
    Vector3D a, b, c;
} Triangle;

// returns 1 is ray intersects triangle 0 otherwise
int TriangleRayIntersection(const Triangle triangle, const Vector3D start, const Vector3D direction, Vector3D *poi,
                            double *z) {
    const Vector3D edge1 = Vector3DSubtraction(triangle.b, triangle.a);
    const Vector3D edge2 = Vector3DSubtraction(triangle.c, triangle.a);
    const Vector3D target = Vector3DSubtraction(start, triangle.a);

    if (fabs(Vector3DDotProduct(Vector3DCrossProduct(edge1, edge2), direction)) < 1e-8) return 0;

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

    if (poi != NULL) {
        *poi = Vector3DAddition(triangle.a, Vector3DAddition(
                Vector3DScalarProduct(x, edge1),
                Vector3DScalarProduct(y, edge2)
        ));
    }

    if (z != NULL) {
        *z = t;
    }

    return (t > 0 && x > 0 && y > 0 && (1 - x - y) > 0);
}

Vector3D TriangleNormal(const Triangle triangle) {
    return Vector3DNormalized(Vector3DCrossProduct(
            Vector3DSubtraction(triangle.b, triangle.a),
            Vector3DSubtraction(triangle.c, triangle.a)
    ));
}

// takes a double from 0 to 1 and outputs a char with that light level
char LightnessToChar(const double l) {
//    const double brightness = (l / 2 + 0.5);
    const int x = (int) (l * l * 12);
    if (x <= 0) return '.';
    if (x >= 12) return '@';
    return ".,-~:;=!*#$@"[x];
}

char TrianglesRaster(const Triangle triangles[], size_t n, const Vector3D start, const Vector3D direction) {
    const Vector3D light = {-4, 6, -10};

    double minZ = 0;
    double z;
    int mini = 0;
    Vector3D poi;
    Vector3D p;

    for (int i = 0; i < n; ++i) {
        const Triangle triangle = triangles[i];
        const int col = TriangleRayIntersection(triangle, start, direction, &p, &z);

        if (col == 0) continue;
        if (z < minZ && minZ) continue;

        minZ = z;
        mini = i;
        poi = p;
    }

    if (minZ == 0) return ' ';

    const Triangle t = triangles[mini];
    const Vector3D normal = TriangleNormal(t);

    const Vector3D lightDirection = Vector3DNormalized(Vector3DSubtraction(poi, light));
    const Vector3D viewDirection = Vector3DNormalized(direction);
    const Vector3D halfwayDirection = Vector3DNormalized(Vector3DAddition(lightDirection, viewDirection));

    const double d = Vector3DDotProduct(normal, Vector3DNormalized(lightDirection));
    const double diffuse = d / 2 + 0.5;
    const double alignment = Vector3DDotProduct(normal, halfwayDirection);
    const double speculative = pow((alignment > 0 ? alignment : 0), 0.5);

    const double l = diffuse * 0.8 + speculative * 0.3;
    return LightnessToChar(l);
}

#define HEIGHT 45
#define WIDTH 80
#define CAR 2.285714285714286

int main(void) {
    static double time = 0;

//    const Vector3D A = {0.0,0.5,0.0};
//    const Vector3D B = {0.58,-0.5,0.0};
//    const Vector3D C = {-0.58,-0.50,0.0};
//    const Vector3D D = {0.0,-0.167,1.054};
//
//    const Triangle triangle1 = {A, B, C};
//    const Triangle triangle2 = {D, B, C};
//    const Triangle triangle3 = {D, A, C};
//    const Triangle triangle4 = {D, B, A};
//    const Triangle triangles[] = {triangle1, triangle2, triangle3, triangle4};

//    Triangle triangles[12] = {
//            {{-0.5, -0.5, -0.5}, {0.5,  -0.5, -0.5}, {0.5,  0.5,  -0.5}},
//            {{0.5,  0.5,  -0.5}, {-0.5, 0.5,  -0.5}, {-0.5, -0.5, -0.5}},
//            {{-0.5, -0.5, 0.5},  {0.5,  -0.5, 0.5},  {0.5,  0.5,  0.5}},
//            {{0.5,  0.5,  0.5},  {-0.5, 0.5,  0.5},  {-0.5, -0.5, 0.5}},
//            {{-0.5, 0.5,  -0.5}, {0.5,  0.5,  -0.5}, {0.5,  0.5,  0.5}},
//            {{0.5,  0.5,  0.5},  {-0.5, 0.5,  0.5},  {-0.5, 0.5,  -0.5}},
//            {{-0.5, -0.5, -0.5}, {0.5,  -0.5, -0.5}, {0.5,  -0.5, 0.5}},
//            {{0.5,  -0.5, 0.5},  {-0.5, -0.5, 0.5},  {-0.5, -0.5, -0.5}},
//            {{-0.5, -0.5, -0.5}, {-0.5, -0.5, 0.5},  {-0.5, 0.5,  0.5}},
//            {{-0.5, 0.5,  0.5},  {-0.5, 0.5,  -0.5}, {-0.5, -0.5, -0.5}},
//            {{0.5,  -0.5, -0.5}, {0.5,  -0.5, 0.5},  {0.5,  0.5,  0.5}},
//            {{0.5,  0.5,  0.5},  {0.5,  0.5,  -0.5}, {0.5,  -0.5, -0.5}}
//    };

    const double phi = (1.0 + sqrt(5.0)) / 2.0;

    Vector3D vertices[] = {
            {0,    1,    phi},
            {0,    -1,   phi},
            {0,    1,    -phi},
            {0,    -1,   -phi},
            {phi,  0,    1},
            {phi,  0,    -1},
            {-phi, 0,    1},
            {-phi, 0,    -1},
            {1,    phi,  0},
            {-1,   phi,  0},
            {1,    -phi, 0},
            {-1,   -phi, 0},
    };

    Triangle triangles[] = {
            {vertices[0], vertices[8], vertices[4]},
            {vertices[0], vertices[4], vertices[1]},
            {vertices[4], vertices[10], vertices[1]},
            {vertices[4], vertices[5], vertices[10]},
            {vertices[5], vertices[4], vertices[8]},
            {vertices[2], vertices[5], vertices[8]},
            {vertices[2], vertices[8], vertices[9]},
            {vertices[8], vertices[9], vertices[2]},
            {vertices[9], vertices[8], vertices[0]},
            {vertices[9], vertices[0], vertices[6]},
            {vertices[6], vertices[0], vertices[1]},
            {vertices[6], vertices[1], vertices[11]},
            {vertices[11], vertices[1], vertices[10]},
            {vertices[11], vertices[10], vertices[3]},
            {vertices[3], vertices[10], vertices[5]},
            {vertices[3], vertices[5], vertices[2]},
            {vertices[3], vertices[2], vertices[7]},
            {vertices[7], vertices[2], vertices[9]},
            {vertices[7], vertices[9], vertices[6]},
            {vertices[7], vertices[6], vertices[11]},
            {vertices[7], vertices[11], vertices[3]},
    };

    while (1) {
        char buffer[(HEIGHT + 1) * (WIDTH + 2 + 1)];
        const Matrix3D R = Matrix3DRotation(0, time, time * 0.7);
        int i = 0;

        for (int y = HEIGHT / 2; y >= -HEIGHT / 2; --y) {
            for (int x = -WIDTH / 2; x <= WIDTH / 2; ++x) {
                size_t n = sizeof(triangles) / sizeof(Triangle);

                Triangle rotated_triangles[n];
                for (int j = 0; j < n; ++j) {
                    rotated_triangles[j] = (Triangle) {
                            Matrix3DVectorProduct(R, triangles[j].a),
                            Matrix3DVectorProduct(R, triangles[j].b),
                            Matrix3DVectorProduct(R, triangles[j].c),
                    };
                }

                const char c = TrianglesRaster(rotated_triangles, n,
                                               (Vector3D) {0, 0, -40},
                                               (Vector3D) {x / CAR, y, 350});
                buffer[i++] = c;
            }
            buffer[i++] = '|';
            buffer[i++] = '\n';
        }
        fwrite(buffer, 1, i, stdout);
        usleep(16666);
        time += 0.1;
    }
}

// TODO: Add inputs into this