

#include <vector>
#include <array>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <iostream>
#include <algorithm>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <Shellapi.h>
#endif

struct Vec3 { float x, y, z; };
struct Vec4 { float x, y, z, w; };
inline Vec3 operator+(const Vec3& a, const Vec3& b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
inline Vec3 operator-(const Vec3& v) { return { -v.x, -v.y, -v.z }; }   // ★ 추가
inline Vec3 operator*(const Vec3& a, float s) { return { a.x * s, a.y * s, a.z * s }; }
inline Vec3 operator/(const Vec3& a, float s) { return { a.x / s, a.y / s, a.z / s }; }
inline float dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }
inline float length(const Vec3& v) { return std::sqrt(dot(v, v)); }
inline Vec3 normalize(const Vec3& v) { return v / length(v); }
inline float clampf(float x, float lo = 0.f, float hi = 1.f) { return x < lo ? lo : (x > hi ? hi : x); }


struct Mat4 {
    float m[4][4]{};
    static Mat4 identity() { Mat4 I; for (int i = 0; i < 4; ++i) I.m[i][i] = 1; return I; }
    Vec4 operator*(const Vec4& v) const {
        Vec4 o{};
        o.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3] * v.w;
        o.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3] * v.w;
        o.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3] * v.w;
        o.w = m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3] * v.w;
        return o;
    }
    Mat4 operator*(const Mat4& b) const {
        Mat4 r{};
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int k = 0; k < 4; ++k)
                    r.m[i][j] += m[i][k] * b.m[k][j];
        return r;
    }
};

Mat4 translate(float tx, float ty, float tz) {
    Mat4 T = Mat4::identity();
    T.m[0][3] = tx; T.m[1][3] = ty; T.m[2][3] = tz; return T;
}
Mat4 scale(float sx, float sy, float sz) {
    Mat4 S = Mat4::identity();
    S.m[0][0] = sx; S.m[1][1] = sy; S.m[2][2] = sz; return S;
}
Mat4 perspective(float l, float r, float b, float t, float n, float f) {
    Mat4 P{};
    P.m[0][0] = 2 * n / (r - l);
    P.m[1][1] = 2 * n / (t - b);
    P.m[0][2] = (r + l) / (r - l);
    P.m[1][2] = (t + b) / (t - b);
    P.m[2][2] = -(f + n) / (f - n);
    P.m[2][3] = -2 * f * n / (f - n);
    P.m[3][2] = -1;
    return P;
}
Mat4 viewport(int nx, int ny) {
    Mat4 V = Mat4::identity();
    V.m[0][0] = nx / 2.0f;  V.m[0][3] = (nx - 1) / 2.0f;
    V.m[1][1] = -ny / 2.0f;  V.m[1][3] = (ny - 1) / 2.0f; // Flip Y
    V.m[2][2] = 0.5f;       V.m[2][3] = 0.5f;
    return V;
}

// -------------------------- Scene generation ---------------------------
int               gNumVertices = 0, gNumTriangles = 0;
std::vector<Vec3> gVertexBuffer;
int* gIndexBuffer = nullptr;
std::vector<Vec3> gVertexColor;    // per‑vertex linear RGB (for Gouraud)

void create_scene() {
    const int   W = 32, H = 16;
    const float PI = 3.14159265358979323846f;
    gNumVertices = (H - 2) * W + 2;
    gNumTriangles = (H - 2) * (W - 1) * 2;
    gVertexBuffer.resize(gNumVertices);
    gVertexColor.resize(gNumVertices);
    gIndexBuffer = new int[3 * gNumTriangles];

    int id = 0;
    for (int j = 1; j < H - 1; ++j)
        for (int i = 0; i < W; ++i, ++id) {
            float theta = float(j) / (H - 1) * PI;
            float phi = float(i) / (W - 1) * 2 * PI;
            gVertexBuffer[id] = { std::sin(theta) * std::cos(phi),
                                   std::cos(theta),
                                  -std::sin(theta) * std::sin(phi) };
        }
    gVertexBuffer[id++] = { 0,  1, 0 }; // north pole
    gVertexBuffer[id++] = { 0, -1, 0 }; // south pole

    id = 0;
    for (int j = 0; j < H - 3; ++j)
        for (int i = 0; i < W - 1; ++i) {
            gIndexBuffer[id++] = j * W + i;
            gIndexBuffer[id++] = (j + 1) * W + i + 1;
            gIndexBuffer[id++] = j * W + i + 1;

            gIndexBuffer[id++] = j * W + i;
            gIndexBuffer[id++] = (j + 1) * W + i;
            gIndexBuffer[id++] = (j + 1) * W + i + 1;
        }
    for (int i = 0; i < W - 1; ++i) {
        gIndexBuffer[id++] = (H - 2) * W;     gIndexBuffer[id++] = i;                 gIndexBuffer[id++] = i + 1;
        gIndexBuffer[id++] = (H - 2) * W + 1; gIndexBuffer[id++] = (H - 3) * W + i + 1; gIndexBuffer[id++] = (H - 3) * W + i;
    }
}

// 2‑D edge function (counter‑clockwise > 0)
inline float edge(const Vec3& a, const Vec3& b, const Vec3& c) {
    return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
}
inline bool inTri(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    float w0 = edge(b, c, p), w1 = edge(c, a, p), w2 = edge(a, b, p);
    return (w0 >= 0 && w1 >= 0 && w2 >= 0) || (w0 <= 0 && w1 <= 0 && w2 <= 0);
}

void saveBMP(const char* fname, int W, int H, const std::vector<std::array<unsigned char, 3>>& pix) {
#ifdef _WIN32
    int rowPitch = ((W * 3 + 3) & ~3);
    int dataSize = rowPitch * H;

    BITMAPFILEHEADER fh{}; BITMAPINFOHEADER ih{};
    fh.bfType = 0x4D42; // "BM"
    fh.bfOffBits = sizeof(fh) + sizeof(ih);
    fh.bfSize = fh.bfOffBits + dataSize;

    ih.biSize = sizeof(ih);
    ih.biWidth = W;  ih.biHeight = H;
    ih.biPlanes = 1;  ih.biBitCount = 24;
    ih.biCompression = BI_RGB;

    std::ofstream ofs(fname, std::ios::binary);
    ofs.write((char*)&fh, sizeof(fh));
    ofs.write((char*)&ih, sizeof(ih));

    std::vector<unsigned char> row(rowPitch);
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            const auto& c = pix[(H - 1 - y) * W + x];
            row[x * 3 + 0] = c[2]; row[x * 3 + 1] = c[1]; row[x * 3 + 2] = c[0];
        }
        ofs.write((char*)row.data(), rowPitch);
    }
#else
    (void)fname; (void)W; (void)H; (void)pix;
#endif
}

// --------------------------- Main --------------------------------------
int main() {
    constexpr int NX = 512, NY = 512;
    create_scene();

    // Transformation matrices
    Mat4 M = translate(0, 0, -7) * scale(2, 2, 2); // world
    Mat4 V = Mat4::identity();                     // view (eye at origin)
    Mat4 P = perspective(-0.1f, 0.1f, -0.1f, 0.1f, -0.1f, -1000.f);
    Mat4 C = viewport(NX, NY) * P * V * M;         // combined

    // Material & light constants
    const Vec3 ka{ 0, 1, 0 }, kd{ 0, 0.5f, 0 }, ks{ 0.5f, 0.5f, 0.5f };
    const float Ia = 0.2f, shininess = 32.f;
    const Vec3  lightPos{ -4, 4, -3 };
    const Vec3  eyePos{ 0, 0, 0 };

    // ----------------------------------------------------------------- //
    // 1)   Per‑vertex lighting (linear RGB)
    // ----------------------------------------------------------------- //
    for (int i = 0; i < gNumVertices; ++i) {
        Vec4 w4 = M * Vec4{ gVertexBuffer[i].x, gVertexBuffer[i].y, gVertexBuffer[i].z, 1 };
        Vec3 w = { w4.x, w4.y, w4.z };

        // Normal: sphere scaled uniformly → just original position direction
        Vec3 N = normalize(Vec3{ gVertexBuffer[i].x, gVertexBuffer[i].y, gVertexBuffer[i].z });
        if (dot(N, eyePos - w) < 0) N = N * -1.f;  // face towards camera

        Vec3 L = normalize(lightPos - w);
        float diff = std::max(0.f, dot(N, L));

        Vec3 R = normalize(-L + N * (2.f * dot(N, L)));
        Vec3 Vv = normalize(eyePos - w);
        float spec = std::pow(std::max(0.f, dot(R, Vv)), shininess);

        gVertexColor[i] = ka * Ia + kd * diff + ks * spec; // linear RGB
    }

    // Frame‑buffer
    std::vector<std::array<unsigned char, 3>> color(NX * NY, { 0,0,0 });
    std::vector<float> depth(NX * NY, FLT_MAX);

    // ----------------------------------------------------------------- //
    // 2)   Triangle traversal & rasterisation (Gouraud shading)
    // ----------------------------------------------------------------- //
    for (int t = 0; t < gNumTriangles; ++t) {
        int ia = gIndexBuffer[3 * t], ib = gIndexBuffer[3 * t + 1], ic = gIndexBuffer[3 * t + 2];

        // Screen‑space vertices (after perspective divide)
        Vec4 p4[3] = {
            C * Vec4{ gVertexBuffer[ia].x, gVertexBuffer[ia].y, gVertexBuffer[ia].z, 1 },
            C * Vec4{ gVertexBuffer[ib].x, gVertexBuffer[ib].y, gVertexBuffer[ib].z, 1 },
            C * Vec4{ gVertexBuffer[ic].x, gVertexBuffer[ic].y, gVertexBuffer[ic].z, 1 }
        };
        Vec3 s[3];
        for (int i = 0; i < 3; ++i) s[i] = { p4[i].x / p4[i].w, p4[i].y / p4[i].w, p4[i].z / p4[i].w };

        // Vertex colours (linear space)
        Vec3 cLinV[3] = { gVertexColor[ia], gVertexColor[ib], gVertexColor[ic] };

        // Bounding box in screen space
        int xmin = std::max(0, std::min({ int(s[0].x), int(s[1].x), int(s[2].x) }));
        int xmax = std::min(NX - 1, std::max({ int(s[0].x), int(s[1].x), int(s[2].x) }));
        int ymin = std::max(0, std::min({ int(s[0].y), int(s[1].y), int(s[2].y) }));
        int ymax = std::min(NY - 1, std::max({ int(s[0].y), int(s[1].y), int(s[2].y) }));

        float denom = edge(s[1], s[2], s[0]); if (denom == 0) continue; // degenerate

        // Raster loop
        for (int y = ymin; y <= ymax; ++y) {
            for (int x = xmin; x <= xmax; ++x) {
                Vec3 p{ float(x) + 0.5f, float(y) + 0.5f, 0 };
                if (!inTri(p, s[0], s[1], s[2])) continue;

                float w0 = edge(s[1], s[2], p) / denom;
                float w1 = edge(s[2], s[0], p) / denom;
                float w2 = 1.0f - w0 - w1; // or edge(s[0], s[1], p)/denom

                float invW[3] = { 1 / p4[0].w, 1 / p4[1].w, 1 / p4[2].w };
                float w0p = w0 * invW[0], w1p = w1 * invW[1], w2p = w2 * invW[2];
                float sum = w0p + w1p + w2p;
                w0p /= sum; w1p /= sum; w2p /= sum;
                // Interpolated depth & colour (linear)
                float z = w0p * s[0].z + w1p * s[1].z + w2p * s[2].z;
                int idx = y * NX + x;
                if (z >= depth[idx]) continue;

                depth[idx] = z;
                Vec3 cLin = cLinV[0] * w0p + cLinV[1] * w1p + cLinV[2] * w2p;

                Vec3 cGam{ std::pow(clampf(cLin.x), 1 / 2.2f),
                           std::pow(clampf(cLin.y), 1 / 2.2f),
                           std::pow(clampf(cLin.z), 1 / 2.2f) };
                color[idx] = { (unsigned char)std::round(clampf(cGam.x) * 255.f),
                                (unsigned char)std::round(clampf(cGam.y) * 255.f),
                                (unsigned char)std::round(clampf(cGam.z) * 255.f) };
            }
        }
    }

    // ----------------------------------------------------------------- //
    // 3)   Save outputs
    // ----------------------------------------------------------------- //
    std::ofstream ppm("Q2_output.ppm", std::ios::binary);
    ppm << "P6\n" << NX << " " << NY << "\n255\n";
    for (auto& c : color) ppm.write((char*)c.data(), 3);

    saveBMP("Q2_output.bmp", NX, NY, color);
#ifdef _WIN32
    ShellExecuteA(nullptr, "open", "Q2_output.bmp", nullptr, nullptr, SW_SHOWNORMAL);
#endif
    return 0;
}
