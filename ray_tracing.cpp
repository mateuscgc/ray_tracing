#include <bits/stdc++.h>
#include <cassert>
#include <GL/glew.h>
#include <glm/vec3.hpp> 
#include <glm/vec4.hpp> 
#include <glm/mat4x4.hpp> 
#include <glm/gtc/matrix_transform.hpp> 
#include <glm/gtc/type_ptr.hpp>
#include <glm/ext.hpp>

#if defined (__APPLE__) || defined(MACOSX)
    #include <GLUT/glut.h>
#else
    #include <GL/glut.h>
#endif

#include "myDataStructures.h"
#include "initShaders.h"

using namespace std;

typedef vector< vector<double> > matrix;
#define point(x) matrix x(4, vector<double>(1))
#define transformation_matrix(x) matrix x(4, vector<double>(4))
#define PI acos(-1)

int SCREEN_WIDTH = 600;
int SCREEN_HEIGHT = 600;

int MAX_RECURSION = 1;

GLuint  mainShader;
GLuint  axisVBO[5];

double X(const matrix& m) { return m[0][0]; }
double Y(const matrix& m) { return m[1][0]; }
double Z(const matrix& m) { return m[2][0]; }

void set_x(matrix& m, double x) { m[0][0] = x; }
void set_y(matrix& m, double y) { m[1][0] = y; }
void set_z(matrix& m, double z) { m[2][0] = z; }

void set_point(matrix& m, double x, double y, double z) {
    m = matrix(4, vector<double>(1));
    set_x(m, x);
    set_y(m, y);
    set_z(m, z);
    m[3][0] = 1;
}

// struct point3D {
//     double x, y, z;
//     point3D(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}
//     point3D& operator-=(const point3D& rhs) {
//         this->x -= rhs.x;
//         this->y -= rhs.y;
//         this->z -= rhs.z;
//         return *this;
//     }
//     const point3D operator-(const point3D& rhs) const {
//         this->x -= rhs.x;
//         this->y -= rhs.y;
//         this->z -= rhs.z;
//         return *this;
//     }
// };

struct Color {
    double r, g, b;
    Color() {}
    Color(double rr, double g, double bb) : r(rr), g(g), b(bb) {}
    Color& operator*=(const Color& rhs) {
        this->r *= rhs.r;
        this->g *= rhs.g;
        this->b *= rhs.b;
        return *this;
    }
    Color& operator*=(double rhs) {
        this->r *= rhs;
        this->g *= rhs;
        this->b *= rhs;
        return *this;
    }
    Color& operator+=(const Color& rhs) {
        this->r += rhs.r;
        this->g += rhs.g;
        this->b += rhs.b;
        return *this;
    }
    const Color operator*(const Color& rhs) {
        this->r *= rhs.r;
        this->g *= rhs.g;
        this->b *= rhs.b;
        return *this;
    }
    const Color operator*(const double rhs) {
        this->r *= rhs;
        this->g *= rhs;
        this->b *= rhs;
        return *this;
    }
    const Color operator-(const Color& rhs) {
        this->r -= rhs.r;
        this->g -= rhs.g;
        this->b -= rhs.b;
        return *this;
    }
    const Color operator+(const Color& rhs) {
        this->r += rhs.r;
        this->g += rhs.g;
        this->b += rhs.b;
        return *this;
    }
    void clamp(double l, double h) {
        r = max(l, min(r, h));
        g = max(l, min(g, h));
        b = max(l, min(b, h));
    }
};

inline Color operator*(double rhs, Color lhs) {
  lhs.r *= rhs;
  lhs.g *= rhs;
  lhs.b *= rhs;
  return lhs;
}

matrix subtractmat(const matrix& a, const matrix& b) {
    matrix ans(a.size(), vector<double>(a[0].size()));
    for(int i = 0; i < ans.size(); i++) {
        for(int j = 0; j < ans[0].size(); j++) {
            ans[i][j] = a[i][j] - b[i][j];
        }
    }
    return ans;
}

matrix summat(const matrix& a, const matrix& b) {
    matrix ans(a.size(), vector<double>(a[0].size()));
    for(int i = 0; i < ans.size(); i++) {
        for(int j = 0; j < ans[0].size(); j++) {
            ans[i][j] = a[i][j] + b[i][j];
        }
    }
    return ans;
}

matrix multiconst(const matrix& a, double b) {
    matrix ans(a.size(), vector<double>(a[0].size()));
    for(int i = 0; i < ans.size(); i++) {
        for(int j = 0; j < ans[0].size(); j++) {
            ans[i][j] = a[i][j] * b;
        }
    }
    return ans;
}

matrix multimat(const matrix& a, const matrix& b) {
    matrix ans(a.size(), vector<double>(b[0].size()));
    for(int i = 0; i < ans.size(); i++) {
        for(int j = 0; j < ans[0].size(); j++) {
            for(int l = 0; l < b.size(); l++) {
                ans[i][j] += (a[i][l]*b[l][j]);
            }
        }
    }
    return ans;
}

double magnitude(const matrix& m) {
    return sqrt(X(m)*X(m) + Y(m)*Y(m) + Z(m)*Z(m));
}

matrix scale_matrix(double rate) {
    transformation_matrix(ans);
    ans[0][0] = rate;
    ans[1][1] = rate;
    ans[2][2] = rate;
    ans[3][3] = 1;
    return ans;
}

matrix normalize(const matrix& m) {
    return multimat(m, scale_matrix(1/magnitude(m)));
}

double dot(const matrix& a, const matrix& b) {
    double ans = 0;
    for(int i = 0; i < 3; i++) {
        ans += a[i][0] * b[i][0];
        cout << "dot " << ans << endl;
    }
    return ans;
}

matrix reflect(const matrix& I, const matrix& N) {
    return subtractmat(I, multiconst(N, 2*dot(N, I)));
}

matrix negative(const matrix& a) {
    matrix ans(a.size(), vector<double>(a[0].size()));
    for(int i = 0; i < ans.size(); i++) {
        for(int j = 0; j < ans[0].size(); j++) {
            ans[i][j] = a[i][j] * -1;
        }
    }
    return ans;
}

matrix clamp(const matrix& a, double l, double h) {
    matrix ans(a.size(), vector<double>(a[0].size()));
    for(int i = 0; i < ans.size(); i++) {
        for(int j = 0; j < ans[0].size(); j++) {
            ans[i][j] = max(l, min(a[i][j], h));
        }
    }
    return ans;
}

struct Sphere {
    // matrix position;
    glm::dvec3 position;
    double radius;
    int id;
    double direct_phong_rate;
    double reflection_rate;
    double transmition_rate;
    glm::dvec3 color;
    // Color color;

    Sphere() {}
    static double a(const glm::dvec3& ray) {
        return ray.x*ray.x + ray.y*ray.y + ray.z*ray.z;
    }

    static double b(const glm::dvec3& ray, const glm::dvec3& ori) {
        return 2*(ray.x*ori.x + ray.y*ori.y + ray.z*ori.z);
    }

    static double c(const glm::dvec3& ori, double radius) {
        return ori.x*ori.x + ori.y*ori.y + ori.z*ori.z - radius*radius;
    }

    static double discriminant(const glm::dvec3& ray, const glm::dvec3& ori, double radius) {
        double a = Sphere::a(ray);
        double b = Sphere::b(ray, ori);
        double c = Sphere::c(ori, radius);
        return b*b - 4*a*c;
    }

    bool intersect(const glm::dvec3& ray, const glm::dvec3& ori, const glm::mat4 view_matrix) {
        // cout << X(multimat(view_matrix, position)) << " " << Y(multimat(view_matrix, position)) << " " << Z(multimat(view_matrix, position)) << endl;
        // cout << X(ori) << " " << Y(ori) << " " << Z(ori) << endl;
        
        // Translate ray origin as sphere whould be on 0,0,0
        // ori = subtractmat(ori, multimat(view_matrix, position));
        glm::dvec3 new_ori = ori - glm::dvec3(view_matrix * glm::dvec4(position, 1.0));
        
        // cout << X(ori) << " " << Y(ori) << " " << Z(ori) << endl;
        // cout << X(ray) << " " << Y(ray) << " " << Z(ray) << endl;
        // cout << radius << endl;
        // cout << Sphere::discriminant(ray, ori, radius) << endl;
        return Sphere::discriminant(ray, new_ori, radius) >= 0;
    }
    pair<double, double> intersections_time(const glm::dvec3& ray, const glm::dvec3& ori, const glm::mat4& view_matrix) {
        // ori = subtractmat(ori, multimat(view_matrix, position));
        glm::dvec3 new_ori = ori - glm::dvec3(view_matrix * glm::dvec4(position, 1.0));
        double a = Sphere::a(ray);
        double b = Sphere::b(ray, new_ori);
        double root1 = (-b + sqrt(Sphere::discriminant(ray, new_ori, radius)))/(2*a);
        double root2 = (-b - sqrt(Sphere::discriminant(ray, new_ori, radius)))/(2*a);
        // cout << root1 << " " << root2 << endl;
        return make_pair(root1, root2);
    }
};

vector<Sphere> spheres;

struct Light {
    // matrix position;
    glm::dvec3 position;
    glm::dvec3 color;
    // Color color;
    int id;
};

vector<Light> lights;

struct Camera {
    glm::dvec3 position;
    glm::dvec3 target;
    glm::dvec3 up_vec;
    glm::mat4 view_matrix;
    // matrix view;

    double field_of_view;
    double focal_length;
};

Camera camera;

pair<Sphere*, double> get_closest_intersection(const glm::dvec3& ray, const glm::dvec3& ori, const glm::mat4& view_matrix, Sphere *avoid) {
    Sphere *ans = nullptr;
    double best_time;
    for(Sphere &s : spheres) {
        if((avoid == nullptr || s.id != avoid->id) && s.intersect(ray, ori, view_matrix)) {
            pair<double, double> intersections = s.intersections_time(ray, ori, view_matrix);

            // Detect if camera is after or inside sphere
            double best = min(intersections.first, intersections.second);
            if(best < 0) continue;
            // cout << "intersection" << endl;

            if(ans == nullptr || best < best_time) {
                best_time = best;
                ans = &s;
            }
        }
    }
    return make_pair(ans, best_time);
}

void create_scene_objects() {
    int id = -1;
    Sphere s;

    s.position = glm::dvec3(-0.5,0.3,5.0);
    // set_point(s.position, 0, 0, 3.0);
    s.radius = 0.5;
    s.id = ++id;
    s.direct_phong_rate = 0.2;
    s.reflection_rate = 0.2;
    s.transmition_rate = 0.6;
    s.color = glm::dvec3(0.7,0.4,0.2);
    // s.color = Color(0.7,0.4,0.2);
    spheres.push_back(s);

    // s.position = glm::dvec3(0.2,0,7.0);
    // // set_point(s.position, 0, 0, 3.0);
    // s.radius = 0.5;
    // s.id = ++id;
    // s.direct_phong_rate = 0.5;
    // s.reflection_rate = 0.3;
    // s.transmition_rate = 0.2;
    // s.color = glm::dvec3(0.2,0.2,0.6);
    // // s.color = Color(0.7,0.4,0.2);
    // spheres.push_back(s);

    // s.position = glm::dvec3(0.2,1.5,7.0);
    // // set_point(s.position, 0, 0, 3.0);
    // s.radius = 0.5;
    // s.id = ++id;
    // s.direct_phong_rate = 0.5;
    // s.reflection_rate = 0.3;
    // s.transmition_rate = 0.2;
    // s.color = glm::dvec3(0.8,0.3,0.2);
    // // s.color = Color(0.7,0.4,0.2);
    // spheres.push_back(s);

    // s.position = glm::dvec3(0.6,1.5,4.0);
    // // set_point(s.position, 0, 0, 3.0);
    // s.radius = 0.5;
    // s.id = ++id;
    // s.direct_phong_rate = 0.5;
    // s.reflection_rate = 0.3;
    // s.transmition_rate = 0.2;
    // s.color = glm::dvec3(0.4,0.7,0.3);
    // // s.color = Color(0.7,0.4,0.2);
    // spheres.push_back(s);

    // ============================================

    Light l;

    l.position = glm::dvec3(0,3.0,0);
    // set_point(l.position, 0.0, 3.0, 0.0);
    l.color = glm::dvec3(1,1,1);
    // l.color = Color(1,1,1);
    l.id = ++id;   
    lights.push_back(l);


    // l.position = glm::dvec3(-3,2.0,4);
    // // set_point(l.position, 0.0, 3.0, 0.0);
    // l.color = glm::dvec3(1,1,1);
    // // l.color = Color(1,1,1);
    // l.id = ++id;   
    // lights.push_back(l);
 
    // ============================================


    // set_point(s.position, 2, 1, 4);
    // s.radius = 0.3;
    // s.id = ++id;   
    // s.direct_phong_rate = 0.5;
    // s.reflection_rate = 0.3;
    // s.transmition_rate = 0.2;
    // s.color = Color(0.7,0.4,0.2);
    // spheres.push_back(s);
    
    // set_point(s.position, 2.2, 0.4, 2.7);
    // s.radius = 1;
    // s.id = ++id;
    // s.direct_phong_rate = 0.2;
    // s.reflection_rate = 0.8;
    // s.transmition_rate = 0.0;
    // s.color = Color(0.0,0.0,0.8);
    // spheres.push_back(s);
 
    // set_point(s.position, 0.8, 1.2, 1.5);
    // s.radius = 0.8;
    // s.id = ++id;
    // s.direct_phong_rate = 0.2;
    // s.reflection_rate = 0.2;
    // s.transmition_rate = 0.6;
    // s.color = Color(0.6,0.8,0.0);
    // spheres.push_back(s);
 
    // set_point(s.position, 3.4, 2.1, 3.4);
    // s.radius = 0.5;
    // s.id = ++id;
    // s.direct_phong_rate = 0.1;
    // s.reflection_rate = 0.5;
    // s.transmition_rate = 0.4;
    // s.color = Color(0.0,0.5,0.5);
    // spheres.push_back(s);

    // set_point(s.position, 0, 0, 0);
    // s.radius = 1;
    // s.id = ++id;
    // s.direct_phong_rate = 0.8;
    // s.reflection_rate = 0.2;
    // s.transmition_rate = 0.0;
    // s.color = Color(1.0,1.0,1.0);
    // spheres.push_back(s);

    // set_point(l.position, 4.0, 5.0, 3.0);
    // l.color = Color(1,1,1);
    // s.id = ++id;
    // lights.push_back(l);

    // set_point(camera.position, 3, 3, 3);
    // set_point(camera.target, 1.5, 0.5, 1.5);
    // set_point(camera.up_vec, 0, 1, 0);

    camera.position = glm::dvec3(0,0,0);
    camera.target = glm::dvec3(0,0,-1);

    // camera.position = glm::dvec3(0,1,2);
    // camera.target = glm::dvec3(0,2,-1);

    // camera.position = glm::dvec3(-1,1.8,2);
    // camera.target = glm::dvec3(-2,3,-2);
    
    camera.up_vec = glm::dvec3(0,1,0);
    camera.view_matrix = glm::lookAt(camera.position, camera.target, camera.up_vec);
    // camera.view = matrix(4, vector<double>(4));

    // const float *pSource = (const float*)glm::value_ptr(camera.view_matrix);
    // for(int i = 0; i < 4; i++) {
    //     for(int j = 0; j < 4; j++) {
    //         camera.view[i][j] = pSource[i*4 + j];
    //         cout << camera.view[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    camera.field_of_view = PI/4; // 45º
    camera.focal_length = 1;
}

glm::dvec3 shade(const glm::dvec3& ray, const glm::dvec3& ray_ori, double intersection_time, Sphere *obj, int recursion_level) {
    if(obj == nullptr) {
        // Background color
        // return Color(0,0,0);
        return glm::dvec3(0,0,0);
    }
    // cout << obj->color.r << endl;

    glm::dvec3 result;

    glm::dvec3 intersection_point = ray_ori + glm::normalize(ray) * intersection_time;
    // cout << "int point " << glm::to_string(intersection_point) << endl;
    glm::dvec3 intersection_normal = intersection_point - glm::dvec3(camera.view_matrix * glm::dvec4(obj->position, 1.0));
    // cout << "int normal " << glm::to_string(intersection_normal) << endl;
    intersection_normal = glm::normalize(intersection_normal);

    glm::dvec3 matAmb = obj->color;
    glm::dvec3 matDif = obj->color;
    glm::dvec3 matSpec = obj->color;

    for(Light l : lights) {
        glm::dvec3 light_ray = glm::normalize(intersection_point - glm::dvec3(camera.view_matrix * glm::dvec4(l.position, 1.0)));
        // cout << "light ray " << glm::to_string(-light_ray) << endl;
        
        // Ambiente
        // Put ambient even with oclusion
        glm::dvec3 ambient = l.color * matAmb;
        // cout << "ambient " << glm::to_string(ambient) << endl;

        // Check if light oclusion
        // Need to check oclusion from obj itself
        pair<Sphere*, double> intersection = get_closest_intersection(light_ray, glm::dvec3(camera.view_matrix * glm::dvec4(l.position, 1.0)), camera.view_matrix, nullptr);
        // if(intersection.first != nullptr) cout << "whhat: " << intersection.first->id << " " << obj->id << endl;
        if(intersection.first != nullptr && intersection.first->id != obj->id) {
            result += obj->direct_phong_rate*ambient;
            continue;
        }

        // Diffuse
        double teta = max(glm::dot(-light_ray, intersection_normal), 0.0);
        glm::dvec3 diffuse = l.color * matDif * teta;
        // cout << "teta " << teta << " diffuse " << glm::to_string(diffuse) << endl;

        // Specular
        glm::dvec3 intersection_to_eye_vec = glm::normalize(ray_ori-intersection_point);
        glm::dvec3 light_reflected_ray = glm::normalize(glm::reflect(light_ray, intersection_normal));
        double omega = max(glm::dot(intersection_to_eye_vec, light_reflected_ray), 0.0);
        glm::dvec3 specular = l.color * matSpec * pow(omega, 2);
        // cout << "omega " << omega << " specular " << glm::to_string(specular) << endl;

        glm::dvec3 light_result = ambient + diffuse + specular;
        glm::clamp(light_result, 0.0, 1.0);
        result += obj->direct_phong_rate*light_result;
    }

    // cout << "-> " << result.r << " " << result.g << " " << result.b << endl;

    // result.clamp(0, 1);
    glm::clamp(result, 0.0, 1.0);

    if(recursion_level < MAX_RECURSION) {
        glm::dvec3 reflected_ray = glm::normalize(glm::reflect(ray, intersection_normal));
        pair<Sphere*, double> intersection = get_closest_intersection(reflected_ray, intersection_point, camera.view_matrix, nullptr);
        result += obj->reflection_rate * shade(reflected_ray, intersection_point, intersection.second, intersection.first, recursion_level+1);

        glm::dvec3 transmited_ray = glm::refract(ray, intersection_normal, 0.3); // Do some distortion
        // Maybe add obj diameter to ray origin or something more precise to avoid get the intersection with the same object (from inside)
        // Or get intersection thats is not with the same object
        intersection = get_closest_intersection(transmited_ray, intersection_point, camera.view_matrix, obj);
        result += obj->transmition_rate * shade(transmited_ray, intersection_point, intersection.second, intersection.first, recursion_level+1);
    }

    // result.clamp(0, 1);
    glm::clamp(result, 0.0, 1.0);

    return result;
}

vector<float> ray_tracing() {
    vector<float> image;

    double aspect_ratio = (double)SCREEN_WIDTH / (double)SCREEN_HEIGHT;

    double width_field_of_view = camera.field_of_view * aspect_ratio;
    double ini_x = -(tan(width_field_of_view/2)*camera.focal_length);
    double end_x = ini_x * -1;
    cout << ini_x << " " << end_x << endl;

    double height_field_of_view = camera.field_of_view;
    double ini_y = -(tan(height_field_of_view/2)*camera.focal_length);
    double end_y = ini_y * -1;
    cout << ini_y << " " << end_y << endl;

    double pixel_width = (end_x - ini_x) / SCREEN_WIDTH;
    double pixel_height = (end_y - ini_y) / SCREEN_HEIGHT;
    // cout << pixel_width << " " << pixel_height << endl; // This is getting different values

    for(int i = 0; i < SCREEN_WIDTH; i++) {
        // cout << i << endl;
        for(int j = 0; j < SCREEN_HEIGHT; j++) {
            double cur_x = ini_x + i*pixel_width + pixel_width/2;
            double cur_y = ini_y + j*pixel_height + pixel_height/2;
            double cur_z = camera.focal_length;

            // cout << cur_x << " " << cur_y << " " << cur_z << endl;

            // point(eye_ray);
            // set_point(eye_ray, cur_x, cur_y, cur_z);
            glm::dvec3 eye_ray(cur_x, cur_y, cur_z);
            eye_ray = glm::normalize(eye_ray);
            // cout << glm::to_string(eye_ray) << endl;

            // point(eye);
            // set_point(eye, 0,0,0);
            glm::dvec3 eye(0, 0, 0);

            pair<Sphere*, double> intersection = get_closest_intersection(eye_ray, eye, camera.view_matrix, nullptr);
            // cout << intersection.first << endl;
            // if(intersection.first != nullptr) cout << "whhat: " << intersection.first->id << endl;
            glm::dvec3 color = shade(eye_ray, eye, intersection.second, intersection.first, 1);
            // cout << fixed << setprecision(2) << "(" << color.r << ", " << color.g << ", " << color.b << ") " << endl;
            // cout << endl;
            image.push_back(color.r);
            image.push_back(color.g);
            image.push_back(color.b);
        }
        // cout << endl;
    }
    return image;
}


/* ************************************************************************* */
/*                                 OPENGL                                    */
/* ************************************************************************* */

void draw(GLenum primitive, GLuint shader) {

int attrV, attrT; 
    
    glBindBuffer(GL_ARRAY_BUFFER, axisVBO[0]);
    attrV = glGetAttribLocation(shader, "aPosition");
    glVertexAttribPointer(attrV, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(attrV);

    glBindBuffer(GL_ARRAY_BUFFER, axisVBO[2]);
    attrT = glGetAttribLocation(shader, "aTexCoord");
    glVertexAttribPointer(attrT, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(attrT);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, axisVBO[1]);
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    glDrawElements(primitive, 6, GL_UNSIGNED_INT, 0);
    
    glDisableVertexAttribArray(attrV);
    glDisableVertexAttribArray(attrT);

    glBindBuffer(GL_ARRAY_BUFFER, 0); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void display(void) {

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

    vector<float> image = ray_tracing();

    // Enable the texture for OpenGL.
    glEnable(GL_TEXTURE_2D);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); //GL_NEAREST = no smoothing
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, SCREEN_WIDTH, SCREEN_HEIGHT, 0, GL_RGB, GL_FLOAT, &image[0]);
    // glTexImage2D(GL_TEXTURE_2D, 0, 4, u2, v2, 0, GL_RGBA, GL_UNSIGNED_BYTE, &image2[0]);
    
    //PROJECTION
    glm::mat4 orthogonal = glm::ortho(  -10.0, 10.0, -10.0, 10.0);
    // glm::mat4 perspective = glm::perspective(9.0f, 1.0f, 0.1f, 30.0f);
    
    //VIEW
    glm::mat4 View = glm::mat4(1.);
    //View = glm::translate(View, glm::vec3(2.0f,4.0f, -25.0f));

    //MODEL
    glm::mat4 Model1 = glm::rotate(glm::mat4(1.0), 0.f, glm::vec3(1, 0, 0));

    glUseProgram(mainShader);
    int projection_loc = glGetUniformLocation( mainShader, "uProjectionMatrix" );
    int view_loc = glGetUniformLocation( mainShader, "uViewMatrix" );
    int model_loc = glGetUniformLocation( mainShader, "uModelMatrix" );
    glUniformMatrix4fv(projection_loc, 1, GL_FALSE, glm::value_ptr(orthogonal));
    glUniformMatrix4fv(view_loc, 1, GL_FALSE, glm::value_ptr(View));
    glUniformMatrix4fv(model_loc, 1, GL_FALSE, glm::value_ptr(Model1));
    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);


ObjectVA    geometry_VA;

    // ========== FACE ============

    int used_for_faces = 6;

    geometry_VA.vFace = (unsigned int *) malloc(6*sizeof(unsigned int));
    if (!geometry_VA.vFace) {
        cout << "QQQ ISSOOOOO 1" << endl;
        exit(-1);
    }

    geometry_VA.vFace[0] = 0;
    geometry_VA.vFace[1] = 1;
    geometry_VA.vFace[2] = 2;
    
    geometry_VA.vFace[3] = 3;
    geometry_VA.vFace[4] = 2;
    geometry_VA.vFace[5] = 1;

    // ========== POINT ============

    //cout << "used_for_points: " << used_for_points << endl;   

    int used_for_points = 3*4;

    geometry_VA.vPoint  = (float *) malloc(used_for_points*sizeof(float));
    if (!geometry_VA.vPoint) {
        cout << "QQQ ISSOOOOO 2" << endl;
        exit(-1);
    }

    // BOTTOM LEFT 
    geometry_VA.vPoint[0] = -10.0;
    geometry_VA.vPoint[1] = -10.0;
    geometry_VA.vPoint[2] = 0.0;

    // BOTTOM RIGHT 
    geometry_VA.vPoint[3] = 10.0;
    geometry_VA.vPoint[4] = -10.0;
    geometry_VA.vPoint[5] = 0.0;

    // TOP RIGHT 
    geometry_VA.vPoint[6] = 10.0;
    geometry_VA.vPoint[7] = 10.0;
    geometry_VA.vPoint[8] = 0.0;

    // TOP LEFT 
    geometry_VA.vPoint[9] = -10.0;
    geometry_VA.vPoint[10] = 10.0;
    geometry_VA.vPoint[11] = 0.0;


    // ========== TEXTURE COORDINATES ============

    int used_for_texture = used_for_points / 3 * 2;

    geometry_VA.vTextCoord = (float *) malloc(used_for_texture*sizeof(float));
    if (!geometry_VA.vTextCoord) {
        cout << "QQQ ISSOOOOO 3" << endl;
        exit(-1);
    }

    // BOTTOM LEFT 
    geometry_VA.vTextCoord[0] = 0.0;
    geometry_VA.vTextCoord[1] = 1.0;

    // BOTTOM RIGHT 
    geometry_VA.vTextCoord[2] = 1.0;
    geometry_VA.vTextCoord[3] = 1.0;

    // TOP RIGHT 
    geometry_VA.vTextCoord[4] = 1.0;
    geometry_VA.vTextCoord[5] = 0.0;

    // TOP LEFT 
    geometry_VA.vTextCoord[6] = 0.0;
    geometry_VA.vTextCoord[7] = 0.0;

    // cout << "BLAAAAAAAA" << endl;

    glBindBuffer(   GL_ARRAY_BUFFER, axisVBO[0]);

    glBufferData(   GL_ARRAY_BUFFER, used_for_points*sizeof(float), 
                    geometry_VA.vPoint, GL_DYNAMIC_DRAW);

    // glBindBuffer(   GL_ARRAY_BUFFER, axisVBO[1]);

    // glBufferData(   GL_ARRAY_BUFFER, used_for_colors*sizeof(float), 
                    // geometry_VA.vColor, GL_DYNAMIC_DRAW);

    glBindBuffer(   GL_ELEMENT_ARRAY_BUFFER, axisVBO[1]);

    glBufferData(   GL_ELEMENT_ARRAY_BUFFER, used_for_faces*sizeof(unsigned int), 
                    geometry_VA.vFace, GL_DYNAMIC_DRAW);

    // glBindBuffer(   GL_ARRAY_BUFFER, axisVBO[3]);

    // glBufferData(   GL_ARRAY_BUFFER, used_for_points*sizeof(float), 
                    // geometry_VA.vNormal, GL_DYNAMIC_DRAW);

    glBindBuffer(   GL_ARRAY_BUFFER, axisVBO[2]);

    glBufferData(   GL_ARRAY_BUFFER, used_for_texture*sizeof(float), 
                    geometry_VA.vTextCoord, GL_DYNAMIC_DRAW);
    
    
    free(geometry_VA.vPoint);
    // free(geometry_VA.vNormal);
    // free(geometry_VA.vColor);
    free(geometry_VA.vFace);
    free(geometry_VA.vTextCoord);

    draw(GL_TRIANGLE_STRIP, mainShader);

    glutSwapBuffers();

}

void initGL(void) {

    glClearColor(0.0, 0.0, 0.0, 0.0);   
    glPointSize(2.0);
    
    if (glewInit() != GLEW_OK) {
        cout << "Unable to initialize GLEW ... exiting" << endl;
        exit(EXIT_FAILURE);
        }
        
    cout << "Status: Using GLEW " << glewGetString(GLEW_VERSION) << endl;   
    
    cout << "Opengl Version: " << glGetString(GL_VERSION) << endl;
    cout << "Opengl Vendor : " << glGetString(GL_VENDOR) << endl;
    cout << "Opengl Render : " << glGetString(GL_RENDERER) << endl;
    cout << "Opengl Shading Language Version : " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;

}

void initShaders(void) {

    // Load shaders and use the resulting shader program
    mainShader = InitShader( "shader.vert", "shader.frag" );
}

void reshape(int w, int h) {

    SCREEN_WIDTH     = w;
    SCREEN_HEIGHT   = h;
    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);      
    glutPostRedisplay();
}

int main(int argc, char *argv[]) {

    create_scene_objects();

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutCreateWindow("Ray Tracing");
    glutReshapeWindow(SCREEN_WIDTH, SCREEN_HEIGHT);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    // glutKeyboardFunc(keyboard);
    // glutMouseFunc(mouse);

    initGL();

    // LoadTexture("wood-texture.png");

    glGenBuffers(3, axisVBO);
    
    initShaders();
    
    glutMainLoop();


    cout << "HERREEEEE" << endl;

    return(0);

    // matrix eye = camera.position;
    // matrix look_at_vec = normalize(subtractmat(camera.look_at - eye));
}