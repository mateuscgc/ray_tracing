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

#define PI acos(-1)

int CAMERA_VIEW_WIDTH = 300;
int CAMERA_VIEW_HEIGHT = 300;

int SCREEN_WIDTH;
int SCREEN_HEIGHT;

int MAX_RECURSION = 3;

GLuint  mainShader;
GLuint  axisVBO[5];

double lastTime;

double last_x = -1;
double last_y = -1;

struct Implicit {
    int id;
    double transparency_rate;
    double refraction_distortion_rate;
    virtual ~Implicit() {}

    virtual bool intersect(const glm::vec3& ray, const glm::vec3& ori) = 0;
    virtual vector<double> intersections_time(const glm::vec3& ray, const glm::vec3& ori, const glm::vec3& camera_position) = 0;
    // shade
    virtual glm::vec3 shadeAmb(const glm::vec3& position) = 0;
    virtual glm::vec3 shadeDif(const glm::vec3& position) = 0;
    virtual glm::vec3 shadeSpec(const glm::vec3& position) = 0;
    virtual glm::vec3 getNormal(const glm::vec3& ori, const glm::vec3& position) = 0;
    virtual void toCameraSystem(const glm::mat4& view_matrix) = 0;
protected:
    Implicit(
        int id,
        double transparency_rate,
        double refraction_distortion_rate
    ) : 
        id(id),
        transparency_rate(transparency_rate),
        refraction_distortion_rate(refraction_distortion_rate)
    {}
};

struct Sphere : public Implicit {
    glm::vec3 position;
    glm::vec3 camera_system_position;
    double radius;
    // int id;
    // double transparency_rate;
    // double refraction_distortion_rate;

    glm::vec3 matAmb;
    glm::vec3 matDif;
    glm::vec3 matSpec;

    // Color color;

    // Sphere() {}
    Sphere(
            glm::vec3 position,
            double radius,
            int id,
            double transparency_rate,
            double refraction_distortion_rate,
            glm::vec3 matAmb,
            glm::vec3 matDif,
            glm::vec3 matSpec
        ) :
            Implicit(id, transparency_rate, refraction_distortion_rate),
            position(position),
            radius(radius),
            matAmb(matAmb),
            matDif(matDif),
            matSpec(matSpec)
    {}

    static double a(const glm::vec3& ray) {
        return ray.x*ray.x + ray.y*ray.y + ray.z*ray.z;
    }

    static double b(const glm::vec3& ray, const glm::vec3& ori) {
        return 2*(ray.x*ori.x + ray.y*ori.y + ray.z*ori.z);
    }

    static double c(const glm::vec3& ori, double radius) {
        return ori.x*ori.x + ori.y*ori.y + ori.z*ori.z - radius*radius;
    }

    static double discriminant(const glm::vec3& ray, const glm::vec3& ori, double radius) {
        double a = Sphere::a(ray);
        double b = Sphere::b(ray, ori);
        double c = Sphere::c(ori, radius);
        return b*b - 4*a*c;
    }
    virtual bool intersect(const glm::vec3& ray, const glm::vec3& ori) {
        glm::vec3 new_ori = ori - camera_system_position;
        return Sphere::discriminant(ray, new_ori, radius) >= 0;
    }
    virtual vector<double> intersections_time(const glm::vec3& ray, const glm::vec3& ori, const glm::vec3& camera_position) {
        glm::vec3 new_ori = ori - camera_system_position;
        double a = Sphere::a(ray);
        double b = Sphere::b(ray, new_ori);
        double root1 = (-b + sqrt(Sphere::discriminant(ray, new_ori, radius)))/(2*a);
        double root2 = (-b - sqrt(Sphere::discriminant(ray, new_ori, radius)))/(2*a);
        return { root1, root2 };
    }

    virtual glm::vec3 shadeAmb(const glm::vec3& position) {
        return matAmb;
    }
    virtual glm::vec3 shadeDif(const glm::vec3& position) {
        return matDif;
    }
    virtual glm::vec3 shadeSpec(const glm::vec3& position) {
        return matSpec;
    }
    virtual glm::vec3 getNormal(const glm::vec3& ori, const glm::vec3& position) {
        return glm::normalize(position - camera_system_position);
    }

    virtual void toCameraSystem(const glm::mat4& view_matrix) {
        camera_system_position = glm::vec3(view_matrix * glm::dvec4(this->position, 1.0));
    }

};

struct Plane : Implicit {
    glm::vec3 normal;
    glm::vec3 camera_system_normal;
    glm::vec3 position;
    glm::vec3 camera_system_position;
    double d;
    double d_camera_system;

    glm::vec3 matAmb[2];
    glm::vec3 matDif[2];
    glm::vec3 matSpec[2];

    double square_side;

    Plane(
            glm::vec3 normal,
            double d,
            int id,
            double transparency_rate,
            double refraction_distortion_rate,
            glm::vec3 matAmb1,
            glm::vec3 matAmb2,
            glm::vec3 matDif1,
            glm::vec3 matDif2,
            glm::vec3 matSpec1,
            glm::vec3 matSpec2,
            double square_side
        ) :
            Implicit(id, transparency_rate, refraction_distortion_rate),
            normal(normal),
            d(d),
            matAmb{matAmb1, matAmb2},
            matDif{matDif1, matDif2},
            matSpec{matSpec1, matSpec2},
            square_side(square_side)
    {
        position = normal * d;
    }

    virtual bool intersect(const glm::vec3& ray, const glm::vec3& ori) {
        // cout << glm::dot(normal, ray) << endl;
        return glm::dot(camera_system_normal, ray) != 0;
    }

    virtual vector<double> intersections_time(const glm::vec3& ray, const glm::vec3& ori, const glm::vec3& camera_position) {
        // glm::vec3 new_ori = ori - camera.position;
        return { -((glm::dot(camera_system_normal, ori) - d_camera_system) / glm::dot(camera_system_normal, ray)) };
        // return { -(normal.x*ori.x + normal.y*ori.y + normal.z*ori.z - d) / (normal.x*ray.x + normal.y*ray.y + normal.z*ray.z) };
    }

    int get_color_index(const glm::vec3& position) {
        int x = (int)(glm::distance(position, glm::vec3(position.x, camera_system_position.y, camera_system_position.z)) / square_side) % 2;
        int z = (int)(glm::distance(position, glm::vec3(camera_system_position.x, camera_system_position.y, position.z)) / square_side) % 2;
        return x ^ z;
    }

    virtual glm::vec3 shadeAmb(const glm::vec3& position) {
        return matAmb[get_color_index(position)];
    }
    virtual glm::vec3 shadeDif(const glm::vec3& position) {
        return matDif[get_color_index(position)];
    }
    virtual glm::vec3 shadeSpec(const glm::vec3& position) {
        return matSpec[get_color_index(position)];
    }
    virtual glm::vec3 getNormal(const glm::vec3& ori, const glm::vec3& position) {
        return camera_system_normal;
    }
    virtual void toCameraSystem(const glm::mat4& view_matrix) {
        camera_system_normal = glm::mat3(view_matrix) * normal;
        camera_system_position = glm::vec3(view_matrix * glm::vec4(position, 1.0));
        d_camera_system = glm::dot(camera_system_position, camera_system_normal);
    }
};

vector< Implicit* >  objects;

struct Light {
    glm::vec3 position;
    glm::vec3 camera_system_position;
    glm::vec3 color;
    int id;

    glm::vec3 atenuation(glm::vec3 target, glm::mat4 view) {
        return this->color / max(1.0f, glm::distance(this->camera_system_position, target) / 4);
    }
};

vector<Light> lights;

struct Camera {
    glm::vec3 position;
    glm::vec3 target;
    glm::vec3 up_vec;
    glm::mat4 view_matrix;
    double field_of_view;

    void updateView() {
        this->view_matrix = glm::lookAt(this->position, this->target, this->up_vec);
    }

    void w() {
        glm::vec3 direction = -glm::normalize(this->target - this->position);
        this->position = glm::vec3(glm::translate(direction/5) * glm::dvec4(position, 1));
        this->target = glm::vec3(glm::translate(direction/5) * glm::dvec4(target, 1));
    }

    void s() {
        glm::vec3 direction = glm::normalize(this->target - this->position);
        this->position = glm::vec3(glm::translate(direction/5) * glm::dvec4(position, 1));
        this->target = glm::vec3(glm::translate(direction/5) * glm::dvec4(target, 1));
    }

    void a() {
        glm::vec3 direction = -glm::cross(normalize(this->target - this->position), this->up_vec);
        this->position = glm::vec3(glm::translate(direction/5) * glm::dvec4(position, 1));
        this->target = glm::vec3(glm::translate(direction/5) * glm::dvec4(target, 1));
    }

    void d() {
        glm::vec3 direction = glm::cross(normalize(this->target - this->position), this->up_vec);
        this->position = glm::vec3(glm::translate(direction/5) * glm::dvec4(position, 1));
        this->target = glm::vec3(glm::translate(direction/5) * glm::dvec4(target, 1));
    } 
};

struct FPSCamera {

    glm::vec3 position;
    float pitch;
    float yaw;
    glm::mat4 view_matrix;
    double field_of_view;

    // Pitch should be in the range of [-90 ... 90] degrees and yaw
    // should be in the range of [0 ... 360] degrees.
    void updateView() {
        // If the pitch and yaw angles are in degrees,
        // they need to be converted to radians. Here
        // I assume the values are already converted to radians.
        float cosPitch = cos(this->pitch);
        float sinPitch = sin(this->pitch);
        float cosYaw = cos(this->yaw);
        float sinYaw = sin(this->yaw);
     
        glm::vec3 xaxis = { cosYaw, 0, -sinYaw };
        glm::vec3 yaxis = { sinYaw * sinPitch, cosPitch, cosYaw * sinPitch };
        glm::vec3 zaxis = { sinYaw * cosPitch, -sinPitch, cosPitch * cosYaw };
     
        // Create a 4x4 view matrix from the right, up, forward and eye position vectors
        this->view_matrix = {
            glm::vec4(       xaxis.x,            yaxis.x,            zaxis.x,      0 ),
            glm::vec4(       xaxis.y,            yaxis.y,            zaxis.y,      0 ),
            glm::vec4(       xaxis.z,            yaxis.z,            zaxis.z,      0 ),
            glm::vec4( -dot( xaxis, this->position ), -dot( yaxis, this->position ), -dot( zaxis, this->position ), 1 )
        };
    }

    void w() {
        this->position.y += 0.2;
    }

    void s() {
        this->position.y += -0.2;
    }

    void a() {
        this->position.x += -0.2;
    }

    void d() {
        this->position.x += 0.2;
    }

    void q() {
        this->position.z += 0.2;
    }

    void e() {
        this->position.z += -0.2;
    }

    void i() {
        this->pitch = glm::clamp(pitch - PI/45, -PI/2, PI/2);
    }

    void k() {
        this->pitch = glm::clamp(pitch + PI/45, -PI/2, PI/2);
    }
    
    void j() {
        this->yaw = yaw - PI/45;
        this->yaw = fmod(yaw, PI*2); 
    }

    void l() {
        this->yaw = yaw + PI/45;
        this->yaw = fmod(yaw, PI*2);
    }
};

FPSCamera camera;

void create_scene_objects() {
    int id = -1;

    // ===================== SPHERES =======================

    // Transparent
    objects.push_back(new Sphere(
        glm::vec3(0,1,5.0),
        0.5,
        ++id,
        0.9,
        0.92,
        glm::vec3(0.7,0.4,0.2),
        glm::vec3(0.7,0.4,0.2),
        glm::vec3(0.8,0.8,0.8)
    ));


    // Blue
    objects.push_back(new Sphere(
        glm::vec3(0.6,1.2,6.0),
        0.4,
        ++id,
        0.1,
        0.95,
        glm::vec3(0.2,0.2,0.6),
        glm::vec3(0.2,0.2,0.6),
        glm::vec3(0.2,0.2,0.2)
    ));

    // White
    objects.push_back(new Sphere(
        glm::vec3(-0.4,1.2,6.0),
        0.3,
        ++id,
        0.0,
        0.75,
        glm::vec3(0.8,0.8,0.8),
        glm::vec3(0.8,0.8,0.8),
        glm::vec3(0.6,0.6,0.6)
    ));

    // Red
    objects.push_back(new Sphere(
        glm::vec3(0.2,1.5,8.0),
        0.5,
        ++id,
        0.3,
        0.9,
        glm::vec3(0.8,0.3,0.2),
        glm::vec3(0.8,0.3,0.2),
        glm::vec3(0.4,0.4,0.4)
    ));

    // Green
    objects.push_back(new Sphere(
        glm::vec3(-0.3,1,2.5),
        0.3,
        ++id,
        0.0,
        0.9,
        glm::vec3(0.4,0.7,0.3),
        glm::vec3(0.4,0.7,0.3),
        glm::vec3(0.4,0.4,0.4)
    ));

    // ===================== PLANES =======================

    objects.push_back(new Plane(
        glm::vec3(0.0, 1.0, 0.0),
        0,
        ++id,
        0,
        1,
        glm::vec3(0.7,0.0,0.0),
        glm::vec3(0.7,0.7,0.0),

        glm::vec3(0.7,0.0,0.0),
        glm::vec3(0.7,0.7,0.0),

        glm::vec3(0.6,0.6,0.6),
        glm::vec3(0.2,0.2,0.2),

        0.1
    ));

    // ======================= LIGHTS =====================

    Light l;

    l.position = glm::vec3(-0.5,2.5,2);
    l.color = glm::vec3(1,1,1);
    l.id = ++id;
    lights.push_back(l);


    // l.position = glm::vec3(-3,2.0,4);
    // l.color = glm::vec3(1,1,1);
    // // l.color = Color(1,1,1);
    // l.id = ++id;   
    // lights.push_back(l);
 
    // ======================== CAMERA ====================

    // camera.position = glm::vec3(0,0,0);
    // camera.target = glm::vec3(0,0,-1);

    // // camera.position = glm::vec3(0,1,2);
    // // camera.target = glm::vec3(0,2,-1);

    // // camera.position = glm::vec3(-1,1.8,2);
    // // camera.target = glm::vec3(-2,3,-2);
    
    // camera.up_vec = glm::vec3(0,1,0);
    // camera.updateView();

    // camera.field_of_view = PI/4; // 45º

    // ======================== FPS CAMERA ====================

    camera.position = glm::vec3(0,1,-2);
    camera.pitch = 0.0;
    camera.yaw = 0.0;
    camera.updateView();

    camera.field_of_view = PI/4; // 45º
}


pair< Implicit*, double> get_closest_intersection(const glm::vec3& ray, const glm::vec3& ori, Implicit *avoid) {
    Implicit *ans = nullptr;
    double best_time = 0;
    for(Implicit *s : objects) {
        if((avoid == nullptr || s->id != avoid->id) && s->intersect(ray, ori)) {
            vector<double> intersections;
            intersections = s->intersections_time(ray, ori, camera.position);

            double best = -1;
            for(double t : intersections) {
                if(t > 0.0001 && (best == -1 || t < best))
                    best = t;
            }
            if(best <= 0) continue;

            if(ans == nullptr || best < best_time) {
                best_time = best;
                ans = s;
            }
        }
    }
    return make_pair(ans, best_time);
}

glm::vec3 shade(const glm::vec3& ray, const glm::vec3& ray_ori, double intersection_time, Implicit *obj, int recursion_level) {
    if(obj == nullptr) {
        // Background color
        // return glm::vec3(0.1,0.75,0.9);
        return glm::vec3(0,0,0);
    }

    glm::vec3 result(0,0,0);

    glm::vec3 intersection_point = ray_ori + glm::normalize(ray) * intersection_time;
    glm::vec3 intersection_normal = obj->getNormal(ray_ori, intersection_point);

    // cout << "ray_ori " << glm::to_string(ray_ori) << endl;
    // cout << "intersection_point " << glm::to_string(intersection_point) << endl;
    glm::vec3 intersection_to_eye_vec = -ray;

    pair<Implicit*, double> intersection;

    glm::vec3 matAmb = obj->shadeAmb(intersection_point);
    glm::vec3 matDif = obj->shadeDif(intersection_point);
    glm::vec3 matSpec = obj->shadeSpec(intersection_point);

    // if(obj->id == 5) {
    //     cout << "intersection_point " << glm::to_string(intersection_point) << endl;
    //     cout << "intersection_normal " << glm::to_string(intersection_normal) << endl;
    //     cout << "intersection_to_eye_vec " << glm::to_string(intersection_to_eye_vec) << endl;
    //     cout << "matAmb " << glm::to_string(matAmb) << endl;
    //     cout << "matDif " << glm::to_string(matDif) << endl;
    //     cout << "matSpec " << glm::to_string(matSpec) << endl;
    // }


    for(Light l : lights) {
        glm::vec3 light_ray = glm::normalize(intersection_point - l.camera_system_position);
        // cout << "light ray " << glm::to_string(-light_ray) << endl;
        
        // Ambiente
        // Put ambient even with oclusion
        // Better ambient
        glm::vec3 ambient = l.color * matAmb;
        // cout << "ambient " << glm::to_string(ambient) << endl;

        // Check if light oclusion
        // Need to check oclusion from obj itself
        intersection = get_closest_intersection(light_ray, l.camera_system_position, nullptr);
        if(intersection.first != nullptr && intersection.first->id != obj->id) {
            result += (1-obj->transparency_rate) * ambient;
            continue;
        }

        // Diffuse
        double teta = max(glm::dot(-light_ray, intersection_normal), 0.0f);
        glm::vec3 diffuse = l.atenuation(intersection_point, camera.view_matrix) * matDif * teta;

        // Specular
        glm::vec3 light_reflected_ray = glm::normalize(glm::reflect(light_ray, intersection_normal));
        double omega = max(glm::dot(intersection_to_eye_vec, light_reflected_ray), 0.0f);
        glm::vec3 specular = l.atenuation(intersection_point, camera.view_matrix) * matSpec * pow(omega, 2);
        // cout << "omega " << omega << " specular " << glm::to_string(specular) << endl;

        glm::vec3 light_result = (1.0-obj->transparency_rate) * (ambient + diffuse + specular);
        glm::clamp(light_result, 0.0f, 1.0f);
        // cout << "light_result " << glm::to_string(light_result) << endl;
        result += light_result;
        // result += ambient;
    }

    // glm::clamp(result, 0.0f, 1.0f);

    if(recursion_level < MAX_RECURSION) {
        if(obj->transparency_rate < 1) {
            glm::vec3 reflected_ray = glm::normalize(glm::reflect(ray, intersection_normal));
            intersection = get_closest_intersection(reflected_ray, intersection_point, nullptr);
            result += (1.0-obj->transparency_rate) * matSpec * shade(reflected_ray, intersection_point, intersection.second, intersection.first, recursion_level+1);
        }

        if(obj->transparency_rate > 0) {
            glm::vec3 transmited_ray = glm::refract(ray, intersection_normal, (float)obj->refraction_distortion_rate); // Do some distortion
            // Maybe add obj diameter to ray origin or something more precise to avoid get the intersection with the same object (from inside)
            // Or get intersection thats is not with the same object
            intersection = get_closest_intersection(transmited_ray, intersection_point, obj);
            result += obj->transparency_rate * shade(transmited_ray, intersection_point, intersection.second, intersection.first, recursion_level+1);
        }
    }

    // glm::clamp(result, 0.0f, 1.0f);

    return result;
}

vector<float> ray_tracing() {
    vector<float> image;

    for(Implicit *s : objects) {
        s->toCameraSystem(camera.view_matrix);
    }

    for(Light &l : lights) {
        l.camera_system_position = glm::vec3(camera.view_matrix * glm::dvec4(l.position, 1.0));
    }

    double pixel_size = tan(camera.field_of_view/2) / CAMERA_VIEW_HEIGHT;
    double CAMERA_VIEW_HALF_WIDTH = CAMERA_VIEW_WIDTH/2;
    double CAMERA_VIEW_HALF_HEIGHT = CAMERA_VIEW_HEIGHT/2;
    double HALF_PIXEL_SIZE = pixel_size/2;

    double cur_z = 1;
    
    for(int j = 0; j < CAMERA_VIEW_HEIGHT; j++) {
        for(int i = 0; i < CAMERA_VIEW_WIDTH; i++) {
            double cur_x = (i - CAMERA_VIEW_HALF_WIDTH)*pixel_size - HALF_PIXEL_SIZE;
            double cur_y = (j - CAMERA_VIEW_HALF_HEIGHT)*pixel_size - HALF_PIXEL_SIZE;

            glm::vec3 eye_ray(cur_x, cur_y, cur_z);
            eye_ray = glm::normalize(eye_ray);

            glm::vec3 eye(0, 0, 0); // Camera system

            pair<Implicit*, double> intersection = get_closest_intersection(eye_ray, eye, nullptr);
            glm::vec3 color = shade(eye_ray, eye, intersection.second, intersection.first, 1);
            image.push_back(color.r);
            image.push_back(color.g);
            image.push_back(color.b);
        }
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

    // cout << "Ray Tracing cost: " << (glutGet(GLUT_ELAPSED_TIME)-lastTime)/1000 << endl;

    // Enable the texture for OpenGL.
    glEnable(GL_TEXTURE_2D);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); //GL_NEAREST = no smoothing
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, CAMERA_VIEW_WIDTH, CAMERA_VIEW_HEIGHT, 0, GL_RGB, GL_FLOAT, &image[0]);
    // glTexImage2D(GL_TEXTURE_2D, 0, 4, u2, v2, 0, GL_RGBA, GL_UNSIGNED_BYTE, &image2[0]);
    
    //PROJECTION
    glm::mat4 orthogonal = glm::ortho(  -10.0, 10.0, -10.0, 10.0);
    
    //VIEW
    glm::mat4 View = glm::mat4(1.);

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
        exit(-1);
    }

    geometry_VA.vFace[0] = 0;
    geometry_VA.vFace[1] = 1;
    geometry_VA.vFace[2] = 3;
    
    geometry_VA.vFace[3] = 1;
    geometry_VA.vFace[4] = 2;
    geometry_VA.vFace[5] = 3;

    // ========== POINT ============

    int used_for_points = 3*4;

    geometry_VA.vPoint  = (float *) malloc(used_for_points*sizeof(float));
    if (!geometry_VA.vPoint) {
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
        exit(-1);
    }

    // BOTTOM LEFT 
    geometry_VA.vTextCoord[0] = 0.0;
    geometry_VA.vTextCoord[1] = 0.0;

    // BOTTOM RIGHT 
    geometry_VA.vTextCoord[2] = 1.0;
    geometry_VA.vTextCoord[3] = 0.0;

    // TOP RIGHT 
    geometry_VA.vTextCoord[4] = 1.0;
    geometry_VA.vTextCoord[5] = 1.0;

    // TOP LEFT 
    geometry_VA.vTextCoord[6] = 0.0;
    geometry_VA.vTextCoord[7] = 1.0;

    glBindBuffer(   GL_ARRAY_BUFFER, axisVBO[0]);

    glBufferData(   GL_ARRAY_BUFFER, used_for_points*sizeof(float), 
                    geometry_VA.vPoint, GL_DYNAMIC_DRAW);

    glBindBuffer(   GL_ELEMENT_ARRAY_BUFFER, axisVBO[1]);

    glBufferData(   GL_ELEMENT_ARRAY_BUFFER, used_for_faces*sizeof(unsigned int), 
                    geometry_VA.vFace, GL_DYNAMIC_DRAW);

    glBindBuffer(   GL_ARRAY_BUFFER, axisVBO[2]);

    glBufferData(   GL_ARRAY_BUFFER, used_for_texture*sizeof(float), 
                    geometry_VA.vTextCoord, GL_DYNAMIC_DRAW);
    
    
    free(geometry_VA.vPoint);
    free(geometry_VA.vFace);
    free(geometry_VA.vTextCoord);

    draw(GL_TRIANGLE_STRIP, mainShader);

    glutSwapBuffers();

    // Measure speed
    double currentTime = glutGet(GLUT_ELAPSED_TIME);
    double frameCost = (currentTime-lastTime)/1000;
    cout << "Frame Cost: " << frameCost << endl;
    cout << "FPS: " << 1/frameCost << endl;
    lastTime = currentTime;

    glutPostRedisplay();
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

    SCREEN_WIDTH    = w;
    SCREEN_HEIGHT   = h;
    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
}

/* ************************************************************************* */
/*                                Controls                                   */
/* ************************************************************************* */

void keyboard (unsigned char key, int x, int y) {

    switch (key) {
        case 27     :   exit(0);
                        break;
        case 'W'    :
        case 'w'    :
                        camera.w();
                        break;
        case 'S'    :
        case 's'    :
                        camera.s();
                        break;

        case 'A'    :
        case 'a'    :
                        camera.a();
                        break;
        case 'D'    :
        case 'd'    :
                        camera.d();
                        break;

        case 'Q'    :
        case 'q'    :
                        camera.q();
                        break;

        case 'E'    :
        case 'e'    :
                        camera.e();
                        break;

        case 'I'    :
        case 'i'    :
                        camera.i();
                        break;
        case 'K'    :
        case 'k'    :
                        camera.k();
                        break;

        case 'J'    :
        case 'j'    :
                        camera.j();
                        break;
        case 'L'    :
        case 'l'    :
                        camera.l();
                        break;

    }
    camera.updateView();
    // camera.view_matrix = glm::lookAt(camera.position, camera.target, camera.up_vec);
}

void mouseMotion(int x, int y) {
    if(last_y != -1) camera.pitch += (y - last_y)*PI/360;
    if(last_x != -1) camera.yaw += (x - last_x)*PI/360;

    last_x = x;
    last_y = y;

    camera.updateView();
}

/* ************************************************************************* */
/*                                  MAIN                                     */
/* ************************************************************************* */

int main(int argc, char *argv[]) {
    ios_base::sync_with_stdio(false); cin.tie(0);

    create_scene_objects();

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutCreateWindow("Ray Tracing");
    glutReshapeWindow(CAMERA_VIEW_WIDTH, CAMERA_VIEW_HEIGHT);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    // glutPassiveMotionFunc(mouseMotion);
    // glutMotionFunc(mouseMotion);
    
    // glutMouseFunc(mouse);

    initGL();

    // LoadTexture("wood-texture.png");

    glGenBuffers(3, axisVBO);
    
    initShaders();


    lastTime = glutGet(GLUT_ELAPSED_TIME);
    glutMainLoop();

    return(0);
}