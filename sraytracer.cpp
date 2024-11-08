#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
using namespace parser;

typedef unsigned char RGB[3];

Vec3f operator+(const Vec3f& a, const Vec3f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Vec3f operator*(const Vec3f& a, const Vec3f& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}
Vec3f operator*(const Vec3f& b, const float a) {
    return {a * b.x, a * b.y, a * b.z};
}
Vec3f operator*(const Vec3i& a, const Vec3f& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}
Vec3f operator*(const float a, const Vec3f& b) {
    return {a * b.x, a * b.y, a * b.z};
}

Vec3f operator-(const Vec3f& a, const Vec3f& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}



float dot(const Vec3f& a, const Vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3f normalize(const Vec3f& v) {
    float len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    if (len == 0) return {0, 0, 0};
    return {v.x / len, v.y / len, v.z / len};
}

Vec3f cross(const Vec3f& a, const Vec3f& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}
float determinant(float a,float b,float c,float d,float e,float f,float g,float h,float i){
    return a*(e*i-h*f) + b*(g*f-d*i) + c*(d*h-e*g);
}


struct HitRecord {
        bool is_intersected;
        float t;
        Vec3f intersection_point;
        Material material;
        Vec3f normal;
};



class Ray {
    
    public: 
    Vec3f origin, direction;
        Ray(const Vec3f& orig, const Vec3f& dir) : origin(orig), direction(normalize(dir)) {}

    Vec3f getDirection() const { return direction; }
    static Ray cameraToPixel(const Camera &camera, int pixel_x, int pixel_y) { //constların yeri hkk bak???
        float l = camera.near_plane.x ;
        float r = camera.near_plane.y;
        float b =camera.near_plane.z;
        float t = camera.near_plane.w;
        Vec3f v = normalize(camera.up);
        Vec3f w = normalize((camera.gaze)*(-1.0));
        Vec3f u = normalize(cross(v,w));
        Vec3f m = camera.position + (camera.gaze)*camera.near_distance ;
        float s_u = (pixel_x + 0.5)* (r-l) / camera.image_width;
        float s_v = (pixel_y + 0.5)* (t-b) / camera.image_height;
        Vec3f q = m + u*l + v*t ;
        Vec3f s = q + u*s_u - v*s_v ;
        return Ray(camera.position, normalize(s - camera.position));
    }


    HitRecord intersectionSphere(Sphere const &sphere, const Scene& scene) {
        HitRecord hit;
        hit.is_intersected = false;
        hit.t=std::numeric_limits<float>::infinity();
        Vec3f center = scene.vertex_data[sphere.center_vertex_id-1];  //the first vertex_id is 1
        Vec3f o_c = this->origin-center;
        Vec3f d = this->direction;
        float A= dot(d,d);
        float B = dot(d,o_c)*2;
        float C = dot(o_c,o_c) - (sphere.radius)*(sphere.radius);
        float discriminant = (B*B - 4*A*C);
        if(discriminant>=0) {
            hit.is_intersected = true;
            float t1 =(-B - sqrt(discriminant)) / (2*A);
            float t2 =(-B + sqrt(discriminant)) / (2*A);
            hit.t = t1 < t2 ? t1 : t2;
            hit.material=scene.materials[sphere.material_id-1]; //burada key gibi veriyoruz sanırım o yüzden +1 yok
            hit.intersection_point = this->origin + d*hit.t ;
            hit.normal = normalize(hit.intersection_point - center);
        }
        return hit;
    }

    HitRecord intersectionTriangle(Triangle const&triangle, const Scene& scene){
        HitRecord hit;
        hit.is_intersected = false;
        hit.t=std::numeric_limits<float>::infinity();
        Vec3f a_b = scene.vertex_data[triangle.indices.v0_id -1] - scene.vertex_data[triangle.indices.v1_id -1] ;
        Vec3f a_c = scene.vertex_data[triangle.indices.v0_id -1] - scene.vertex_data[triangle.indices.v2_id -1] ;
        Vec3f a_o = scene.vertex_data[triangle.indices.v0_id -1] - this->origin ;
        float det_A = determinant(a_b.x,a_b.y,a_b.z,a_c.x,a_c.y,a_c.z,this->direction.x,this->direction.y,this->direction.z);
        float beta = determinant(a_o.x,a_o.y,a_o.z,a_c.x,a_c.y,a_c.z,this->direction.x,this->direction.y,this->direction.z) / det_A ;
        float gama = determinant(a_b.x,a_b.y,a_b.z,a_o.x,a_o.y,a_o.z,this->direction.x,this->direction.y,this->direction.z) / det_A ;
        float t = determinant(a_b.x,a_b.y,a_b.z,a_c.x,a_c.y,a_c.z,a_o.x,a_o.y,a_o.z) / det_A ;
        if(fabs(det_A) < 1e-6) return hit;
        if (beta+gama <= 1 && 0 <= beta && 0 <= gama && t>0){ //tmin<=t<=tmax ???
            hit.is_intersected = true;
            hit.material = scene.materials[triangle.material_id-1];
            hit.t = t;
            hit.intersection_point = this->origin + (this->direction)*(hit.t) ;
            hit.normal = normalize(cross(scene.vertex_data[triangle.indices.v2_id -1] - scene.vertex_data[triangle.indices.v1_id -1] ,a_b));
        }
        return hit;
    }

    HitRecord intersectionMesh(Mesh const&mesh, const Scene& scene){
        HitRecord hit;
        HitRecord current_hit;
        hit.is_intersected = false;
        hit.t=std::numeric_limits<float>::infinity();
        Triangle mesh_triangle;

        for(int i = 0; i<mesh.faces.size(); i++){
            mesh_triangle = {mesh.material_id, mesh.faces[i]};
            current_hit = intersectionTriangle(mesh_triangle,scene);
            if(current_hit.is_intersected==true) {
                if(current_hit.t < hit.t) {hit = current_hit; }
            } 
        }
        return hit;
    }

    HitRecord find_closest_hit(Scene& scene){
        HitRecord hit;
        HitRecord current_hit;
        hit.is_intersected = false;
        hit.t=std::numeric_limits<float>::infinity();
        for(int i=0; i< scene.spheres.size() ; i++){
            current_hit =intersectionSphere(scene.spheres[i],scene);
            if(current_hit.is_intersected==true) {
                if(current_hit.t < hit.t) {hit = current_hit; }
            } 
        }

        for(int i=0; i< scene.triangles.size() ; i++){
            current_hit =intersectionTriangle(scene.triangles[i],scene);
            if(current_hit.is_intersected==true) {
                if(current_hit.t < hit.t) {hit = current_hit; }
            } 
        }

       for(int i=0; i< scene.meshes.size() ; i++){
            current_hit =intersectionMesh(scene.meshes[i],scene);
            if(current_hit.is_intersected==true) {
                if(current_hit.t < hit.t) {hit = current_hit; }
            } 
        }

        return hit;
    }

 
}

;

Vec3f compute_color( Ray& ray, int depth,  Scene& scene);
Vec3f apply_shading( Ray& ray, int depth,  HitRecord& hit,  Scene& scene){

    Vec3f color = {0.0f, 0.0f, 0.0f};
    color = color + (scene.ambient_light * hit.material.ambient);
    Vec3f eyeVector = normalize(ray.origin - hit.intersection_point);

    if (hit.material.is_mirror && depth < scene.max_recursion_depth){
        float cosTheta = dot(hit.normal, eyeVector);
        Vec3f refRayDirection = eyeVector * (-1.0f) + hit.normal * (cosTheta * (2.0f));
        Vec3f refRayOrigin = hit.intersection_point + (refRayDirection * (scene.shadow_ray_epsilon));
        Ray reflectedRay(refRayDirection, refRayOrigin);
        HitRecord reflectHit = reflectedRay.find_closest_hit(scene);
        if (reflectHit.is_intersected && reflectHit.t>0.0f){
            Vec3f refRayColor = compute_color(reflectedRay, depth+1, scene);
            color = color + (refRayColor * hit.material.mirror);
        }
    }

    for ( PointLight& pointLight : scene.point_lights){

        Vec3f objToLight = pointLight.position - hit.intersection_point;
        if (dot(objToLight,hit.normal) < 0) continue;
        

            Vec3f normalOToL = normalize(objToLight);
            Vec3f shadowRayOrigin = hit.intersection_point + normalOToL * (scene.shadow_ray_epsilon);
            
            Ray shadowRay = Ray(shadowRayOrigin,normalOToL);
            HitRecord shadowHit = shadowRay.find_closest_hit(scene);
            
            float distToLightSource = sqrt(objToLight.x * objToLight.x + objToLight.y * objToLight.y + objToLight.z * objToLight.z);
            float minT =shadowHit.t;
            Vec3f vec = shadowHit.intersection_point - hit.intersection_point;
            float distance = sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
            if ((minT > 0 && distance > distToLightSource) || minT == std::numeric_limits<float>::infinity()) {
                
                Vec3f irradiance = {0.0f, 0.0f, 0.0f};
                if (distToLightSource>0){
                    irradiance.x = pointLight.intensity.x / pow(distToLightSource,2);
                    irradiance.y = pointLight.intensity.y / pow(distToLightSource,2);
                    irradiance.z = pointLight.intensity.z / pow(distToLightSource,2);
                }


                Vec3f diffuseComponent = {0.0f, 0.0f, 0.0f};
                float clampedDotProduct = dot(hit.normal, normalOToL);
                diffuseComponent.x = hit.material.diffuse.x * irradiance.x * clampedDotProduct;
                diffuseComponent.y = hit.material.diffuse.y * irradiance.y * clampedDotProduct;
                diffuseComponent.z = hit.material.diffuse.z * irradiance.z * clampedDotProduct;
                printf("Burası diffuse:%f %f %f %f\n", clampedDotProduct, diffuseComponent.x, diffuseComponent.y, diffuseComponent.z);

                Vec3f normalHalfwayVector = normalize(normalOToL + eyeVector);
                Vec3f specularComponent = {0.0f, 0.0f, 0.0f};
                float cosAlpha = std::max(0.0f, dot(hit.normal, normalHalfwayVector));
                float b = pow(cosAlpha, hit.material.phong_exponent);
                specularComponent.x = hit.material.specular.x * irradiance.x * b;
                specularComponent.y = hit.material.specular.y * irradiance.y * b;
                specularComponent.z = hit.material.specular.z * irradiance.z * b;

                color = (color + specularComponent) + diffuseComponent;
                


            }

    }
    return color;

}


Vec3f compute_color( Ray& ray, int depth,  Scene& scene){

    HitRecord hit = ray.find_closest_hit(scene);
    Vec3f zero = {0.0f, 0.0f, 0.0f};

    if (depth>scene.max_recursion_depth){
        return zero;
    }

    else if (hit.is_intersected==false) {
        zero.x = scene.background_color.x;
        zero.y = scene.background_color.y;
        zero.z = scene.background_color.z;
        return zero;
    }

    else if (hit.is_intersected){
        
        return apply_shading(ray, depth, hit, scene);
    }
    
    else {
        return zero;
    }

}





int main()
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;
    scene.loadFromXml("inputs/mirror_spheres.xml");


    int width = scene.cameras[0].image_width;
    int height = scene.cameras[0].image_height;

    unsigned char* image = new unsigned char [width * height * 3];

    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {   
            Ray ray = Ray::cameraToPixel(scene.cameras[0],x,y);
            Vec3f color = compute_color(ray,0,scene);
            
            int roundedX = static_cast<int>(std::round(color.x));
            int roundedY = static_cast<int>(std::round(color.y));
            int roundedZ = static_cast<int>(std::round(color.z));

            image[i++] = roundedX > 255 ? 255 : roundedX;
            image[i++] = roundedY > 255 ? 255 : roundedY;
            image[i++] = roundedZ > 255 ? 255 : roundedZ;
        }
    }
    
    write_ppm("monkey.ppm", image, width, height);

}
