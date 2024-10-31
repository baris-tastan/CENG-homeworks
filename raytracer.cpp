#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
using namespace parser;

typedef unsigned char RGB[3];

Vec3f operator+(const Vec3f& a, const Vec3f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3f operator-(const Vec3f& a, const Vec3f& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3f operator*(const Vec3f& v, float s) {
    return {v.x * s, v.y * s, v.z * s};
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
    Vec3f origin, direction;
    public: 
        Ray(const Vec3f& orig, const Vec3f& dir) : origin(orig), direction(normalize(dir)) {}

    //we can use this as: Ray ray = Ray::cameraToPixel(camera, pixel_x, pixel_y);
    //The cameraToPixel method is static because it creates a Ray without needing an instance of the Ray class. 
    //It uses the camera's parameters and pixel coordinates to calculate the origin and direction of the ray.
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
        
        hit.normal = normalize(cross(scene.vertex_data[triangle.indices.v2_id -1] - scene.vertex_data[triangle.indices.v1_id -1] ,a_b));

        if (dot(hit.normal,this.direction) > 0) return hit;
        else (beta+gama <= 1 && 0 <= beta && 0 <= gama && t>0){ 

            hit.is_intersected = true;
            hit.material = scene.materials[triangle.material_id-1];
            hit.t = t;
            hit.intersection_point = this->origin + (this->direction)*(hit.t) ;
            
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

 
};

float calculateDistance(Vec3f a, Vec3f b){
    return sqrt(pow(a.x - b.x ,2)+pow(a.y - b.y,2)+pow(a.z - b.z,2));
}

Vec3f applyShading(Ray ray, int depth, HitRecord hit, Scene& scene);
Vec3f compute_color(Ray ray,int depth, Scene& scene);

bool inShadow(Vec3f intersection_point, PointLight I, Scene& scene){ //kendisiyle çakışması ve arasında olmasının garantisi
        //printf("burası");
        Vec3f shadowRayDir = normalize(I.position - intersection_point);
        Ray shadowRay = Ray((intersection_point+ shadowRayDir * (scene.shadow_ray_epsilon)),shadowRayDir); //offseti sildim
        float objToLightDistance, objToObjDistance;
        objToLightDistance = calculateDistance(intersection_point,I.position);
        HitRecord hit;
        hit.is_intersected = false;

        for(int i=0; i<scene.spheres.size(); i++){

            hit=shadowRay.intersectionSphere(scene.spheres[i],scene);

            if(hit.is_intersected && calculateDistance(hit.intersection_point,intersection_point)>scene.shadow_ray_epsilon
            && calculateDistance(hit.intersection_point,I.position) < objToLightDistance
            &&calculateDistance(hit.intersection_point,intersection_point) < objToLightDistance ) {
            return true;}
            
        }

        for(int i=0; i< scene.triangles.size(); i++){

            hit=shadowRay.intersectionTriangle(scene.triangles[i],scene);

            if(hit.is_intersected && calculateDistance(hit.intersection_point,intersection_point)>scene.shadow_ray_epsilon
            &&calculateDistance(hit.intersection_point,I.position) < objToLightDistance
            &&calculateDistance(hit.intersection_point,intersection_point) < objToLightDistance ) {
            return true;}
            
        }

        for(int i=0; i<scene.meshes.size(); i++){

            hit =shadowRay.intersectionMesh(scene.meshes[i],scene);

            if(hit.is_intersected && calculateDistance(hit.intersection_point,intersection_point)>scene.shadow_ray_epsilon
            && calculateDistance(hit.intersection_point,I.position) < objToLightDistance
            &&calculateDistance(hit.intersection_point,intersection_point) < objToLightDistance ) {
            return true;}
        }

        return false;
}

Vec3f applyShading(Ray ray, int depth, HitRecord hit, Scene& scene){
    Vec3f color= {0,0,0};
    //ambient contribution
    color.x += (scene.ambient_light).x * (hit.material.ambient).x;
    color.y += (scene.ambient_light).y * (hit.material.ambient).y;
    color.z += (scene.ambient_light).z * (hit.material.ambient).z;
    
    //Doing the above is mathematically correct, but if you then try to raytrace from that ray, you may hit the same object you are reflecting off.
    // One way to fight that problem is to push the ray positin a small amount away from the surface to make sure it misses it. 
    //You can do that for example like this: ReflectRayLocation = ReflectRayLocation + ReflectRayDirection * 0.01

    if(hit.material.is_mirror && depth <= scene.max_recursion_depth){ //reflectance bak doğru hesaplanmış mı diye
        Vec3f reflection_ray_direction = ray.getDirection() *(-1) + hit.normal*(2*dot(ray.getDirection(),hit.normal));
        reflection_ray_direction = normalize(reflection_ray_direction);
        Ray reflection_ray = Ray(hit.intersection_point + reflection_ray_direction * scene.shadow_ray_epsilon , reflection_ray_direction);

        Vec3f computed_color={0,0,0};
         //gereksiz gibi constructorda vardı zaten
        computed_color = computed_color + compute_color(reflection_ray , depth+1,scene);
        color.x += computed_color.x * hit.material.mirror.x ;
        color.y += computed_color.y * hit.material.mirror.y ;
        color.z += computed_color.z * hit.material.mirror.z ;
        color.x = color.x > 255 ? 255 : color.x;
        color.y = color.y >255 ? 255 :color.y;
        color.z = color.z > 255 ? 255 : color.z;
    }
    

    for(int i=0; i<scene.point_lights.size(); i++){

        if(!inShadow(hit.intersection_point,scene.point_lights[i],scene)){
                Vec3f objToLight = normalize(scene.point_lights[i].position-hit.intersection_point);
                float clampedDotProduct = fmax(0,dot(objToLight,normalize(hit.normal)));
                float diffuseRed = scene.point_lights[i].intensity.x * clampedDotProduct * hit.material.diffuse.x / pow(calculateDistance(scene.point_lights[i].position,hit.intersection_point),2);
                float diffuseGreen = scene.point_lights[i].intensity.y * clampedDotProduct * hit.material.diffuse.y / pow(calculateDistance(scene.point_lights[i].position,hit.intersection_point),2);
                float diffuseBlue = scene.point_lights[i].intensity.z * clampedDotProduct * hit.material.diffuse.z / pow(calculateDistance(scene.point_lights[i].position,hit.intersection_point),2);
                color.x+=diffuseRed;
                color.y+=diffuseGreen;
                color.z+=diffuseBlue;

                Vec3f viewDirection =  ray.getDirection() * (-1.0);
                Vec3f halfwayVector = normalize(objToLight + viewDirection);
                float clampedDotProduct2 = fmax(0, dot(halfwayVector, hit.normal)); //normalize et çok önemli
                float specularRed = hit.material.specular.x * pow(clampedDotProduct2, hit.material.phong_exponent)* scene.point_lights[i].intensity.x  / pow(calculateDistance(scene.point_lights[i].position,hit.intersection_point),2);
                float specularGreen = hit.material.specular.y * pow(clampedDotProduct2, hit.material.phong_exponent) * scene.point_lights[i].intensity.y  / pow(calculateDistance(scene.point_lights[i].position,hit.intersection_point),2);
                float specularBlue = hit.material.specular.z * pow(clampedDotProduct2, hit.material.phong_exponent) * scene.point_lights[i].intensity.z / pow(calculateDistance(scene.point_lights[i].position,hit.intersection_point),2);
                color.x+=specularRed;
                color.y+=specularGreen;
                color.z+=specularBlue;

        }
        

    }
    return color;
}


Vec3f compute_color(Ray ray,int depth, Scene& scene){ //en son mainin içinde Vec3i'ye çevirip clample
        if(depth > scene.max_recursion_depth){
            return {0,0,0};
        }

        if(ray.find_closest_hit(scene).is_intersected == true){
            return applyShading(ray,depth,ray.find_closest_hit(scene),scene);
        }

        else if(depth==0) { //Vec3i Vec3f'ye dönüştürüldü. Umarım loss olmamıştır.
            Vec3f background;
            background.x = scene.background_color.x;
            background.y = scene.background_color.y;
            background.z = scene.background_color.z;
            return {background};
        }

        else return {0,0,0};
}




    //compute_color(Ray r, int depth) içinde,, o rayi scenedeki tüm objelerle kesiştirmeye çalışacağız 
    //(all spheres loop, all triangles loop,...) hit record tutacağız bir yandan ve bu kesişimin (hitrecord.is_intersected==true) ise
    // ve bir önceki kesişimden daha yakınsa(Hitrecord.t < current_min_t) hit recordu güncelle.
    //tüm objeleri döndüğünde elinde kalan hitrecord en yakın hit olmuş olacak.
    //daha sonra bu ray'i ve hitrecordunu shading fonskiyonuna gönder. shading fonksiyonu recursive, depthi artırarar compute_color'u geri çağırabilir.


    //mainin içinde, cameradan tüm pixelleri dön, cameratopixel'i çağır, rayi elde et. Depthi 0 vererek compute_color'u her ray için çağır.
    //Compute color bir renk dönecek, onu da 0,255 arasına sıkıştır. SON.







int main()
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;
    scene.loadFromXml("inputs/monkey.xml");


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
