#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <random>
#include <omp.h>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#ifndef M_PI
#define M_PI 3.14159265358979323856
#endif

static std::default_random_engine engine[32];
static std::uniform_real_distribution<double> uniform(0, 1);

double sqr(double x) { return x * x; };

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	double norm2() const {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
	double operator[](int i) const { return data[i]; };
	double& operator[](int i) { return data[i]; };
	double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray {
public:
	Ray(const Vector& origin, const Vector& unit_direction) : O(origin), u(unit_direction) {};
	Vector O, u;
};

class Object {
public:
	Object(const Vector& albedo, bool mirror = false, bool transparent = false) : albedo(albedo), mirror(mirror), transparent(transparent) {};

	virtual bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const = 0;

	Vector albedo;
	bool mirror, transparent;
};

class Sphere : public Object {
public:
	Sphere(const Vector& center, double radius, const Vector& albedo, bool mirror = false, bool transparent = false) : ::Object(albedo, mirror, transparent), C(center), R(radius) {};

	// returns true iif there is an intersection between the ray and the sphere
	// if there is an intersection, also computes the point of intersection P, 
	// t>=0 the distance between the ray origin and P (i.e., the parameter along the ray)
	// and the unit normal N
	bool intersect(const Ray& ray, Vector& P, double &t, Vector& N) const {
		 // TODO (lab 1) : compute the intersection (just true/false at the begining of lab 1, then P, t and N as well)
		Vector OC = ray.O - C;

		double b = dot(ray.u, OC);
		double norme2_OC = dot(OC, OC);

		double delta = b * b - (norme2_OC - R * R);

		if (delta < 0) return false;

		double sqrt_delta = sqrt(delta);

		double t_moins = -b - sqrt_delta;
		double t_plus = -b + sqrt_delta;

		if (t_plus < 0) return false;

		t = (t_moins >= 0) ? t_moins : t_plus;

		P = ray.O + t * ray.u;

		N = P - C;
		N.normalize();

		return true;
	}

	double R;
	Vector C;
};




// Class only used in labs 3 and 4 
class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1) {
		vtx[0] = vtxi; vtx[1] = vtxj; vtx[2] = vtxk;
		uv[0] = uvi; uv[1] = uvj; uv[2] = uvk;
		n[0] = ni; n[1] = nj; n[2] = nk;
		this->group = group;
	};
	int vtx[3]; // indices within the vertex coordinates array
	int uv[3];  // indices within the uv coordinates array
	int n[3];   // indices within the normals array
	int group;  // face group
};

// Class only used in labs 3 and 4 


class BoundingBox {
public:
	BoundingBox() : Bmin(1e9, 1e9, 1e9), Bmax(-1e9, -1e9, -1e9) {}

	void extend(const Vector& P) {
		for (int i = 0; i < 3; i++) {
			Bmin[i] = std::min(Bmin[i], P[i]);
			Bmax[i] = std::max(Bmax[i], P[i]);
		}
	}

	bool intersect(const Ray& ray) const {
		double t_xmin = (Bmin[0] - ray.O[0]) / ray.u[0];
		double t_xmax = (Bmax[0] - ray.O[0]) / ray.u[0];
		if (t_xmin > t_xmax) std::swap(t_xmin, t_xmax);

		double t_ymin = (Bmin[1] - ray.O[1]) / ray.u[1];
		double t_ymax = (Bmax[1] - ray.O[1]) / ray.u[1];
		if (t_ymin > t_ymax) std::swap(t_ymin, t_ymax);

		double t_zmin = (Bmin[2] - ray.O[2]) / ray.u[2];
		double t_zmax = (Bmax[2] - ray.O[2]) / ray.u[2];
		if (t_zmin > t_zmax) std::swap(t_zmin, t_zmax);

		double t_entree = std::max(t_xmin, std::max(t_ymin, t_zmin));
		double t_sortie = std::min(t_xmax, std::min(t_ymax, t_zmax));

		return (t_sortie >= t_entree && t_sortie >= 0.0);
	}

	Vector Bmin, Bmax;
};


class TriangleMesh : public Object {
public:
	TriangleMesh(const Vector& albedo, bool mirror = false, bool transparent = false) : ::Object(albedo, mirror, transparent) {};

	// first scale and then translate the current object
	void scale_translate(double s, const Vector& t) {

		bbox = BoundingBox();
		for (int i = 0; i < vertices.size(); i++) {
			vertices[i] = vertices[i] * s + t;
			bbox.extend(vertices[i]);
		}

	}

	// read an .obj file
	void readOBJ(const char* obj) {
		std::ifstream f(obj);
		if (!f) return;

		std::map<std::string, int> mtls;
		int curGroup = -1, maxGroup = -1;

		// OBJ indices are 1-based and can be negative (relative), this normalizes them
		auto resolveIdx = [](int i, int size) {
			return i < 0 ? size + i : i - 1;
		};

		auto setFaceVerts = [&](TriangleIndices& t, int i0, int i1, int i2) {
			t.vtx[0] = resolveIdx(i0, vertices.size());
			t.vtx[1] = resolveIdx(i1, vertices.size());
			t.vtx[2] = resolveIdx(i2, vertices.size());
		};
		auto setFaceUVs = [&](TriangleIndices& t, int j0, int j1, int j2) {
			t.uv[0] = resolveIdx(j0, uvs.size());
			t.uv[1] = resolveIdx(j1, uvs.size());
			t.uv[2] = resolveIdx(j2, uvs.size());
		};
		auto setFaceNormals = [&](TriangleIndices& t, int k0, int k1, int k2) {
			t.n[0] = resolveIdx(k0, normals.size());
			t.n[1] = resolveIdx(k1, normals.size());
			t.n[2] = resolveIdx(k2, normals.size());
		};

		std::string line;
		while (std::getline(f, line)) {
			// Trim trailing whitespace
			line.erase(line.find_last_not_of(" \r\t\n") + 1);
			if (line.empty()) continue;

			const char* s = line.c_str();

			if (line.rfind("usemtl ", 0) == 0) {
				std::string matname = line.substr(7);
				auto result = mtls.emplace(matname, maxGroup + 1);
				if (result.second) {
					curGroup = ++maxGroup;
				} else {
					curGroup = result.first->second;
				}
			} else if (line.rfind("vn ", 0) == 0) {
				Vector v;
				sscanf(s, "vn %lf %lf %lf", &v[0], &v[1], &v[2]);
				normals.push_back(v);
			} else if (line.rfind("vt ", 0) == 0) {
				Vector v;
				sscanf(s, "vt %lf %lf", &v[0], &v[1]);
				uvs.push_back(v);
			} else if (line.rfind("v ", 0) == 0) {
				Vector pos, col;
				if (sscanf(s, "v %lf %lf %lf %lf %lf %lf", &pos[0], &pos[1], &pos[2], &col[0], &col[1], &col[2]) == 6) {
					for (int i = 0; i < 3; i++) col[i] = std::min(1.0, std::max(0.0, col[i]));
					vertexcolors.push_back(col);
				} else {
					sscanf(s, "v %lf %lf %lf", &pos[0], &pos[1], &pos[2]);
				}
				vertices.push_back(pos);
			}
			else if (line[0] == 'f') {
				int i[4], j[4], k[4], offset, nn;
				const char* cur = s + 1;
				TriangleIndices t;
				t.group = curGroup;

				// Try each face format: v/vt/vn, v/vt, v//vn, v
				if ((nn = sscanf(cur, "%d/%d/%d %d/%d/%d %d/%d/%d%n", &i[0], &j[0], &k[0], &i[1], &j[1], &k[1], &i[2], &j[2], &k[2], &offset)) == 9) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceUVs(t, j[0], j[1], j[2]); 
					setFaceNormals(t, k[0], k[1], k[2]);
				} else if ((nn = sscanf(cur, "%d/%d %d/%d %d/%d%n", &i[0], &j[0], &i[1], &j[1], &i[2], &j[2], &offset)) == 6) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceUVs(t, j[0], j[1], j[2]);
				} else if ((nn = sscanf(cur, "%d//%d %d//%d %d//%d%n", &i[0], &k[0], &i[1], &k[1], &i[2], &k[2], &offset)) == 6) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceNormals(t, k[0], k[1], k[2]);
				} else if ((nn = sscanf(cur, "%d %d %d%n", &i[0], &i[1], &i[2], &offset)) == 3) {
					setFaceVerts(t, i[0], i[1], i[2]);
				}
				else continue;

				indices.push_back(t);
				cur += offset;

				// Fan triangulation for polygon faces (4+ vertices)
				while (*cur && *cur != '\n') {
					TriangleIndices t2;
					t2.group = curGroup;
					if ((nn = sscanf(cur, " %d/%d/%d%n", &i[3], &j[3], &k[3], &offset)) == 3) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceUVs(t2, j[0], j[2], j[3]); 
						setFaceNormals(t2, k[0], k[2], k[3]);
					} else if ((nn = sscanf(cur, " %d/%d%n", &i[3], &j[3], &offset)) == 2) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceUVs(t2, j[0], j[2], j[3]);
					} else if ((nn = sscanf(cur, " %d//%d%n", &i[3], &k[3], &offset)) == 2) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceNormals(t2, k[0], k[2], k[3]);
					} else if ((nn = sscanf(cur, " %d%n", &i[3], &offset)) == 1) {
						setFaceVerts(t2, i[0], i[2], i[3]);
					} else { 
						cur++; 
						continue; 
					}

					indices.push_back(t2);
					cur += offset;
					i[2] = i[3]; j[2] = j[3]; k[2] = k[3];
				}
			}
		}
		bbox = BoundingBox();

		for (int i = 0; i < vertices.size(); i++) {
			bbox.extend(vertices[i]);

		}
	}


	int construire_bvh_recursif(int debut, int fin) {
		int idx_noeud = (int)bvh_min.size();

		Vector boite_min_locale(1e18, 1e18, 1e18);
		Vector boite_max_locale(-1e18, -1e18, -1e18);

		for (int k = debut; k < fin; k++) {
			for (int s = 0; s < 3; s++) {
				const Vector& sommet = vertices[indices[k].vtx[s]];
				for (int axe = 0; axe < 3; axe++) {
					if (sommet[axe] < boite_min_locale[axe]) boite_min_locale[axe] = sommet[axe];
					if (sommet[axe] > boite_max_locale[axe]) boite_max_locale[axe] = sommet[axe];
				}
			}
		}

		bvh_min.push_back(boite_min_locale);
		bvh_max.push_back(boite_max_locale);
		bvh_gauche.push_back(-1);
		bvh_droite.push_back(-1);
		bvh_debut.push_back(debut);
		bvh_fin.push_back(fin);

		if (fin - debut < 5) return idx_noeud;

		Vector diag = boite_max_locale - boite_min_locale;
		int axe_coupe = 0;
		if (diag[1] > diag[0]) axe_coupe = 1;
		if (diag[2] > diag[axe_coupe]) axe_coupe = 2;

		double milieu = 0.5 * (boite_min_locale[axe_coupe] + boite_max_locale[axe_coupe]);

		int pivot = debut;
		for (int k = debut; k < fin; k++) {
			Vector bary = (vertices[indices[k].vtx[0]] + vertices[indices[k].vtx[1]] + vertices[indices[k].vtx[2]]) / 3.0;
			if (bary[axe_coupe] < milieu) {
				std::swap(indices[k], indices[pivot]);
				pivot++;
			}
		}

		if (pivot <= debut || pivot >= fin) pivot = (debut + fin) / 2;

		int idx_g = construire_bvh_recursif(debut, pivot);
		int idx_d = construire_bvh_recursif(pivot, fin);

		bvh_gauche[idx_noeud] = idx_g;
		bvh_droite[idx_noeud] = idx_d;

		return idx_noeud;
	}

	void construire_bvh() {
		bvh_min.clear();
		bvh_max.clear();
		bvh_gauche.clear();
		bvh_droite.clear();
		bvh_debut.clear();
		bvh_fin.clear();

		if (!indices.empty()) construire_bvh_recursif(0, (int)indices.size());
	}

	bool intersecte_boite(const Ray& ray, int idx_noeud) const {
		double tmin = -1e18;
		double tmax = 1e18;

		for (int axe = 0; axe < 3; axe++) {
			double t0 = (bvh_min[idx_noeud][axe] - ray.O[axe]) / ray.u[axe];
			double t1 = (bvh_max[idx_noeud][axe] - ray.O[axe]) / ray.u[axe];
			if (t0 > t1) std::swap(t0, t1);
			if (t0 > tmin) tmin = t0;
			if (t1 < tmax) tmax = t1;
		}

		return !(tmax < 0 || tmin > tmax);
	}

	bool intersecte_bvh(const Ray& ray, Vector& P, double& t, Vector& N, int idx_noeud) const {
		if (!intersecte_boite(ray, idx_noeud)) return false;

		if (bvh_gauche[idx_noeud] == -1) {
			bool trouve = false;

			for (int k = bvh_debut[idx_noeud]; k < bvh_fin[idx_noeud]; k++) {
				const Vector& A = vertices[indices[k].vtx[0]];
				const Vector& B = vertices[indices[k].vtx[1]];
				const Vector& C = vertices[indices[k].vtx[2]];

				Vector e1 = B - A;
				Vector e2 = C - A;
				Vector normale_tri = cross(e1, e2);

				double denom = dot(ray.u, normale_tri);
				if (std::fabs(denom) < 1e-10) continue;

				Vector AO = A - ray.O;
				Vector AOxu = cross(AO, ray.u);

				double beta = dot(e2, AOxu) / denom;
				double gamma = -dot(e1, AOxu) / denom;
				double alpha = 1.0 - beta - gamma;
				double t_local = dot(AO, normale_tri) / denom;

				if (alpha < 0 || beta < 0 || gamma < 0 || t_local < 1e-6) continue;

				if (t_local < t) {
					t = t_local;
					P = ray.O + t * ray.u;
					N = normale_tri;
					N.normalize();
					if (dot(ray.u, N) > 0) N = -1.0 * N;
					trouve = true;
				}
			}

			return trouve;
		}

		bool hit_gauche = intersecte_bvh(ray, P, t, N, bvh_gauche[idx_noeud]);
		bool hit_droite = intersecte_bvh(ray, P, t, N, bvh_droite[idx_noeud]);
		return hit_gauche || hit_droite;
	}
	

	// TODO ray-mesh intersection (labs 3 and 4)

	bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const {

		if (bvh_min.empty()) return false;
		t = 1e18;
		return intersecte_bvh(ray, P, t, N, 0);

	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	std::vector<Vector> bvh_min, bvh_max;

	std::vector<int> bvh_gauche, bvh_droite;

	std::vector<int> bvh_debut, bvh_fin;
	BoundingBox bbox;
};




class Scene {
public:
	Scene() {};
	void addObject(const Object* obj) {
		objects.push_back(obj);
	}

	// returns true iif there is an intersection between the ray and any object in the scene
    // if there is an intersection, also computes the point of the *nearest* intersection P, 
    // t>=0 the distance between the ray origin and P (i.e., the parameter along the ray)
    // and the unit normal N. 
	// Also returns the index of the object within the std::vector objects in object_id
	bool intersect(const Ray& ray, Vector& P, double& t, Vector& N, int &object_id) const  {

		t = 1e18;
		bool intersection_trouvee = false;
		for (int i = 0; i < (int)objects.size(); i++) {
			Vector P_local, N_local;
			double t_local;
			if (objects[i]->intersect(ray, P_local, t_local, N_local) && t_local < t) {
				t = t_local;
				P = P_local;
				N = N_local;
				object_id = i;
				intersection_trouvee = true;
			}
		}
		return intersection_trouvee;

	}


	// return the radiance (color) along ray
	Vector getColor(const Ray& ray, int recursion_depth) {
		if (recursion_depth >= max_light_bounce) return Vector(0, 0, 0);

		Vector P, N;
		double t;
		int object_id;

		if (intersect(ray, P, t, N, object_id)) {
			Vector vers_lumiere = light_position - P;
			Vector dir_lumiere = vers_lumiere / vers_lumiere.norm();

			Ray rayon_ombre(P + 1e-6 * N, dir_lumiere);
			Vector P_ombre, N_ombre;
			double t_ombre;
			int id_ombre;

			Vector couleur_directe(0, 0, 0);

			if (intersect(rayon_ombre, P_ombre, t_ombre, N_ombre, id_ombre)) {
				if (t_ombre >= vers_lumiere.norm() - 1e-6) {
					double attenuation = light_intensity / (4 * M_PI * vers_lumiere.norm2());
					Vector couleur_matiere = objects[object_id]->albedo / M_PI;
					double cos_theta = std::max(0.0, dot(N, dir_lumiere));
					couleur_directe = attenuation * couleur_matiere * cos_theta;
				}
			}
			else {
				double attenuation = light_intensity / (4 * M_PI * vers_lumiere.norm2());
				Vector couleur_matiere = objects[object_id]->albedo / M_PI;
				double cos_theta = std::max(0.0, dot(N, dir_lumiere));
				couleur_directe = attenuation * couleur_matiere * cos_theta;
			}

			int tid = omp_get_thread_num();
			double r1 = uniform(engine[tid]);
			double r2 = uniform(engine[tid]);

			Vector tangent1;
			if (std::fabs(N[0]) <= std::fabs(N[1]) && std::fabs(N[0]) <= std::fabs(N[2]))
				tangent1 = cross(N, Vector(1, 0, 0));
			else if (std::fabs(N[1]) <= std::fabs(N[2]))
				tangent1 = cross(N, Vector(0, 1, 0));
			else
				tangent1 = cross(N, Vector(0, 0, 1));

			tangent1.normalize();
			Vector tangent2 = cross(N, tangent1);

			double cos_t = std::sqrt(1 - r2);
			double sin_t = std::sqrt(r2);
			double phi = 2 * M_PI * r1;

			Vector dir_rebond =
				sin_t * std::cos(phi) * tangent1 +
				sin_t * std::sin(phi) * tangent2 +
				cos_t * N;
			dir_rebond.normalize();

			Ray rayon_indirect(P + 1e-6 * N, dir_rebond);
			Vector couleur_indirecte = objects[object_id]->albedo * getColor(rayon_indirect, recursion_depth + 1);

			return couleur_directe + couleur_indirecte;
		}

		return Vector(0, 0, 0);
	}

	std::vector<const Object*> objects;

	Vector camera_center, light_position;
	double fov, gamma, light_intensity;
	int max_light_bounce;
};


int main() {
	int W = 512;
	int H = 512;

	for (int i = 0; i<32; i++) {
		engine[i].seed(i);
	}

	Sphere center_sphere(Vector(0, 0, 0), 10., Vector(0.8, 0.8, 0.8));
	Sphere wall_left(Vector(-1000, 0, 0), 940, Vector(0.5, 0.8, 0.1));
	Sphere wall_right(Vector(1000, 0, 0), 940, Vector(0.9, 0.2, 0.3));
	Sphere wall_front(Vector(0, 0, -1000), 940, Vector(0.1, 0.6, 0.7));
	Sphere wall_behind(Vector(0, 0, 1000), 940, Vector(0.8, 0.2, 0.9));
	Sphere ceiling(Vector(0, 1000, 0), 940, Vector(0.3, 0.5, 0.3));
	Sphere floor(Vector(0, -1000, 0), 990, Vector(0.6, 0.5, 0.7));

TriangleMesh chat(Vector(0.8, 0.6, 0.4));

	chat.readOBJ("cat.obj");


	chat.scale_translate(0.6, Vector(0, -10, 0));

	chat.construire_bvh();

	Scene scene;
	scene.camera_center = Vector(0, 0, 55);
	scene.light_position = Vector(-10,20,40);
	scene.light_intensity = 3E7;
	scene.fov = 60 * M_PI / 180.;
	scene.gamma = 2.2;    // TODO (lab 1) : play with gamma ; typically, gamma = 2.2
	scene.max_light_bounce = 5;

	scene.addObject(&chat);

	scene.addObject(&wall_left);
	scene.addObject(&wall_right);
	scene.addObject(&wall_front);
	scene.addObject(&wall_behind);
	scene.addObject(&ceiling);
	scene.addObject(&floor);


	std::vector<unsigned char> image(W * H * 3, 0);

#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector color;

			int tid = omp_get_thread_num();
			int nb_samples = 64;

			for (int s = 0; s < nb_samples; s++) {
				double dx = uniform(engine[tid]) - 0.5;
				double dy = uniform(engine[tid]) - 0.5;

				Vector dir_camera(
					j - W / 2 + 0.5 + dx,
					H / 2 - i - 0.5 + dy,
					-(W / (2 * std::tan(scene.fov / 2)))
				);
				dir_camera.normalize();

				Ray ray(scene.camera_center, dir_camera);
				color = color + scene.getColor(ray, 0);
			}

			color = color / nb_samples;

			image[(i * W + j) * 3 + 0] = std::min(255.0, std::max(0.0, 255.0 * std::pow(color[0] / 255.0, 1.0 / scene.gamma)));
			image[(i * W + j) * 3 + 1] = std::min(255.0, std::max(0.0, 255.0 * std::pow(color[1] / 255.0, 1.0 / scene.gamma)));
			image[(i * W + j) * 3 + 2] = std::min(255.0, std::max(0.0, 255.0 * std::pow(color[2] / 255.0, 1.0 / scene.gamma)));
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}

// modification to commit lab 4 end of the session, 6th  of May 2026