// AOBench C++ AMP / Windows
//
// original AOBench: http://code.google.com/p/aobench/
// C++ AMP RNG Library: http://amprng.codeplex.com/


#define AO_USE_AMP_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <assert.h>


#define _USE_MATH_DEFINES
#include <amp.h>
#include <amp_math.h>
#include <math.h>

#include <windows.h>
#include <mmsystem.h>

#include "amp_tinymt_rng.h"


#define WIDTH		256
#define HEIGHT	   256
#define NSUBSAMPLES  2
#define NAO_SAMPLES  8

#if defined( AO_USE_AMP_ )
#define AMP_SUFFIX restrict(amp)
#else
#define AMP_SUFFIX
#endif

using namespace concurrency;
using namespace concurrency::precise_math;

typedef struct _vec
{
	double x;
	double y;
	double z;
} vec;


typedef struct _Isect
{
	double t;
	vec	p;
	vec	n;
	int	hit; 
} Isect;

typedef struct _Sphere
{
	vec	center;
	double radius;

} Sphere;

typedef struct _Plane
{
	vec	p;
	vec	n;

} Plane;

typedef struct _Ray
{
	vec	org;
	vec	dir;
} Ray;

struct Scene{
	Sphere spheres[3];
	Plane  plane;
};

double vdot(vec v0, vec v1) AMP_SUFFIX
{
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

void vcross(vec *c, vec v0, vec v1) AMP_SUFFIX
{
	
	c->x = v0.y * v1.z - v0.z * v1.y;
	c->y = v0.z * v1.x - v0.x * v1.z;
	c->z = v0.x * v1.y - v0.y * v1.x;
}

void vnormalize(vec *c) AMP_SUFFIX
{
	double length = sqrt(vdot((*c), (*c)));

	if (fabs(length) > 1.0e-17) {
		c->x /= length;
		c->y /= length;
		c->z /= length;
	}
}

void
ray_sphere_intersect(Isect *isect, const Ray *ray, const Sphere *sphere) AMP_SUFFIX
{
	vec rs;

	rs.x = ray->org.x - sphere->center.x;
	rs.y = ray->org.y - sphere->center.y;
	rs.z = ray->org.z - sphere->center.z;

	double B = vdot(rs, ray->dir);
	double C = vdot(rs, rs) - sphere->radius * sphere->radius;
	double D = B * B - C;

	if (D > 0.0) {
		double t = -B - sqrt(D);
		
		if ((t > 0.0) && (t < isect->t)) {
			isect->t = t;
			isect->hit = 1;
			
			isect->p.x = ray->org.x + ray->dir.x * t;
			isect->p.y = ray->org.y + ray->dir.y * t;
			isect->p.z = ray->org.z + ray->dir.z * t;

			isect->n.x = isect->p.x - sphere->center.x;
			isect->n.y = isect->p.y - sphere->center.y;
			isect->n.z = isect->p.z - sphere->center.z;

			vnormalize(&(isect->n));
		}
	}
}

void
ray_plane_intersect(Isect *isect, const Ray *ray, const Plane *plane) AMP_SUFFIX
{
	double d = -vdot(plane->p, plane->n);
	double v = vdot(ray->dir, plane->n);

	if (fabs(v) < 1.0e-17) return;

	double t = -(vdot(ray->org, plane->n) + d) / v;

	if ((t > 0.0) && (t < isect->t)) {
		isect->t = t;
		isect->hit = 1;
		
		isect->p.x = ray->org.x + ray->dir.x * t;
		isect->p.y = ray->org.y + ray->dir.y * t;
		isect->p.z = ray->org.z + ray->dir.z * t;

		isect->n = plane->n;
	}
}

void
orthoBasis(vec *basis, vec n) AMP_SUFFIX
{
	basis[2] = n;
	basis[1].x = 0.0; basis[1].y = 0.0; basis[1].z = 0.0;

	if ((n.x < 0.6) && (n.x > -0.6)) {
		basis[1].x = 1.0;
	} else if ((n.y < 0.6) && (n.y > -0.6)) {
		basis[1].y = 1.0;
	} else if ((n.z < 0.6) && (n.z > -0.6)) {
		basis[1].z = 1.0;
	} else {
		basis[1].x = 1.0;
	}

	vcross(&basis[0], basis[1], basis[2]);
	vnormalize(&basis[0]);

	vcross(&basis[1], basis[2], basis[0]);
	vnormalize(&basis[1]);
}

#if defined( AO_USE_AMP_ )
void ambient_occlusion_amp(vec *col, const Isect *isect, const Scene& scene, const tinymt_collection<2>& randomClass, concurrency::index<2>& idx ) AMP_SUFFIX
{
	int	i, j;
	int	ntheta = NAO_SAMPLES;
	int	nphi   = NAO_SAMPLES;
	double eps = 0.0001;

	vec p;

	p.x = isect->p.x + eps * isect->n.x;
	p.y = isect->p.y + eps * isect->n.y;
	p.z = isect->p.z + eps * isect->n.z;

	vec basis[3];
	orthoBasis(basis, isect->n);

	double occlusion = 0.0;

	for (j = 0; j < ntheta; j++) {
		for (i = 0; i < nphi; i++) {
			double theta = sqrt(randomClass[idx].next_single());
			double phi   = 2.0 * M_PI * randomClass[idx].next_single();

			double x = cos(phi) * theta;
			double y = sin(phi) * theta;
			double z = sqrt(1.0 - theta * theta);

			// local -> global
			double rx = x * basis[0].x + y * basis[1].x + z * basis[2].x;
			double ry = x * basis[0].y + y * basis[1].y + z * basis[2].y;
			double rz = x * basis[0].z + y * basis[1].z + z * basis[2].z;

			Ray ray;

			ray.org = p;
			ray.dir.x = rx;
			ray.dir.y = ry;
			ray.dir.z = rz;

			Isect occIsect;
			occIsect.t   = 1.0e+17;
			occIsect.hit = 0;

			ray_sphere_intersect(&occIsect, &ray, &scene.spheres[0]); 
			ray_sphere_intersect(&occIsect, &ray, &scene.spheres[1]); 
			ray_sphere_intersect(&occIsect, &ray, &scene.spheres[2]); 
			ray_plane_intersect (&occIsect, &ray, &scene.plane); 

			if (occIsect.hit) occlusion += 1.0;
			
		}
	}

	occlusion = (ntheta * nphi - occlusion) / (double)(ntheta * nphi);

	col->x = occlusion;
	col->y = occlusion;
	col->z = occlusion;
}

#else
double drand48( void ){
	return (double)rand() / (double)RAND_MAX;
}

void ambient_occlusion(vec *col, const Isect *isect, const Scene& scene )
{
	int	i, j;
	int	ntheta = NAO_SAMPLES;
	int	nphi   = NAO_SAMPLES;
	double eps = 0.0001;

	vec p;

	p.x = isect->p.x + eps * isect->n.x;
	p.y = isect->p.y + eps * isect->n.y;
	p.z = isect->p.z + eps * isect->n.z;

	vec basis[3];
	orthoBasis(basis, isect->n);

	double occlusion = 0.0;

	for (j = 0; j < ntheta; j++) {
		for (i = 0; i < nphi; i++) {
			double theta = sqrt(drand48());
			double phi   = 2.0 * M_PI * drand48();

			double x = cos(phi) * theta;
			double y = sin(phi) * theta;
			double z = sqrt(1.0 - theta * theta);

			// local -> global
			double rx = x * basis[0].x + y * basis[1].x + z * basis[2].x;
			double ry = x * basis[0].y + y * basis[1].y + z * basis[2].y;
			double rz = x * basis[0].z + y * basis[1].z + z * basis[2].z;

			Ray ray;

			ray.org = p;
			ray.dir.x = rx;
			ray.dir.y = ry;
			ray.dir.z = rz;

			Isect occIsect;
			occIsect.t   = 1.0e+17;
			occIsect.hit = 0;

			ray_sphere_intersect(&occIsect, &ray, &scene.spheres[0]); 
			ray_sphere_intersect(&occIsect, &ray, &scene.spheres[1]); 
			ray_sphere_intersect(&occIsect, &ray, &scene.spheres[2]); 
			ray_plane_intersect (&occIsect, &ray, &scene.plane); 

			if (occIsect.hit) occlusion += 1.0;
			
		}
	}

	occlusion = (ntheta * nphi - occlusion) / (double)(ntheta * nphi);

	col->x = occlusion;
	col->y = occlusion;
	col->z = occlusion;
}
#endif

unsigned char
clamp(double f)
{
  int i = (int)(f * 255.5);

  if (i < 0) i = 0;
  if (i > 255) i = 255;

  return (unsigned char)i;
}

void
render(unsigned char *img, const Scene& scene, int w, int h, int nsubsamples)
{

	struct RGB{
		double r_;
		double g_;
		double b_;
	};

	std::vector<RGB> imgBuffer(w * h);

	concurrency::extent<2> ext(w, h);
	concurrency::array_view< RGB, 2 > view( ext, imgBuffer);
	concurrency::index<2> idx(w, h);
	tinymt_collection<2>	tinyRand( ext );
	
	#if defined( AO_USE_AMP_ )
	view.discard_data();
	auto loopBody = [=](concurrency::index<2> idx) AMP_SUFFIX {

		for (int v = 0; v < nsubsamples; v++) {
			for (int u = 0; u < nsubsamples; u++) {
				double px = (idx[0] + (u / (double)nsubsamples) - (w / 2.0)) / (w / 2.0);
				double py = -(idx[1] + (v / (double)nsubsamples) - (h / 2.0)) / (h / 2.0);

				Ray ray;

				ray.org.x = 0.0;
				ray.org.y = 0.0;
				ray.org.z = 0.0;

				ray.dir.x = px;
				ray.dir.y = py;
				ray.dir.z = -1.0;
				vnormalize(&(ray.dir));

				Isect isect;
				isect.t = 1.0e+17;
				isect.hit = 0;

				ray_sphere_intersect(&isect, &ray, &scene.spheres[0]);
				ray_sphere_intersect(&isect, &ray, &scene.spheres[1]);
				ray_sphere_intersect(&isect, &ray, &scene.spheres[2]);
				ray_plane_intersect(&isect, &ray, &scene.plane);

				if (isect.hit) {
					vec col;
					ambient_occlusion_amp(&col, &isect, scene, tinyRand, idx );

					view[idx].r_ += col.x;
					view[idx].g_ += col.y;
					view[idx].b_ += col.z;
				}

			}
		}
	};

	concurrency::parallel_for_each(ext, loopBody);
	view.synchronize();
	#else

	for( idx[0] = 0; idx[0] < ext[0]; ++idx[0] ){
		for (idx[1] = 0; idx[1] < ext[1]; ++idx[1]){

			for (int v = 0; v < nsubsamples; v++) {
				for (int u = 0; u < nsubsamples; u++) {
					double px = (idx[0] + (u / (double)nsubsamples) - (w / 2.0)) / (w / 2.0);
					double py = -(idx[1] + (v / (double)nsubsamples) - (h / 2.0)) / (h / 2.0);

					Ray ray;

					ray.org.x = 0.0;
					ray.org.y = 0.0;
					ray.org.z = 0.0;

					ray.dir.x = px;
					ray.dir.y = py;
					ray.dir.z = -1.0;
					vnormalize(&(ray.dir));

					Isect isect;
					isect.t = 1.0e+17;
					isect.hit = 0;

					ray_sphere_intersect(&isect, &ray, &scene.spheres[0]);
					ray_sphere_intersect(&isect, &ray, &scene.spheres[1]);
					ray_sphere_intersect(&isect, &ray, &scene.spheres[2]);
					ray_plane_intersect(&isect, &ray, &scene.plane);

					if (isect.hit) {
						vec col;
						ambient_occlusion(&col, &isect, scene);

						view[idx].r_ += col.x;
						view[idx].g_ += col.y;
						view[idx].b_ += col.z;
					}
				}
			}
			
		}
	}
	
	#endif

	for (idx[0] = 0; idx[0] < ext[0]; idx[0]++) {
		for (idx[1] = 0; idx[1] < ext[1]; idx[1]++) {
			view[idx].r_ /= (double)(nsubsamples * nsubsamples);
			view[idx].g_ /= (double)(nsubsamples * nsubsamples);
			view[idx].b_ /= (double)(nsubsamples * nsubsamples);
		
			img[3 * (idx[1] * w + idx[0]) + 0] = clamp(view[idx].r_);
			img[3 * (idx[1] * w + idx[0]) + 1] = clamp(view[idx].g_);
			img[3 * (idx[1] * w + idx[0]) + 2] = clamp(view[idx].b_);
		}
	}

}

void
init_scene( Scene& scene )
{
	scene.spheres[0].center.x = -2.0;
	scene.spheres[0].center.y =  0.0;
	scene.spheres[0].center.z = -3.5;
	scene.spheres[0].radius = 0.5;

	scene.spheres[1].center.x = -0.5;
	scene.spheres[1].center.y =  0.0;
	scene.spheres[1].center.z = -3.0;
	scene.spheres[1].radius = 0.5;

	scene.spheres[2].center.x =  1.0;
	scene.spheres[2].center.y =  0.0;
	scene.spheres[2].center.z = -2.2;
	scene.spheres[2].radius = 0.5;

	scene.plane.p.x = 0.0;
	scene.plane.p.y = -0.5;
	scene.plane.p.z = 0.0;

	scene.plane.n.x = 0.0;
	scene.plane.n.y = 1.0;
	scene.plane.n.z = 0.0;

}

void
saveppm(const char *fname, int w, int h, unsigned char *img)
{
	FILE *fp;

	fp = fopen(fname, "wb");
	assert(fp);

	fprintf(fp, "P6\n");
	fprintf(fp, "%d %d\n", w, h);
	fprintf(fp, "255\n");
	fwrite(img, w * h * 3, 1, fp);
	fclose(fp);
}

void showAccelerator( void ){
	std::vector<accelerator> accs = accelerator::get_all();

	std::cout << "=======Devices======" << std::endl;
	for( auto itr = accs.begin(); itr != accs.end(); ++itr ){
		std::wcout <<  itr->get_description() << std::endl;
		std::cout << "Shared memory support:";
		if( itr->supports_cpu_shared_memory ){
			std::cout << "[O]" << std::endl;
		}else{
			std::cout << "[X]" << std::endl;
		}
		std::cout << "---------------------" << std::endl;
	}

	accelerator defaultAcce;
	std::wcout << "System Default Accelerator:" << defaultAcce.get_description() << std::endl;
	
}

int main(int argc, char* argv[])
{
	unsigned char *img = (unsigned char *)malloc(WIDTH * HEIGHT * 3);


	Scene scene;

	init_scene( scene );

	int loopMax = 10;

	#if defined( AO_USE_AMP_ )
	showAccelerator();
	#endif
	
	DWORD startTime = timeGetTime();
	
	for( int i = 0; i < loopMax; ++i ){
		render(img, scene, WIDTH, HEIGHT, NSUBSAMPLES);
		saveppm("ao.ppm", WIDTH, HEIGHT, img);
		printf( "Loop %d\n", i );
	}
	
	DWORD endTime = timeGetTime();

	printf( "total : %i ms\n", ( endTime - startTime ) );
	
	
	return 0;
}

