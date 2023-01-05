#include "raylib.h"
#include <raymath.h>
#include "rlgl.h"
#include <math.h>
#include <float.h>
#include <vector>

#include <iostream>
using namespace std;

#if defined(PLATFORM_DESKTOP)
#define GLSL_VERSION            330
#else   // PLATFORM_RPI, PLATFORM_ANDROID, PLATFORM_WEB
#define GLSL_VERSION            100
#endif

#define EPSILON 1.e-6f

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

#pragma region Structs

struct ReferenceFrame {
	Vector3 origin;
	Vector3 i, j, k;
	Quaternion q;
	ReferenceFrame()
	{
		origin = { 0,0,0 };
		i = { 1,0,0 };
		j = { 0,1,0 };
		k = { 0,0,1 };
		q = QuaternionIdentity();
	}
	ReferenceFrame(Vector3 origin, Quaternion q)
	{
		this->q = q;
		this->origin = origin;
		i = Vector3RotateByQuaternion({ 1,0,0 }, q);
		j = Vector3RotateByQuaternion({ 0,1,0 }, q);
		k = Vector3RotateByQuaternion({ 0,0,1 }, q);
	}
	void Translate(Vector3 vect)
	{
		this->origin = Vector3Add(this->origin, vect);
	}
	void RotateByQuaternion(Quaternion qRot)
	{
		q = QuaternionMultiply(qRot, q);
		i = Vector3RotateByQuaternion({ 1,0,0 }, q);
		j = Vector3RotateByQuaternion({ 0,1,0 }, q);
		k = Vector3RotateByQuaternion({ 0,0,1 }, q);
	}
};

struct Quad {
	ReferenceFrame ref;
	Vector3 extents;
};

struct Plane
{
	Vector3 n;
	float d;
};

struct Disk {
	ReferenceFrame ref;
	float radius;
};

struct Sphere {
	ReferenceFrame ref;
	float radius;
};

struct Box {
	ReferenceFrame ref;
	Vector3 extents;
};

struct Cylinder {
	ReferenceFrame ref;
	float halfHeight;
	float radius;
};

struct Capsule {
	ReferenceFrame ref;
	float halfHeight;
	float radius;
};

struct RoundedBox {
	ReferenceFrame ref;
	Vector3 extents;
	float radius;
};

struct Polar
{
	float rho;
	float theta;
};

struct Cylindrical
{
	float rho;
	float theta;
	float y;
};

struct Spherical
{
	float rho;
	float theta;
	float phi;
};

#pragma endregion

#pragma region Conversion

/*CartesianToPolar*/
Polar CartesianToPolar(Vector2 cart, bool keepThetaPositive = true)
{
	Polar polar = { sqrt(cart.x * cart.x + cart.y * cart.y), atan2f(cart.y,cart.x) };
	if (keepThetaPositive && polar.theta < 0)
	{
		polar.theta += 2 * PI;
	}
	return polar;
}

/*PolarToCartesian*/
Vector2 PolarToCartesian(Polar polar)
{
	return Vector2Scale({ cosf(polar.theta),sinf(polar.theta) }, polar.rho);
}

/*CartesianToCylindrical*/
Cylindrical CartesianToCylindrical(Vector3 cart, bool keepThetaPositive = true)
{
	Cylindrical cylindrical = { sqrt(pow(cart.x, 2) + pow(cart.z, 2)), atan2f(cart.x,cart.z), cart.y };
	if (keepThetaPositive && cylindrical.theta < 0)
	{
		cylindrical.theta += 2 * PI;
	}

	return cylindrical;
}

/*CylindricalToCartesian*/
Vector3 CylindricalToCartesian(Cylindrical cylindrical)
{
	return Vector3({ cylindrical.rho * sinf(cylindrical.theta), cylindrical.y, cylindrical.rho * cosf(cylindrical.theta) });
}

/*CartesianToSpherical*/
Spherical CartesianToSpherical(Vector3 cart, bool keepThetaPositive = true)
{
	float rho = sqrt(pow(cart.x, 2) + pow(cart.y, 2) + pow(cart.z, 2));
	float phi = acosf(cart.y / rho);

	if (abs(rho) < EPSILON) phi = 0;
	else phi = acosf(cart.y / rho);

	Spherical spherical = { rho, atan2f(cart.x,cart.z), phi };
	if (keepThetaPositive && spherical.theta < 0)
	{
		spherical.theta += 2 * PI;
	}

	return spherical;
}

/*SphericalToCartesian*/
Vector3 SphericalToCartesian(Spherical spherical)
{
	return Vector3({ spherical.rho * sinf(spherical.phi) * sinf(spherical.theta), spherical.rho * cosf(spherical.phi), spherical.rho * sinf(spherical.phi) * cosf(spherical.theta) });
}

#pragma endregion

#pragma region Camera
void MyUpdateOrbitalCamera(Camera* camera, float deltaTime)
{
	static Spherical sphPos = { 10,PI / 4.f,PI / 4.f };
	Spherical sphSpeed = { 2 ,0.04f,0.04f };
	float rhoMin = 4;
	float rhoMax = 40;
	Vector2 mousePos;
	Spherical sphDelta;

	/*Calcul du vecteur de déplacement de la souris par différence entre la position courante de la souris et la position précédente*/
	static Vector2 prevMousePos = { 0, 0 };
	mousePos = GetMousePosition();
	Vector2 mouseVect = Vector2Subtract(mousePos, prevMousePos);
	prevMousePos = GetMousePosition();
	//cout << "x: " << mouseVect.x << ", y: " << mouseVect.y << "\n";

	/*Calcul du vecteur de déplacement de la caméra en coordonnées sphériques*/
	Spherical cam = CartesianToSpherical(camera->position);
	cam.rho -= GetMouseWheelMove();
	/*Calcul de la nouvelle position de la caméra en coordonnées sphériques*/
	if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON))
	{
		cam.theta -= mouseVect.x / 100;
		cam.phi -= mouseVect.y / 100;
	}
	camera->position = SphericalToCartesian(cam);
}
#pragma endregion

#pragma region Draw

#pragma region quad
//QUAD
void MyDrawPolygonQuad(Quad quad, Color color = LIGHTGRAY)
{
	int numVertex = 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();
	// BEGINNING OF SPACE TRANSFORMATION INDUCED BY THE LOCAL REFERENCE FRAME
	// methods should be called in this order: rlTranslatef, rlRotatef & rlScalef
	// so that transformations occur in the opposite order: scale, then rotation, then translation
	rlPushMatrix();
	//TRANSLATION
	rlTranslatef(quad.ref.origin.x, quad.ref.origin.y, quad.ref.origin.z);
	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(quad.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	//SCALING
	rlScalef(quad.extents.x, 1, quad.extents.z);
	// END OF SPACE TRANSFORMATION INDUCED BY THE LOCAL REFERENCE FRAME
	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);
	rlVertex3f(1, 0, 1);
	rlVertex3f(1, 0, -1);
	rlVertex3f(-1, 0, -1);
	rlVertex3f(1, 0, 1);
	rlVertex3f(-1, 0, -1);
	rlVertex3f(-1, 0, 1);
	rlEnd();
	//EVERY rlPushMatrix method call should be followed by a rlPopMatrix method call
	rlPopMatrix();
}

void MyDrawWireframeQuad(Quad quad, Color color = DARKGRAY)
{
	int numVertex = 10;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();
	rlPushMatrix();
	rlTranslatef(quad.ref.origin.x, quad.ref.origin.y, quad.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(quad.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(quad.extents.x, 1, quad.extents.z);
	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);
	rlVertex3f(1, 0, 1);
	rlVertex3f(1, 0, -1);
	rlVertex3f(1, 0, -1);
	rlVertex3f(-1, 0, -1);
	rlVertex3f(-1, 0, -1);
	rlVertex3f(1, 0, 1);
	rlVertex3f(-1, 0, -1);
	rlVertex3f(-1, 0, 1);
	rlVertex3f(-1, 0, 1);
	rlVertex3f(1, 0, 1);
	rlEnd();
	rlPopMatrix();
}
void MyDrawQuad(Quad quad, bool drawPolygon = true, bool drawWireframe = true,
	Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY)
{
	if (drawPolygon) MyDrawPolygonQuad(quad, polygonColor);
	if (drawWireframe)MyDrawWireframeQuad(quad, wireframeColor);
}

#pragma endregion

#pragma region box

void MyDrawPolygonBox(Box box, Color color = LIGHTGRAY)
{
	int numVertex = 36;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(box.ref.origin.x, box.ref.origin.y, box.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(box.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(box.extents.x, box.extents.y, box.extents.z);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	// Top
	rlVertex3f(-1, 1, -1);
	rlVertex3f(-1, 1, 1);
	rlVertex3f(1, 1, 1);

	rlVertex3f(-1, 1, -1);
	rlVertex3f(1, 1, 1);
	rlVertex3f(1, 1, -1);

	// Bottom
	rlVertex3f(-1, -1, -1);
	rlVertex3f(1, -1, -1);
	rlVertex3f(1, -1, 1);
	rlVertex3f(1, -1, 1);
	rlVertex3f(-1, -1, 1);
	rlVertex3f(-1, -1, -1);

	// Front
	rlVertex3f(-1, 1, 1);
	rlVertex3f(-1, -1, 1);
	rlVertex3f(1, -1, 1);
	rlVertex3f(-1, 1, 1);
	rlVertex3f(1, -1, 1);
	rlVertex3f(1, 1, 1);

	// Back
	rlVertex3f(-1, 1, -1);
	rlVertex3f(1, 1, -1);
	rlVertex3f(1, -1, -1);
	rlVertex3f(-1, 1, -1);
	rlVertex3f(1, -1, -1);
	rlVertex3f(-1, -1, -1);

	// Left
	rlVertex3f(-1, -1, 1);
	rlVertex3f(-1, 1, 1);
	rlVertex3f(-1, 1, -1);
	rlVertex3f(-1, -1, -1);
	rlVertex3f(-1, -1, 1);
	rlVertex3f(-1, 1, -1);

	// Right
	rlVertex3f(1, -1, 1);
	rlVertex3f(1, -1, -1);
	rlVertex3f(1, 1, -1);
	rlVertex3f(1, 1, 1);
	rlVertex3f(1, -1, 1);
	rlVertex3f(1, 1, -1);

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeBox(Box box, Color color = DARKGRAY)
{
	int numVertex = 36;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(box.ref.origin.x, box.ref.origin.y, box.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(box.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(box.extents.x, box.extents.y, box.extents.z);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	// Front
	rlVertex3f(-1, 1, 1);
	rlVertex3f(1, 1, 1);

	rlVertex3f(1, 1, 1);
	rlVertex3f(1, -1, 1);

	rlVertex3f(1, -1, 1);
	rlVertex3f(-1, -1, 1);

	rlVertex3f(-1, -1, 1);
	rlVertex3f(-1, 1, 1);

	rlVertex3f(-1, 1, 1);
	rlVertex3f(1, -1, 1);

	// Top
	rlVertex3f(-1, 1, 1);
	rlVertex3f(-1, 1, -1);

	rlVertex3f(-1, 1, -1);
	rlVertex3f(1, 1, -1);

	rlVertex3f(1, 1, -1);
	rlVertex3f(1, 1, 1);

	rlVertex3f(-1, 1, -1);
	rlVertex3f(1, 1, 1);

	// Back
	rlVertex3f(1, 1, -1);
	rlVertex3f(1, -1, -1);

	rlVertex3f(1, -1, -1);
	rlVertex3f(-1, -1, -1);

	rlVertex3f(-1, -1, -1);
	rlVertex3f(-1, 1, -1);

	rlVertex3f(-1, 1, -1);
	rlVertex3f(1, -1, -1);

	// Bottom
	rlVertex3f(-1, -1, 1);
	rlVertex3f(-1, -1, -1);

	rlVertex3f(1, -1, 1);
	rlVertex3f(1, -1, -1);

	rlVertex3f(1, -1, 1);
	rlVertex3f(-1, -1, -1);

	// Right
	rlVertex3f(1, -1, 1);
	rlVertex3f(1, 1, -1);

	// Left
	rlVertex3f(-1, -1, 1);
	rlVertex3f(-1, 1, -1);

	rlEnd();
	rlPopMatrix();
}

void MyDrawBox(Box box, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY)
{
	if (drawPolygon) MyDrawPolygonBox(box, polygonColor);
	if (drawWireframe) MyDrawWireframeBox(box, wireframeColor);
}

#pragma endregion

#pragma region disk

void MyDrawPolygonDisk(Disk disk, int nSectors, Color color = LIGHTGRAY)
{
	int numVertex = nSectors * 3;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = 2 * PI / nSectors;

	rlPushMatrix();
	rlTranslatef(disk.ref.origin.x, disk.ref.origin.y, disk.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(disk.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(disk.radius, 1, disk.radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	Cylindrical cylPoint;
	cylPoint.rho = 1;

	for (int i = 0; i < nSectors; i++)
	{
		cylPoint.theta = sectorAngle * i;
		cylPoint.y = 0;

		Vector3 point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		Vector3 point2 = CylindricalToCartesian(cylPoint);

		rlVertex3f(0, 0, 0);
		rlVertex3f(point1.x, point1.y, point1.z);
		rlVertex3f(point2.x, point2.y, point2.z);
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeDisk(Disk disk, int nSectors, Color color = DARKGRAY)
{
	int numVertex = nSectors * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = 2 * PI / nSectors;

	rlPushMatrix();
	rlTranslatef(disk.ref.origin.x, disk.ref.origin.y, disk.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(disk.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(disk.radius, 1, disk.radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	Cylindrical cylPoint;
	cylPoint.rho = 1;

	for (int i = 0; i < nSectors; i++)
	{
		cylPoint.theta = sectorAngle * i;
		cylPoint.y = 0;

		Vector3 point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		Vector3 point2 = CylindricalToCartesian(cylPoint);

		rlVertex3f(0, 0, 0);
		rlVertex3f(point1.x, point1.y, point1.z);

		rlVertex3f(point1.x, point1.y, point1.z);
		rlVertex3f(point2.x, point2.y, point2.z);

		rlVertex3f(point2.x, point2.y, point2.z);
		rlVertex3f(0, 0, 0);
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawDisk(Disk disk, int nSectors, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY)
{
	if (drawPolygon) MyDrawPolygonDisk(disk, nSectors, polygonColor);
	if (drawWireframe) MyDrawWireframeDisk(disk, nSectors, wireframeColor);
}

#pragma endregion

#pragma region sphere

void MyDrawPolygonSphere(Sphere sphere, int nMeridians, int nParallels, Color color = LIGHTGRAY)
{
	int numVertex = 6 * (nParallels)*nMeridians;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(sphere.ref.origin.x, sphere.ref.origin.y, sphere.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(sphere.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(sphere.radius, sphere.radius, sphere.radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float parallelAngle = (PI / nParallels);
	float meridianAngle = (2 * PI / nMeridians);

	Spherical sphPoint;
	Vector3 point;
	sphPoint.rho = 1;

	for (int j = 0; j < nParallels; j++)
	{
		for (int i = 0; i <= nMeridians; i++)
		{
			// Premier triangle
			sphPoint.theta = meridianAngle * i;
			sphPoint.phi = parallelAngle * j;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi += parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			sphPoint.phi -= parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			// Second triangle
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi -= parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeSphere(Sphere sphere, int nMeridians, int nParallels, Color color = DARKGRAY)
{
	int numVertex = 6 * (nParallels)*nMeridians;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(sphere.ref.origin.x, sphere.ref.origin.y, sphere.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(sphere.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(sphere.radius, sphere.radius, sphere.radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float parallelAngle = (PI / nParallels);
	float meridianAngle = (2 * PI / nMeridians);

	Spherical sphPoint;
	Vector3 point;
	sphPoint.rho = 1;

	for (int j = 0; j < nParallels; j++)
	{
		for (int i = 0; i <= nMeridians; i++)
		{
			sphPoint.theta = meridianAngle * i;
			sphPoint.phi = parallelAngle * j;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi += parallelAngle;
			sphPoint.theta -= meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);

			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi -= parallelAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawSphere(Sphere sphere, int nMeridians, int nParallels, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY)
{
	if (drawPolygon) MyDrawPolygonSphere(sphere, nMeridians, nParallels, polygonColor);
	if (drawWireframe) MyDrawWireframeSphere(sphere, nMeridians, nParallels, wireframeColor);
}

#pragma endregion

#pragma region cylinder

void MyDrawPolygonCylinder(Cylinder cylinder, int nSectors, bool drawCaps = false, Color color = LIGHTGRAY) {
	int numVertex = nSectors * 6;
	if (drawCaps) numVertex *= 2;

	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = 2 * PI / nSectors;

	rlPushMatrix();
	rlTranslatef(cylinder.ref.origin.x, cylinder.ref.origin.y, cylinder.ref.origin.z);
	rlScalef(cylinder.radius, cylinder.halfHeight, cylinder.radius);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(cylinder.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	for (int i = 0; i < nSectors; i++)
	{
		Cylindrical cylPoint;
		cylPoint.rho = cylinder.radius;
		cylPoint.theta = sectorAngle * i;
		cylPoint.y = -1;

		Vector3 point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		Vector3 point2 = CylindricalToCartesian(cylPoint);

		// Premier triangle
		rlVertex3f(point2.x, 1, point2.z);
		rlVertex3f(point1.x, 1, point1.z);
		rlVertex3f(point1.x, -1, point1.z);

		// Second triangle
		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point2.x, -1, point2.z);
		rlVertex3f(point2.x, 1, point2.z);

		// Fermeture du cylindre
		if (drawCaps) {
			rlVertex3f(0, 1, 0);
			rlVertex3f(point1.x, 1, point1.z);
			rlVertex3f(point2.x, 1, point2.z);

			rlVertex3f(0, -1, 0);
			rlVertex3f(point2.x, -1, point2.z);
			rlVertex3f(point1.x, -1, point1.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeCylinder(Cylinder cylinder, int nSectors, bool drawCaps = false, Color color = LIGHTGRAY) {
	int numVertex = nSectors * 6;
	if (drawCaps) numVertex *= 2;

	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = 2 * PI / nSectors;

	rlPushMatrix();
	rlTranslatef(cylinder.ref.origin.x, cylinder.ref.origin.y, cylinder.ref.origin.z);
	rlScalef(cylinder.radius, cylinder.halfHeight, cylinder.radius);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(cylinder.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	for (int i = 0; i < nSectors; i++)
	{
		Cylindrical cylPoint;
		cylPoint.rho = cylinder.radius;
		cylPoint.theta = sectorAngle * i;
		cylPoint.y = -1;

		Vector3 point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		Vector3 point2 = CylindricalToCartesian(cylPoint);

		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point2.x, -1, point2.z);

		rlVertex3f(point1.x, 1, point1.z);
		rlVertex3f(point2.x, 1, point2.z);

		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point1.x, 1, point1.z);

		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point2.x, 1, point2.z);

		if (drawCaps) {
			rlVertex3f(0, 1, 0);
			rlVertex3f(point1.x, 1, point1.z);

			rlVertex3f(0, 1, 0);
			rlVertex3f(point2.x, 1, point2.z);

			rlVertex3f(0, -1, 0);
			rlVertex3f(point1.x, -1, point1.z);

			rlVertex3f(0, -1, 0);
			rlVertex3f(point2.x, -1, point2.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawCylinder(Cylinder cylinder, int nSectors, bool drawCaps = false, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY) {
	if (drawPolygon) {
		MyDrawPolygonCylinder(cylinder, nSectors, drawCaps, polygonColor);

	}
	if (drawWireframe) {
		MyDrawWireframeCylinder(cylinder, nSectors, drawCaps, wireframeColor);
	}
}

#pragma endregion

#pragma region capsule


void MyDrawPolygonCapsule(Capsule capsule, int nSectors, int nParallels, Color color = LIGHTGRAY) {
	int numVertex = nSectors * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = 2 * PI / nSectors;

	rlPushMatrix();
	rlTranslatef(capsule.ref.origin.x, capsule.ref.origin.y, capsule.ref.origin.z);
	rlScalef(capsule.radius, capsule.halfHeight + capsule.radius, capsule.radius);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(capsule.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	for (int i = 0; i < nSectors; i++)
	{
		Cylindrical cylPoint;
		cylPoint.rho = capsule.radius;
		cylPoint.theta = sectorAngle * i;
		cylPoint.y = -1;

		Vector3 point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		Vector3 point2 = CylindricalToCartesian(cylPoint);

		// Premier triangle
		rlVertex3f(point2.x, 1, point2.z);
		rlVertex3f(point1.x, 1, point1.z);
		rlVertex3f(point1.x, -1, point1.z);

		// Second triangle
		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point2.x, -1, point2.z);
		rlVertex3f(point2.x, 1, point2.z);
	}

	float parallelAngle = (PI / nParallels);
	float meridianAngle = (2 * PI / nSectors);

	Vector3 point;
	Spherical sphTop = CartesianToSpherical({ 0,1,0 });
	Spherical sphBottom = CartesianToSpherical({ 0,-1,0 });
	sphTop.rho = 1;
	sphBottom.rho = 1;

	// Hemisphère Nord
	int j = 0;
	for (; j < nParallels / 2; j++) {
		for (int i = 0; i <= nSectors; i++)
		{
			// Premier triangle
			sphTop.theta = meridianAngle * i;
			sphTop.phi = parallelAngle * j;
			point = SphericalToCartesian(sphTop);
			rlVertex3f(point.x, point.y, point.z);
			sphTop.phi += parallelAngle;
			point = SphericalToCartesian(sphTop);
			rlVertex3f(point.x, point.y, point.z);
			sphTop.theta += meridianAngle;
			sphTop.phi -= parallelAngle;
			point = SphericalToCartesian(sphTop);
			rlVertex3f(point.x, point.y, point.z);

			// Second triangle
			rlVertex3f(point.x, point.y, point.z);
			sphTop.theta += meridianAngle;
			point = SphericalToCartesian(sphTop);
			rlVertex3f(point.x, point.y, point.z);
			sphTop.phi -= parallelAngle;
			point = SphericalToCartesian(sphTop);
			rlVertex3f(point.x, point.y, point.z);
		}
	}

	for (; j < nParallels; j++) {

	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeCapsule(Capsule capsule, int nSectors, int nParallels, Color color = LIGHTGRAY) {
	int numVertex = nSectors * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = 2 * PI / nSectors;

	rlPushMatrix();
	rlTranslatef(capsule.ref.origin.x, capsule.ref.origin.y, capsule.ref.origin.z);
	rlScalef(capsule.radius, capsule.halfHeight, capsule.radius);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(capsule.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	for (int i = 0; i < nSectors; i++)
	{
		Cylindrical cylPoint;
		cylPoint.rho = capsule.radius;
		cylPoint.theta = sectorAngle * i;
		cylPoint.y = -1;

		Vector3 point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		Vector3 point2 = CylindricalToCartesian(cylPoint);

		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point2.x, -1, point2.z);

		rlVertex3f(point1.x, 1, point1.z);
		rlVertex3f(point2.x, 1, point2.z);

		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point1.x, 1, point1.z);

		rlVertex3f(point1.x, -1, point1.z);
		rlVertex3f(point2.x, 1, point2.z);
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawCapsule(Capsule capsule, int nSectors, int nParallels, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY) {
	if (drawPolygon) {
		MyDrawPolygonCapsule(capsule, nSectors, nParallels, polygonColor);

	}
	if (drawWireframe) {
		MyDrawWireframeCapsule(capsule, nSectors, nParallels, wireframeColor);
	}
}

#pragma endregion

#pragma endregion

int main(int argc, char* argv[])
{
	// ICI LES TEST 

	// Initialization
	//--------------------------------------------------------------------------------------
	float screenSizeCoef = .9f;
	const int screenWidth = 1920 * screenSizeCoef;
	const int screenHeight = 1080 * screenSizeCoef;

	InitWindow(screenWidth, screenHeight, "ESIEE - E3FI - 2022 - 2023 - Maths 3D");

	SetTargetFPS(60);

	// CAMERA 
	Vector3 cameraPos = { 8, 15, 14 };
	Camera camera = { 0 };
	camera.position = cameraPos;
	camera.target = { 0, 0, 0 };
	camera.up = { 0, 1, 0 };
	camera.fovy = 120;
	camera.type = CAMERA_PERSPECTIVE;
	SetCameraMode(camera, CAMERA_CUSTOM);  // Set an orbital camera mode

	//--------------------------------------------------------------------------------------

	// Main game loop
	while (!WindowShouldClose())    // Detect window close button or ESC key
	{
		// Update
		//----------------------------------------------------------------------------------
		// TODO: Update your variables here
		//----------------------------------------------------------------------------------

		float deltaTime = GetFrameTime();
		float time = (float)GetTime();

		MyUpdateOrbitalCamera(&camera, deltaTime);

		// Draw
		//----------------------------------------------------------------------------------
		BeginDrawing();

		ClearBackground(RAYWHITE);

		BeginMode3D(camera);
		{
			//3D REFERENTIAL
			DrawGrid(20, 1);        // Draw a grid
			DrawLine3D({ 0 }, { 0,10,0 }, DARKGRAY);
			DrawSphere({ 10,0,0 }, .2f, RED);
			DrawSphere({ 0,10,0 }, .2f, GREEN);
			DrawSphere({ 0,0,10 }, .2f, BLUE);

			ReferenceFrame ref;

			// DISK DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ 0,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI/2));
			//Disk disk = { ref, 5 };
			//MyDrawDisk(disk, 15, true, true);

			//BOX DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ 0,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize({ 1,1,1 }), 0));
			//Box box = { ref, {5,3,2} };
			//MyDrawBox(box, true, true, BLUE, BLACK);

			//QUAD DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ 0,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize({ 1,1,1 }), PI/3));
			//Quad quad = { ref,{1,535,3} }; // x 
			//MyDrawQuad(quad, true, false);

			// SPHERE DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ 10,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), PI/2));
			//Sphere sphere = { ref, 5 };
			//MyDrawSphere(sphere, 15, 15, true, true);

			// CYLINDER DISPLAY TEST
			//ref = ReferenceFrame(
				// { 0,0,0 },
				// QuaternionFromAxisAngle(Vector3Normalize({ 1,1,1 }), 0));
			//Cylinder cylinder = { ref, 5, 3 };
			//MyDrawCylinder(cylinder, 30, true, true, true);

			// CAPSULE DISPLAY TEST
			//ref = ReferenceFrame(
			//	 { 0,0,0 },
			//	 QuaternionFromAxisAngle(Vector3Normalize({ 1,1,1 }), 0));
			//Capsule capsule = { ref, 10, 2 };
			//MyDrawCapsule(capsule, 15, 15, true, true);

		}
		EndMode3D();

		EndDrawing();
		//----------------------------------------------------------------------------------
	}

	// De-Initialization
	//--------------------------------------------------------------------------------------   
	CloseWindow();        // Close window and OpenGL context
	//--------------------------------------------------------------------------------------

	return 0;
}