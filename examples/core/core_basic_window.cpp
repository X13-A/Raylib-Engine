#pragma region init

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
#define GRAVITY {0, -9.81, 0}

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

#pragma endregion

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

struct Disk {
	ReferenceFrame ref;
	float radius;
};

struct Sphere {
	ReferenceFrame ref;
	float radius;
	Color color = SKYBLUE;
};

struct SphereCorner {
	ReferenceFrame ref;
	float radius;
};

struct Hemisphere {
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

struct CylinderQuarter {
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
	Color color = LIGHTGRAY;
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

struct Segment {
	Vector3 a;
	Vector3 b;
};

struct Line {
	Vector3 pt;
	Vector3 dir;
};

struct Plane {
	Vector3 normal;
	float d;
};

#pragma endregion

#pragma region Utils

Vector3 LocalToGlobalVect(Vector3 localVect, ReferenceFrame localRef)
{
	localVect = Vector3RotateByQuaternion(localVect, localRef.q);
	return Vector3Add(localVect, localRef.origin);
}

Vector3 GlobalToLocalVect(Vector3 globalVect, ReferenceFrame localRef)
{
	Vector3 localPos = Vector3Subtract(globalVect, localRef.origin);
	return Vector3RotateByQuaternion(localPos, QuaternionInvert(localRef.q));
}

Vector3 LocalToGlobalPos(Vector3 localPos, ReferenceFrame localRef)
{
	localPos = Vector3RotateByQuaternion(localPos, localRef.q);
	return Vector3Add(localPos, localRef.origin);
}

Vector3 GlobalToLocalPos(Vector3 globalPos, ReferenceFrame localRef)
{
	Vector3 localPos = Vector3Subtract(globalPos, localRef.origin);
	return Vector3RotateByQuaternion(localPos, QuaternionInvert(localRef.q));
}

Vector3 ProjectedPointOnLine(Vector3 linePt, Vector3 lineUnitDir, Vector3 pt) {
	Vector3 v = Vector3Subtract(pt, linePt);  // Vecteur allant de linePt à pt
	float d = Vector3DotProduct(v, lineUnitDir);  // Distance de pt au point projeté
	return Vector3Add(linePt, Vector3Scale(lineUnitDir, d));
}

bool IsPointInsideBox(Box box, Vector3 globalPt)
{
	Vector3 localPt = GlobalToLocalPos(globalPt, box.ref);
	bool w = localPt.x >= -box.extents.x && localPt.x <= box.extents.x;
	bool h = localPt.y >= -box.extents.y && localPt.y <= box.extents.y;
	bool d = localPt.z >= -box.extents.z && localPt.z <= box.extents.z;
	return (w && h && d);
}

#pragma endregion

#pragma region Conversion

Polar CartesianToPolar(Vector2 cart, bool keepThetaPositive = true)
{
	Polar polar = { sqrt(cart.x * cart.x + cart.y * cart.y), atan2f(cart.y,cart.x) };
	if (keepThetaPositive && polar.theta < 0)
	{
		polar.theta += 2 * PI;
	}
	return polar;
}

Vector2 PolarToCartesian(Polar polar)
{
	return Vector2Scale({ cosf(polar.theta),sinf(polar.theta) }, polar.rho);
}

Cylindrical CartesianToCylindrical(Vector3 cart, bool keepThetaPositive = true)
{
	Cylindrical cylindrical = { sqrt(pow(cart.x, 2) + pow(cart.z, 2)), atan2f(cart.x,cart.z), cart.y };
	if (keepThetaPositive && cylindrical.theta < 0)
	{
		cylindrical.theta += 2 * PI;
	}

	return cylindrical;
}

Vector3 CylindricalToCartesian(Cylindrical cylindrical)
{
	return Vector3({ cylindrical.rho * sinf(cylindrical.theta), cylindrical.y, cylindrical.rho * cosf(cylindrical.theta) });
}

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

Vector3 SphericalToCartesian(Spherical spherical)
{
	return Vector3({ spherical.rho * sinf(spherical.phi) * sinf(spherical.theta), spherical.rho * cosf(spherical.phi), spherical.rho * sinf(spherical.phi) * cosf(spherical.theta) });
}

#pragma endregion

#pragma region Camera

void MyUpdateOrbitalCamera(Camera* camera, float deltaTime)
{
	// Initialisation
	static Spherical sphPos = { 10,PI / 4.f,PI / 4.f };
	Spherical sphSpeed = { 2 ,0.04f,0.04f };
	float rhoMin = 4;
	float rhoMax = 40;
	Vector2 mousePos;
	Spherical sphDelta;

	// Calcul du vecteur de déplacement de la souris par différence entre la position courante de la souris et la position précédente
	static Vector2 prevMousePos = { 0, 0 };
	mousePos = GetMousePosition(); // Récupère la position de la souris
	Vector2 mouseVect = Vector2Subtract(mousePos, prevMousePos); // Différence entre la position actuelle et précédente de la souris
	prevMousePos = GetMousePosition(); // Mise à jour de la position précédente

	// Conversion de la position en coordonnées sphériques
	Spherical cam = CartesianToSpherical(camera->position);
	cam.rho -= GetMouseWheelMove(); // Zoom in / out

	// Calcul de la nouvelle position de la caméra en coordonnées sphériques
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

void MyDrawQuad(Quad quad, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY)
{
	if (drawPolygon) MyDrawPolygonQuad(quad, polygonColor);
	if (drawWireframe)MyDrawWireframeQuad(quad, wireframeColor);
}

#pragma endregion

#pragma region box

void MyDrawPolygonBox(Box box, Color color = LIGHTGRAY)
{
	// Vérifie si le tampon de rendu a suffisamment d'espace pour dessiner
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
	int numVertex = 6 * nParallels * nMeridians;
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
		for (int i = 0; i < nMeridians; i++)
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
	int numVertex = 6 * nParallels * nMeridians;
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
		for (int i = 0; i < nMeridians; i++)
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

#pragma region Sphere corner

void MyDrawPolygonSphereCorner(SphereCorner sphereCorner, int nMeridians, int nParallels, Color color = LIGHTGRAY)
{
	int numVertex = 6 * nParallels * nMeridians;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(sphereCorner.ref.origin.x, sphereCorner.ref.origin.y, sphereCorner.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(sphereCorner.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(sphereCorner.radius, sphereCorner.radius, sphereCorner.radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float parallelAngle = ((PI / 2) / nParallels);
	float meridianAngle = (2 * (PI / 4) / nMeridians);

	Spherical sphPoint;
	Vector3 point;
	sphPoint.rho = 1;

	for (int j = 0; j < nParallels; j++)
	{
		for (int i = 0; i < nMeridians; i++)
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

			//// Second triangle
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.phi += parallelAngle;
			sphPoint.theta -= meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeSphereCorner(SphereCorner sphereCorner, int nMeridians, int nParallels, Color color = DARKGRAY)
{
	int numVertex = (6 * nParallels * nMeridians) + (2 * nMeridians) + (2 * nParallels);
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(sphereCorner.ref.origin.x, sphereCorner.ref.origin.y, sphereCorner.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(sphereCorner.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(sphereCorner.radius, sphereCorner.radius, sphereCorner.radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float parallelAngle = ((PI / 2) / nParallels);
	float meridianAngle = (2 * (PI / 4) / nMeridians);

	Spherical sphPoint;
	Vector3 point;
	sphPoint.rho = 1;

	for (int j = 0; j < nParallels; j++)
	{
		for (int i = 0; i < nMeridians; i++)
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

			if (i == nMeridians - 1) {
				sphPoint.theta = meridianAngle * i;
				sphPoint.phi = parallelAngle * j;
				sphPoint.theta += meridianAngle;
				point = SphericalToCartesian(sphPoint);
				rlVertex3f(point.x, point.y, point.z);

				sphPoint.phi += parallelAngle;
				point = SphericalToCartesian(sphPoint);
				rlVertex3f(point.x, point.y, point.z);
			}

			if (j == nParallels - 1) {
				sphPoint.theta = meridianAngle * i;
				sphPoint.phi = parallelAngle * j;
				sphPoint.phi += parallelAngle;
				point = SphericalToCartesian(sphPoint);
				rlVertex3f(point.x, point.y, point.z);

				sphPoint.theta += meridianAngle;
				point = SphericalToCartesian(sphPoint);
				rlVertex3f(point.x, point.y, point.z);
			}
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawSphereCorner(SphereCorner sphereCorner, int nMeridians, int nParallels, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY)
{
	if (drawPolygon) MyDrawPolygonSphereCorner(sphereCorner, nMeridians, nParallels, polygonColor);
	if (drawWireframe) MyDrawWireframeSphereCorner(sphereCorner, nMeridians, nParallels, wireframeColor);
}

#pragma endregion

#pragma region Hemisphere

void MyDrawWireframeHemisphere(Hemisphere hemisphere, int nMeridians, int nParallels, Color color = DARKGRAY)
{
	int numVertex = (6 * nParallels * nMeridians) + (2 * nMeridians);
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(hemisphere.ref.origin.x, hemisphere.ref.origin.y, hemisphere.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(hemisphere.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(hemisphere.radius, hemisphere.radius, hemisphere.radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float parallelAngle = (PI / 2) / nParallels;
	float meridianAngle = (2 * PI) / nMeridians;

	Spherical sphPoint;
	Vector3 point;
	sphPoint.rho = 1;

	for (int j = 0; j < nParallels; j++)
	{
		for (int i = 0; i < nMeridians; i++)
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

			if (j == nParallels - 1) {
				sphPoint.theta = meridianAngle * i;
				sphPoint.phi = parallelAngle * j;
				sphPoint.phi += parallelAngle;
				point = SphericalToCartesian(sphPoint);
				rlVertex3f(point.x, point.y, point.z);

				sphPoint.theta += meridianAngle;
				point = SphericalToCartesian(sphPoint);
				rlVertex3f(point.x, point.y, point.z);
			}
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawPolygonHemisphere(Hemisphere hemisphere, int nMeridians, int nParallels, Color color = DARKGRAY)
{
	int numVertex = 6 * nParallels * nMeridians;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(hemisphere.ref.origin.x, hemisphere.ref.origin.y, hemisphere.ref.origin.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(hemisphere.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(hemisphere.radius, hemisphere.radius, hemisphere.radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float parallelAngle = (PI / 2) / nParallels;
	float meridianAngle = (2 * PI) / nMeridians;

	Spherical sphPoint;
	Vector3 point;
	sphPoint.rho = 1;

	for (int j = 0; j < nParallels; j++)
	{
		for (int i = 0; i < nMeridians; i++)
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
			sphPoint.phi += parallelAngle;
			sphPoint.theta -= meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
			sphPoint.theta += meridianAngle;
			point = SphericalToCartesian(sphPoint);
			rlVertex3f(point.x, point.y, point.z);
		}
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawHemisphere(Hemisphere hemisphere, int nMeridians, int nParallels, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY)
{
	if (drawPolygon) MyDrawPolygonHemisphere(hemisphere, nMeridians, nParallels, polygonColor);
	if (drawWireframe) MyDrawWireframeHemisphere(hemisphere, nMeridians, nParallels, wireframeColor);
}

#pragma endregion

#pragma region Cylinder

void MyDrawPolygonCylinder(Cylinder cylinder, int nSectors, bool drawCaps = false, Color color = LIGHTGRAY) {
	int numVertex = nSectors * 6;
	if (drawCaps) numVertex *= 2;

	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = 2 * PI / nSectors;

	rlPushMatrix();
	rlTranslatef(cylinder.ref.origin.x, cylinder.ref.origin.y, cylinder.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(cylinder.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(cylinder.radius, cylinder.halfHeight, cylinder.radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	Cylindrical cylPoint;
	cylPoint.rho = 1;
	Vector3 point1;
	Vector3 point2;

	for (int i = 0; i < nSectors; i++)
	{
		cylPoint.theta = sectorAngle * i;
		cylPoint.y = -1;

		point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		point2 = CylindricalToCartesian(cylPoint);

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

void MyDrawWireframeCylinder(Cylinder cylinder, int nSectors, bool drawCaps = false, Color color = DARKGRAY) {
	int numVertex = nSectors * 8;
	if (drawCaps) numVertex *= 2;

	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = 2 * PI / nSectors;

	rlPushMatrix();
	rlTranslatef(cylinder.ref.origin.x, cylinder.ref.origin.y, cylinder.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(cylinder.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(cylinder.radius, cylinder.halfHeight, cylinder.radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	Cylindrical cylPoint;
	cylPoint.rho = 1;
	Vector3 point1;
	Vector3 point2;

	for (int i = 0; i < nSectors; i++)
	{
		cylPoint.theta = sectorAngle * i;
		cylPoint.y = -1;

		point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		point2 = CylindricalToCartesian(cylPoint);

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

void MyDrawCylinder(Cylinder cylinder, int nSectors, bool drawCaps = false, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY)
{
	if (drawPolygon) MyDrawPolygonCylinder(cylinder, nSectors, drawCaps, polygonColor);
	if (drawWireframe) MyDrawWireframeCylinder(cylinder, nSectors, drawCaps, wireframeColor);
}

#pragma endregion

#pragma region cylinderCorner

void MyDrawPolygonCylinderQuarter(CylinderQuarter cylinderQuarter, int nSectors, bool drawCaps = false, Color color = LIGHTGRAY) {
	int numVertex = nSectors * 6;
	if (drawCaps) numVertex *= 2;

	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = (PI / 2) / nSectors;

	rlPushMatrix();
	rlTranslatef(cylinderQuarter.ref.origin.x, cylinderQuarter.ref.origin.y, cylinderQuarter.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(cylinderQuarter.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(cylinderQuarter.radius, cylinderQuarter.halfHeight, cylinderQuarter.radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	Cylindrical cylPoint;
	cylPoint.rho = 1;
	cylPoint.y = -1;

	Vector3 point1;
	Vector3 point2;

	for (int i = 0; i < nSectors; i++)
	{
		cylPoint.theta = sectorAngle * i;

		point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		point2 = CylindricalToCartesian(cylPoint);

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

void MyDrawWireframeCylinderQuarter(CylinderQuarter cylinderQuarter, int nSectors, bool drawCaps = false, Color color = DARKGRAY) {
	int numVertex = nSectors * 8;
	if (drawCaps) numVertex *= 2;

	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	float sectorAngle = (PI / 2) / nSectors;

	rlPushMatrix();
	rlTranslatef(cylinderQuarter.ref.origin.x, cylinderQuarter.ref.origin.y, cylinderQuarter.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(cylinderQuarter.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(cylinderQuarter.radius, cylinderQuarter.halfHeight, cylinderQuarter.radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	Cylindrical cylPoint;
	cylPoint.rho = 1;
	cylPoint.y = -1;
	Vector3 point1;
	Vector3 point2;

	for (int i = 0; i < nSectors; i++)
	{
		cylPoint.theta = sectorAngle * i;

		point1 = CylindricalToCartesian(cylPoint);
		cylPoint.theta = sectorAngle * (i + 1);
		point2 = CylindricalToCartesian(cylPoint);

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

void MyDrawCylinderQuarter(CylinderQuarter cylinderQuarter, int nSectors, bool drawCaps = false, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY)
{
	if (drawPolygon) MyDrawPolygonCylinderQuarter(cylinderQuarter, nSectors, drawCaps, polygonColor);
	if (drawWireframe) MyDrawWireframeCylinderQuarter(cylinderQuarter, nSectors, drawCaps, wireframeColor);
}

#pragma endregion

#pragma region Capsule

void MyDrawPolygonCapsule(Capsule capsule, int nSectors, int nParallels, Color color = LIGHTGRAY) {
	int numVertexCylinder = nSectors * 6;
	int numVertexHemispheres = 2 * (6 * nParallels * nSectors);
	int numVertex = numVertexCylinder + numVertexHemispheres;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	rlTranslatef(capsule.ref.origin.x, capsule.ref.origin.y, capsule.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(capsule.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	ReferenceFrame ref = ReferenceFrame(
		{ 0, 0, 0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	Cylinder cylinder = { ref, capsule.halfHeight, capsule.radius };
	MyDrawPolygonCylinder(cylinder, nSectors, false, color);

	ref = ReferenceFrame(
		{ 0, capsule.halfHeight, 0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), PI / 2));
	Hemisphere hemisphere = { ref, capsule.radius };
	MyDrawPolygonHemisphere(hemisphere, nSectors, nParallels, color);

	rlRotatef(PI * RAD2DEG, 1, 0, 0);
	ref = ReferenceFrame(
		{ 0, capsule.halfHeight, 0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	hemisphere = { ref, capsule.radius };
	MyDrawPolygonHemisphere(hemisphere, nSectors, nParallels, color);

	rlPopMatrix();
}

void MyDrawWireframeCapsule(Capsule capsule, int nSectors, int nParallels, Color color = DARKGRAY) {
	int numVertexCylinder = nSectors * 6;
	int numVertexHemispheres = 2 * ((6 * nParallels * nSectors) + (2 * nSectors));
	int numVertex = numVertexCylinder + numVertexHemispheres;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	rlTranslatef(capsule.ref.origin.x, capsule.ref.origin.y, capsule.ref.origin.z);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(capsule.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	ReferenceFrame ref = ReferenceFrame(
		{ 0, 0, 0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	Cylinder cylinder = { ref, capsule.halfHeight, capsule.radius };
	MyDrawWireframeCylinder(cylinder, nSectors, false, color);

	ref = ReferenceFrame(
		{ 0, capsule.halfHeight, 0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), PI / 2));
	Hemisphere hemisphere = { ref, capsule.radius };
	MyDrawWireframeHemisphere(hemisphere, nSectors, nParallels, color);

	rlRotatef(PI * RAD2DEG, 1, 0, 0);
	ref = ReferenceFrame(
		{ 0, capsule.halfHeight, 0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	hemisphere = { ref, capsule.radius };
	MyDrawWireframeHemisphere(hemisphere, nSectors, nParallels, color);
	rlRotatef(PI * RAD2DEG, 1, 0, 0);

	rlPopMatrix();
}

void MyDrawCapsule(Capsule capsule, int nSectors, int nParallels, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY)
{
	if (drawPolygon) MyDrawPolygonCapsule(capsule, nSectors, nParallels, polygonColor);
	if (drawWireframe) MyDrawWireframeCapsule(capsule, nSectors, nParallels, wireframeColor);
}

#pragma endregion

#pragma region RoundedBox

void MyDrawPolygonRoundedBox(RoundedBox roundedBox, int nSectors, Color color = LIGHTGRAY)
{
	int numVertex = 6 * 6; // Quads
	numVertex += 12 * nSectors * 6; // Cylinder quarters
	numVertex += 8 * 6 * nSectors * nSectors; // Sphere corners

	if (rlCheckBufferLimit(numVertex)) rlglDraw();
	rlPushMatrix();
	rlTranslatef(roundedBox.ref.origin.x, roundedBox.ref.origin.y, roundedBox.ref.origin.z);

	Vector3 vect;
	float angle;

	QuaternionToAxisAngle(roundedBox.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlColor4ub(color.r, color.g, color.b, color.a);

	#pragma region Pavé

	ReferenceFrame top_ref = ReferenceFrame(
		{ 0,roundedBox.extents.y + roundedBox.radius,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	Quad top_quad = { top_ref,{roundedBox.extents.x,0,roundedBox.extents.z} };
	MyDrawPolygonQuad(top_quad, color);

	ReferenceFrame bottom_ref = ReferenceFrame(
		{ 0,-roundedBox.extents.y - roundedBox.radius,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI));
	Quad bottom_quad = { bottom_ref,{roundedBox.extents.x,0,roundedBox.extents.z} };
	MyDrawPolygonQuad(bottom_quad, color);

	ReferenceFrame front_ref = ReferenceFrame(
		{ 0,0,roundedBox.extents.z + roundedBox.radius },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	Quad front_quad = { front_ref,{roundedBox.extents.x,0,roundedBox.extents.y} };
	MyDrawPolygonQuad(front_quad, color);

	ReferenceFrame back_ref = ReferenceFrame(
		{ 0,0,-roundedBox.extents.z - roundedBox.radius },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2));
	Quad back_quad = { back_ref,{roundedBox.extents.x,0,roundedBox.extents.y} };
	MyDrawPolygonQuad(back_quad, color);

	ReferenceFrame left_ref = ReferenceFrame(
		{ -roundedBox.extents.x - roundedBox.radius,0,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2)));
	Quad left_quad = { left_ref,{roundedBox.extents.z,0,roundedBox.extents.y} };
	MyDrawPolygonQuad(left_quad, color);

	ReferenceFrame right_ref = ReferenceFrame(
		{ roundedBox.extents.x + roundedBox.radius,0,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2)));
	Quad right_quad = { right_ref,{roundedBox.extents.z,0,roundedBox.extents.y} };
	MyDrawPolygonQuad(right_quad, color);

	#pragma endregion

	#pragma region Quart de cylindre vertical

	ReferenceFrame first_vertical_ref(
		{ roundedBox.extents.x,0,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	CylinderQuarter firstVerticalCylinderQuarter = { first_vertical_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(firstVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_vertical_ref(
		{ -roundedBox.extents.x,0,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), 3 * PI / 2));
	CylinderQuarter secondVerticalCylinderQuarter = { second_vertical_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(secondVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_vertical_ref(
		{ -roundedBox.extents.x,0,-roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	CylinderQuarter thirdVerticalCylinderQuarter = { third_vertical_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(thirdVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_vertical_ref(
		{ roundedBox.extents.x,0,-roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2));
	CylinderQuarter fourthVerticalCylinderQuarter = { fourth_vertical_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(fourthVerticalCylinderQuarter, nSectors, false, color);

	#pragma endregion

	#pragma region Quart de cylindre horizontal face supérieure

	ReferenceFrame first_superior_horizontal_ref(
		{ 0,roundedBox.extents.y,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2));
	CylinderQuarter firstSuperiorHorizontalCylinderQuarter = { first_superior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(firstSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_superior_horizontal_ref(
		{ 0,roundedBox.extents.y,-roundedBox.extents.z },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2)));
	CylinderQuarter secondSuperiorHorizontalCylinderQuarter = { second_superior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(secondSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_superior_horizontal_ref(
		{ roundedBox.extents.x,roundedBox.extents.y,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2));
	CylinderQuarter thirdSuperiorHorizontalCylinderQuarter = { third_superior_horizontal_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(thirdSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_superior_horizontal_ref(
		{ -roundedBox.extents.x,roundedBox.extents.y,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2)));
	CylinderQuarter fourthSuperiorHorizontalCylinderQuarter = { fourth_superior_horizontal_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(fourthSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	#pragma endregion

	#pragma region Quart de cylindre horizontal face inférieur

	ReferenceFrame first_inferior_horizontal_ref(
		{ 0,-roundedBox.extents.y,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2));
	CylinderQuarter firstInferiorHorizontalCylinderQuarter = { first_inferior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(firstInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_inferior_horizontal_ref(
		{ 0,-roundedBox.extents.y,-roundedBox.extents.z },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2)));
	CylinderQuarter secondInferiorHorizontalCylinderQuarter = { second_inferior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(secondInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_inferior_horizontal_ref(
		{ roundedBox.extents.x,-roundedBox.extents.y,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	CylinderQuarter thirdInferiorHorizontalCylinderQuarter = { third_inferior_horizontal_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(thirdInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_inferior_horizontal_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.y,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2)));
	CylinderQuarter fourthInferiorHorizontalCylinderQuarter = { fourth_inferior_horizontal_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(fourthInferiorHorizontalCylinderQuarter, nSectors, false, color);

	#pragma endregion

	#pragma region Coin face supérieur

	ReferenceFrame first_superior_corner_ref(
		{ roundedBox.extents.x,roundedBox.extents.y,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	SphereCorner firstSuperiorCorner = { first_superior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(firstSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame second_superior_corner_ref(
		{ roundedBox.extents.x,roundedBox.extents.y,-roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2));
	SphereCorner secondSuperiorCorner = { second_superior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(secondSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame third_superior_corner_ref(
		{ -roundedBox.extents.x,roundedBox.extents.y,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2));
	SphereCorner thirdSuperiorCorner = { third_superior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(thirdSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame fourth_superior_corner_ref(
		{ -roundedBox.extents.x,roundedBox.extents.y,-roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	SphereCorner fourthSuperiorCorner = { fourth_superior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(fourthSuperiorCorner, nSectors, nSectors, color);

	#pragma endregion

	#pragma region Coin face inférieur

	ReferenceFrame first_inferior_corner_ref(
		{ roundedBox.extents.x,-roundedBox.extents.y,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	SphereCorner firstInferiorCorner = { first_inferior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(firstInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame second_inferior_corner_ref(
		{ roundedBox.extents.x,-roundedBox.extents.y,-roundedBox.extents.z },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2)));
	SphereCorner secondInferiorCorner = { second_inferior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(secondInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame third_inferior_corner_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.y,roundedBox.extents.z },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI)));
	SphereCorner thirdInferiorCorner = { third_inferior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(thirdInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame fourth_inferior_corner_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.y,-roundedBox.extents.z },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI),
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2)));
	SphereCorner fourthInferiorCorner = { fourth_inferior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(fourthInferiorCorner, nSectors, nSectors, color);

	#pragma endregion

	rlEnd();
	rlPopMatrix();
}

void MyDrawWireframeRoundedBox(RoundedBox roundedBox, int nSectors, Color color = DARKGRAY)
{
	int numVertex = 6 * 10; // Quads
	numVertex += 12 * nSectors * 8; // Cylinder quarters
	numVertex += 8 * ((6 * nSectors * nSectors) + (2 * nSectors) + (2 * nSectors)); // Sphere corners

	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	rlTranslatef(roundedBox.ref.origin.x, roundedBox.ref.origin.y, roundedBox.ref.origin.z);

	Vector3 vect;
	float angle;

	QuaternionToAxisAngle(roundedBox.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlColor4ub(color.r, color.g, color.b, color.a);

	#pragma region Pavé

	ReferenceFrame top_ref = ReferenceFrame(
		{ 0,roundedBox.extents.y + roundedBox.radius,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	Quad top_quad = { top_ref,{roundedBox.extents.x,0,roundedBox.extents.z} };
	MyDrawWireframeQuad(top_quad, color);

	ReferenceFrame bottom_ref = ReferenceFrame(
		{ 0,-roundedBox.extents.y - roundedBox.radius,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI));
	Quad bottom_quad = { bottom_ref,{roundedBox.extents.x,0,roundedBox.extents.z} };
	MyDrawWireframeQuad(bottom_quad, color);

	ReferenceFrame front_ref = ReferenceFrame(
		{ 0,0,roundedBox.extents.z + roundedBox.radius },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	Quad front_quad = { front_ref,{roundedBox.extents.x,0,roundedBox.extents.y} };
	MyDrawWireframeQuad(front_quad, color);

	ReferenceFrame back_ref = ReferenceFrame(
		{ 0,0,-roundedBox.extents.z - roundedBox.radius },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2));
	Quad back_quad = { back_ref,{roundedBox.extents.x,0,roundedBox.extents.y} };
	MyDrawWireframeQuad(back_quad, color);

	ReferenceFrame left_ref = ReferenceFrame(
		{ -roundedBox.extents.x - roundedBox.radius,0,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2)));
	Quad left_quad = { left_ref,{roundedBox.extents.z,0,roundedBox.extents.y} };
	MyDrawWireframeQuad(left_quad, color);

	ReferenceFrame right_ref = ReferenceFrame(
		{ roundedBox.extents.x + roundedBox.radius,0,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2)));
	Quad right_quad = { right_ref,{roundedBox.extents.z,0,roundedBox.extents.y} };
	MyDrawWireframeQuad(right_quad, color);

	#pragma endregion

	#pragma region Quart de cylindre vertical

	ReferenceFrame first_vertical_ref(
		{ roundedBox.extents.x,0,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	CylinderQuarter firstVerticalCylinderQuarter = { first_vertical_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(firstVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_vertical_ref(
		{ -roundedBox.extents.x,0,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), 3 * PI / 2));
	CylinderQuarter secondVerticalCylinderQuarter = { second_vertical_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(secondVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_vertical_ref(
		{ -roundedBox.extents.x,0,-roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	CylinderQuarter thirdVerticalCylinderQuarter = { third_vertical_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(thirdVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_vertical_ref(
		{ roundedBox.extents.x,0,-roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2));
	CylinderQuarter fourthVerticalCylinderQuarter = { fourth_vertical_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(fourthVerticalCylinderQuarter, nSectors, false, color);

	#pragma endregion

	#pragma region Quart de cylindre horizontal face supérieure

	ReferenceFrame first_superior_horizontal_ref(
		{ 0,roundedBox.extents.y,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2));
	CylinderQuarter firstSuperiorHorizontalCylinderQuarter = { first_superior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(firstSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_superior_horizontal_ref(
		{ 0,roundedBox.extents.y,-roundedBox.extents.z },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2)));
	CylinderQuarter secondSuperiorHorizontalCylinderQuarter = { second_superior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(secondSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_superior_horizontal_ref(
		{ roundedBox.extents.x,roundedBox.extents.y,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2));
	CylinderQuarter thirdSuperiorHorizontalCylinderQuarter = { third_superior_horizontal_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(thirdSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_superior_horizontal_ref(
		{ -roundedBox.extents.x,roundedBox.extents.y,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2)));
	CylinderQuarter fourthSuperiorHorizontalCylinderQuarter = { fourth_superior_horizontal_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(fourthSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	#pragma endregion

	#pragma region Quart de cylindre horizontal face inférieur

	ReferenceFrame first_inferior_horizontal_ref(
		{ 0,-roundedBox.extents.y,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2));
	CylinderQuarter firstInferiorHorizontalCylinderQuarter = { first_inferior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(firstInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_inferior_horizontal_ref(
		{ 0,-roundedBox.extents.y,-roundedBox.extents.z },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2)));
	CylinderQuarter secondInferiorHorizontalCylinderQuarter = { second_inferior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(secondInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_inferior_horizontal_ref(
		{ roundedBox.extents.x,-roundedBox.extents.y,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	CylinderQuarter thirdInferiorHorizontalCylinderQuarter = { third_inferior_horizontal_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(thirdInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_inferior_horizontal_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.y,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2)));
	CylinderQuarter fourthInferiorHorizontalCylinderQuarter = { fourth_inferior_horizontal_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(fourthInferiorHorizontalCylinderQuarter, nSectors, false, color);

	#pragma endregion

	#pragma region Coin face supérieur

	ReferenceFrame first_superior_corner_ref(
		{ roundedBox.extents.x,roundedBox.extents.y,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	SphereCorner firstSuperiorCorner = { first_superior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(firstSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame second_superior_corner_ref(
		{ roundedBox.extents.x,roundedBox.extents.y,-roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2));
	SphereCorner secondSuperiorCorner = { second_superior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(secondSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame third_superior_corner_ref(
		{ -roundedBox.extents.x,roundedBox.extents.y,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2));
	SphereCorner thirdSuperiorCorner = { third_superior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(thirdSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame fourth_superior_corner_ref(
		{ -roundedBox.extents.x,roundedBox.extents.y,-roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	SphereCorner fourthSuperiorCorner = { fourth_superior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(fourthSuperiorCorner, nSectors, nSectors, color);

	#pragma endregion

	#pragma region Coin face inférieur

	ReferenceFrame first_inferior_corner_ref(
		{ roundedBox.extents.x,-roundedBox.extents.y,roundedBox.extents.z },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	SphereCorner firstInferiorCorner = { first_inferior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(firstInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame second_inferior_corner_ref(
		{ roundedBox.extents.x,-roundedBox.extents.y,-roundedBox.extents.z },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2)));
	SphereCorner secondInferiorCorner = { second_inferior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(secondInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame third_inferior_corner_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.y,roundedBox.extents.z },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI)));
	SphereCorner thirdInferiorCorner = { third_inferior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(thirdInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame fourth_inferior_corner_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.y,-roundedBox.extents.z },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI),
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2)));
	SphereCorner fourthInferiorCorner = { fourth_inferior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(fourthInferiorCorner, nSectors, nSectors, color);

	#pragma endregion

	rlEnd();
	rlPopMatrix();
}

void MyDrawRoundedBox(RoundedBox roundedBox, int nSectors, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY)
{
	if (drawPolygon) MyDrawPolygonRoundedBox(roundedBox, nSectors, polygonColor);
	if (drawWireframe) MyDrawWireframeRoundedBox(roundedBox, nSectors, wireframeColor);
}

#pragma endregion

#pragma region segment

void MyDrawSegment(Segment segment, Color colorA = RED, Color colorB = GREEN, Color color = BLACK) {
	Vector3 a = segment.a;
	Vector3 b = segment.b;
	
	DrawLine3D(a, b, color);
	MyDrawPolygonSphere({ {a,QuaternionIdentity()},.15f }, 16, 8, colorA);
	MyDrawPolygonSphere({ {b,QuaternionIdentity()},.15f }, 16, 8, colorB);

}

#pragma endregion

#pragma endregion

#pragma region Intersections

bool IntersectLinePlane(Line line, Plane plane, float& t, Vector3& interPt, Vector3& interNormal)
{
	// no intersection if line is parallel to the plane
	float dotProd = Vector3DotProduct(plane.normal, line.dir);
	if (fabsf(dotProd) < EPSILON) return false;

	// intersection: t, interPt & interNormal
	t = (plane.d - Vector3DotProduct(plane.normal, line.pt)) / dotProd;
	interPt = Vector3Add(line.pt, Vector3Scale(line.dir, t)); // OM = OA+tAB
	interNormal = Vector3Scale(plane.normal, Vector3DotProduct(Vector3Subtract(line.pt, interPt), plane.normal) < 0 ? -1.f : 1.f);

	return true;
}

bool IntersectSegmentPlane(Segment segment, Plane plane, float& t, Vector3& interPt, Vector3& interNormal)
{
	// créer une droite à partir du segment
	Line line;
	line.pt = segment.a;
	line.dir = Vector3Subtract(segment.b, segment.a);

	// Vérifie s'il y a une intersection entre la droite et le plan
	if (!IntersectLinePlane(line, plane, t, interPt, interNormal)) return false;

	// Vérifie si le point d'intersection est sur le segment
	if (t < 0 || t > 1) return false;

	return true;
}

bool IntersectSegmentQuad(Segment segment, Quad quad, float& t, Vector3& interPt, Vector3& interNormal)
{
	// Coordonnées locales des points a et b
	Vector3 a = GlobalToLocalPos(segment.a, quad.ref);
	Vector3 b = GlobalToLocalPos(segment.b, quad.ref);

	// Si le segment ne traverse pas la plane représentant le quad, on retourne false
	if (!((a.y >= 0 && b.y <= 0) || (b.y >= 0 && a.y <= 0)))
	{
		return false;
	}
	if ((b.y - a.y) == 0) return false;

	// Calcul du t représentant le pourcentage du segment parcouru avant l'intersection sur le plan
	t = -a.y / (b.y - a.y);

	// Si l'intersection n'est pas sur le segment, on retourne false
	if (t < 0 || t > 1) return false;

	// Calcul des coordonnées du point d'intersection, si le point n'est pas compris dans les dimensions du quad, on retourne false
	interPt = Vector3Add(a, Vector3Scale(Vector3Subtract(b, a), t));
	if (fabs(interPt.x) > fabs(quad.extents.x) || fabs(interPt.z) > fabs(quad.extents.z)) return false;

	// Calcul des valeurs renvoyées
	interPt = LocalToGlobalPos(interPt, quad.ref);
	interNormal = Vector3RotateByQuaternion({ 0,1,0 }, quad.ref.q);

	return true;
}

bool IntersectSegmentDisk(Segment segment, Disk disk, float& t, Vector3& interPt, Vector3& interNormal)
{
	// Même principe que pour le quad
	Vector3 a = GlobalToLocalPos(segment.a, disk.ref);
	Vector3 b = GlobalToLocalPos(segment.b, disk.ref);
	if (!((a.y >= 0 && b.y <= 0) || (b.y >= 0 && a.y <= 0)))
	{
		return false;
	}
	t = -a.y / (b.y - a.y);

	// La seule différence est qu'on vérifie si le point d'intersection est dans le radius (plutôt que dans les extents)
	interPt = Vector3Add(a, Vector3Scale(Vector3Subtract(b, a), t));
	if (Vector3Distance(interPt, { 0,0,0 }) > disk.radius) return false;

	interPt = LocalToGlobalPos(interPt, disk.ref);
	interNormal = Vector3RotateByQuaternion({ 0,1,0 }, disk.ref.q);

	return true;
}

bool IntersectSegmentSphere(Segment seg, Sphere s, float& t, Vector3& interPt, Vector3& interNormal)
{
	// Origine du segment
	Vector3 Ro = seg.a;

	// Direction du segment
	Vector3 Rd = Vector3Subtract(seg.b, seg.a);

	// Origine de la sphère
	Vector3 S = s.ref.origin;

	// Point le plus proche de l'origine de la sphère qui est sur le segment
	Vector3 tp = ProjectedPointOnLine(Ro, Vector3Normalize(Rd), S);

	// Calcul de la distance entre tp et la surface de la sphère
	float r = s.radius;
	float y = Vector3Distance(S, tp);
	if (y > r) return false;
	float x = sqrt((r * r - y * y));

	// Calcul des pourcentages menant aux points d'intersections
	float t0 = Vector3Distance(Ro, tp) / Vector3Distance(seg.a, seg.b);
	float t1 = t0 - x / Vector3Distance(seg.a, seg.b);
	float t2 = t0 + x / Vector3Distance(seg.a, seg.b);

	// Choix du bon t
	char bestT = '0';
	if (t1 < 0 || t1 > 1) t1 = INFINITY;
	if (t2 < 0 || t2 > 1) t2 = INFINITY;

	if (t1 < t2 && t1 < 1) bestT = '1';
	if (t2 < t1 && t2 < 1) bestT = '2';
	if (bestT == '0') return false;

	// Calcul du point d'intersection
	interPt = bestT == '1' ? Vector3Add(Ro, Vector3Scale(Rd, t1)) : Vector3Add(Ro, Vector3Scale(Rd, t2));
	t = bestT == '1' ? t1 : t2;

	interNormal = Vector3Normalize(Vector3Subtract(interPt, S));
	return true;
}

bool IntersectSegmentBox(Segment seg, Box box, float& t, Vector3& interPt, Vector3& interNormal)
{
	// Front quad
	ReferenceFrame front_ref = box.ref;
	front_ref.Translate(Vector3RotateByQuaternion({ 0,0,box.extents.z }, front_ref.q));
	front_ref.q = QuaternionMultiply(front_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, PI / 2));
	Quad front_quad = { front_ref, {box.extents.x, 0, box.extents.y} };

	// Back quad
	ReferenceFrame back_ref = box.ref;
	back_ref.Translate(Vector3RotateByQuaternion({ 0,0,-box.extents.z }, back_ref.q));
	back_ref.q = QuaternionMultiply(back_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, -PI / 2));
	Quad back_quad = { back_ref, {box.extents.x, 0, box.extents.y} };

	// Left quad
	ReferenceFrame left_ref = box.ref;
	left_ref.Translate(Vector3RotateByQuaternion({ -box.extents.x,0,0 }, left_ref.q));
	left_ref.q = QuaternionMultiply(left_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI / 2));
	Quad left_quad = { left_ref, {box.extents.y, 0, box.extents.z} };

	// Right quad
	ReferenceFrame right_ref = box.ref;
	right_ref.Translate(Vector3RotateByQuaternion({ box.extents.x,0,0 }, right_ref.q));
	right_ref.q = QuaternionMultiply(right_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, -PI / 2));
	Quad right_quad = { right_ref, {box.extents.y, 0, box.extents.z} };

	// Top quad
	ReferenceFrame top_ref = box.ref;
	top_ref.Translate(Vector3RotateByQuaternion({ 0,box.extents.y,0 }, top_ref.q));
	top_ref.q = QuaternionMultiply(top_ref.q, QuaternionFromAxisAngle({ 0,0,0 }, 0));
	Quad top_quad = { top_ref, {box.extents.x, 0, box.extents.z} };

	// Bottom quad
	ReferenceFrame bottom_ref = box.ref;
	bottom_ref.Translate(Vector3RotateByQuaternion({ 0,-box.extents.y,0 }, bottom_ref.q));
	bottom_ref.q = QuaternionMultiply(bottom_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI));
	Quad bottom_quad = { bottom_ref, {box.extents.x, 0, box.extents.z} };

	Quad faces[6] = { front_quad, back_quad, left_quad, right_quad, top_quad, bottom_quad };

	float tempT;
	Vector3 tempInterPt;
	Vector3 tempInterNormal;
	float minDist = INFINITY;
	float dist;
	bool intersect = false;

	// Vérifie si il y intersection sur chaque face de la box
	for (int i = 0; i < 6; i++) {
		if (IntersectSegmentQuad(seg, faces[i], tempT, tempInterPt, tempInterNormal))
		{
			intersect = true;
			dist = Vector3Distance(seg.a, tempInterPt);
			if (dist < minDist)
			{
				minDist = dist;
				t = tempT;
				interPt = tempInterPt;
				interNormal = tempInterNormal;
			}
		}
	}

	return intersect;
}

bool IntersectSegmentInfiniteCylinder(Segment segment, Cylinder cylinder, float& t, Vector3& interPt, Vector3& interNormal)
{
	Vector3 aLocal = GlobalToLocalPos(segment.a, cylinder.ref);
	Vector3 bLocal = GlobalToLocalPos(segment.b, cylinder.ref);

	float p1x = aLocal.x;
	float p2x = bLocal.x;

	float p1z = aLocal.z;
	float p2z = bLocal.z;

	//Polynôme du second degré
	float a = pow((p2x - p1x), 2) + pow((p2z - p1z), 2);
	float b = 2 * p1x * (p2x - p1x) + 2 * p1z * (p2z - p1z);
	float c = pow(p1x, 2) + pow(p1z, 2) - pow(cylinder.radius, 2);
	float delta = pow(b, 2) - 4 * a * c;

	if (delta > 0)
	{
		float t1 = (-b + sqrtf(delta)) / (2 * a);
		float t2 = (-b - sqrtf(delta)) / (2 * a);

		// Choix du bon t
		char bestT = '0';
		if (t1 < 0 || t1 > 1) t1 = INFINITY;
		if (t2 < 0 || t2 > 1) t2 = INFINITY;

		if (t1 < t2 && t1 < 1) bestT = '1';
		if (t2 < t1 && t2 < 1) bestT = '2';
		if (bestT == '0') return false;

		t = bestT == '1' ? t1 : t2;

		interPt = Vector3Add(segment.a, Vector3Scale(Vector3Subtract(segment.b, segment.a), t));
		Vector3 interPtLocal = GlobalToLocalPos(interPt, cylinder.ref);
		interNormal = Vector3Subtract(interPtLocal, { 0, interPtLocal.y, 0 });
		interNormal = LocalToGlobalPos(Vector3Add(interPtLocal, interNormal), cylinder.ref);
		interNormal = Vector3Normalize(Vector3Subtract(interNormal, interPt));
		return true;
	}

	// Cas très rare (delta == 0)
	if (fabs(delta) <= EPSILON)
	{
		t = -b / (2 * a);
		return true;
	}

	return false;
}

bool IntersectSegmentCylinder(Segment segment, Cylinder cylinder, float& t, Vector3& interPt, Vector3& interNormal)
{
	// Calcul des "bouchons" du cylindre
	ReferenceFrame top_ref = cylinder.ref;
	top_ref.Translate(Vector3RotateByQuaternion({ 0, cylinder.halfHeight,0 }, top_ref.q));
	Disk top_disk = { top_ref, cylinder.radius };

	ReferenceFrame bottom_ref = cylinder.ref;
	bottom_ref.Translate(Vector3RotateByQuaternion({ 0, -cylinder.halfHeight,0 }, bottom_ref.q));
	bottom_ref.q = QuaternionMultiply(bottom_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI));
	Disk bottom_disk = { bottom_ref, cylinder.radius };

	Disk caps[2] = { top_disk, bottom_disk };

	float tempT;
	Vector3 tempInterPt;
	Vector3 tempInterNormal;
	float minDist = INFINITY;
	float dist;
	bool intersect = false;

	// Vérification sur les deux bouchons
	for (int i = 0; i < 2; i++) {
		if (IntersectSegmentDisk(segment, caps[i], tempT, tempInterPt, tempInterNormal))
		{
			intersect = true;
			dist = Vector3Distance(segment.a, tempInterPt);
			if (dist < minDist)
			{
				minDist = dist;
				t = tempT;
				interPt = tempInterPt;
				interNormal = tempInterNormal;
			}
		}
	}

	// Vérification sur le cylindre infini
	if (IntersectSegmentInfiniteCylinder(segment, cylinder, tempT, tempInterPt, tempInterNormal)) {
		Vector3 interPtLocal = GlobalToLocalPos(tempInterPt, cylinder.ref);
		if (fabs(interPtLocal.y) < fabs(cylinder.halfHeight))
		{
			intersect = true;
			dist = Vector3Distance(segment.a, tempInterPt);
			if (dist < minDist)
			{
				minDist = dist;
				t = tempT;
				interPt = tempInterPt;
				interNormal = tempInterNormal;
			}
		}
	}
	if (t < 0 || t > 1) return false;
	return intersect;
}

bool IntersectSegmentCapsule(Segment segment, Capsule capsule, float& t, Vector3& interPt, Vector3& interNormal)
{
	// Calcul des sphères aux extrémitées
	ReferenceFrame top_ref = capsule.ref;
	top_ref.Translate(Vector3RotateByQuaternion({ 0, capsule.halfHeight,0 }, top_ref.q));
	Sphere top_sphere = { top_ref, capsule.radius };

	ReferenceFrame bottom_ref = capsule.ref;
	bottom_ref.Translate(Vector3RotateByQuaternion({ 0, -capsule.halfHeight,0 }, bottom_ref.q));
	bottom_ref.q = QuaternionMultiply(bottom_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI));
	Sphere bottom_sphere = { bottom_ref, capsule.radius };

	Sphere spheres[2] = { top_sphere, bottom_sphere };

	float tempT;
	Vector3 tempInterPt;
	Vector3 tempInterNormal;
	float minDist = INFINITY;
	float dist;
	bool intersect = false;

	// Vérification sur les deux sphères
	for (int i = 0; i < 2; i++) {
		if (IntersectSegmentSphere(segment, spheres[i], tempT, tempInterPt, tempInterNormal))
		{
			intersect = true;
			dist = Vector3Distance(segment.a, tempInterPt);
			if (dist < minDist)
			{
				minDist = dist;
				t = tempT;
				interPt = tempInterPt;
				interNormal = tempInterNormal;
			}
		}
	}

	ReferenceFrame cylinder_ref = capsule.ref;
	Cylinder cylinder = { capsule.ref, capsule.halfHeight, capsule.radius };

	// Vérification sur le cylindre
	if (IntersectSegmentCylinder(segment, cylinder, tempT, tempInterPt, tempInterNormal)) {
		intersect = true;
		dist = Vector3Distance(segment.a, tempInterPt);
		if (dist < minDist)
		{
			minDist = dist;
			t = tempT;
			interPt = tempInterPt;
			interNormal = tempInterNormal;
		}
	}

	return intersect;
}

bool IntersectSegmentRoundedBox(Segment seg, RoundedBox rndBox, float& t, Vector3& interPt, Vector3& interNormal)
{
	#pragma region pavé

	// Front quad
	ReferenceFrame front_ref = rndBox.ref;
	front_ref.Translate(Vector3RotateByQuaternion({ 0,0,rndBox.extents.z + rndBox.radius }, front_ref.q));
	front_ref.q = QuaternionMultiply(front_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, PI / 2));
	Quad front_quad = { front_ref, {rndBox.extents.x, 0, rndBox.extents.y} };

	// Back quad
	ReferenceFrame back_ref = rndBox.ref;
	back_ref.Translate(Vector3RotateByQuaternion({ 0,0,-rndBox.extents.z - rndBox.radius }, back_ref.q));
	back_ref.q = QuaternionMultiply(back_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, -PI / 2));
	Quad back_quad = { back_ref, {rndBox.extents.x, 0, rndBox.extents.y} };

	// Left quad
	ReferenceFrame left_ref = rndBox.ref;
	left_ref.Translate(Vector3RotateByQuaternion({ -rndBox.extents.x - rndBox.radius,0,0 }, left_ref.q));
	left_ref.q = QuaternionMultiply(left_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI / 2));
	Quad left_quad = { left_ref, {rndBox.extents.y, 0, rndBox.extents.z} };

	// Right quad
	ReferenceFrame right_ref = rndBox.ref;
	right_ref.Translate(Vector3RotateByQuaternion({ rndBox.extents.x + rndBox.radius,0,0 }, right_ref.q));
	right_ref.q = QuaternionMultiply(right_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, -PI / 2));
	Quad right_quad = { right_ref, {rndBox.extents.y, 0, rndBox.extents.z} };

	// Top quad
	ReferenceFrame top_ref = rndBox.ref;
	top_ref.Translate(Vector3RotateByQuaternion({ 0,rndBox.extents.y + rndBox.radius,0 }, top_ref.q));
	top_ref.q = QuaternionMultiply(top_ref.q, QuaternionFromAxisAngle({ 0,0,0 }, 0));
	Quad top_quad = { top_ref, {rndBox.extents.x, 0, rndBox.extents.z} };

	// Bottom quad
	ReferenceFrame bottom_ref = rndBox.ref;
	bottom_ref.Translate(Vector3RotateByQuaternion({ 0,-rndBox.extents.y - rndBox.radius,0 }, bottom_ref.q));
	bottom_ref.q = QuaternionMultiply(bottom_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI));
	Quad bottom_quad = { bottom_ref, {rndBox.extents.x, 0, rndBox.extents.z} };

	#pragma endregion

	#pragma region vertical cylinders

	ReferenceFrame left_front_vertical_ref = rndBox.ref;
	left_front_vertical_ref.Translate(Vector3RotateByQuaternion({ -rndBox.extents.x ,0,rndBox.extents.z }, left_front_vertical_ref.q));
	left_front_vertical_ref.q = QuaternionMultiply(left_front_vertical_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, 3 * PI / 2));
	Cylinder left_front_vertical_cylinder = { left_front_vertical_ref, rndBox.extents.y, rndBox.radius };

	ReferenceFrame right_front_vertical_ref = rndBox.ref;
	right_front_vertical_ref.Translate(Vector3RotateByQuaternion({ rndBox.extents.x ,0,rndBox.extents.z }, right_front_vertical_ref.q));
	Cylinder right_front_vertical_cylinder = { right_front_vertical_ref, rndBox.extents.y, rndBox.radius };

	ReferenceFrame right_back_vertical_ref = rndBox.ref;
	right_back_vertical_ref.Translate(Vector3RotateByQuaternion({ rndBox.extents.x ,0,-rndBox.extents.z }, right_back_vertical_ref.q));
	right_back_vertical_ref.q = QuaternionMultiply(right_back_vertical_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, PI / 2));
	Cylinder right_back_vertical_cylinder = { right_back_vertical_ref, rndBox.extents.y, rndBox.radius };

	ReferenceFrame left_back_vertical_ref = rndBox.ref;
	left_back_vertical_ref.Translate(Vector3RotateByQuaternion({ -rndBox.extents.x ,0,-rndBox.extents.z }, left_back_vertical_ref.q));
	left_back_vertical_ref.q = QuaternionMultiply(left_back_vertical_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, PI));
	Cylinder left_back_vertical_cylinder = { left_back_vertical_ref, rndBox.extents.y, rndBox.radius };

	#pragma endregion

	#pragma region corners

	// Top
	ReferenceFrame top_front_left_corner_ref = rndBox.ref;
	top_front_left_corner_ref.Translate(Vector3RotateByQuaternion({ -rndBox.extents.x, rndBox.extents.y, rndBox.extents.z }, top_front_left_corner_ref.q));
	top_front_left_corner_ref.q = QuaternionMultiply(top_front_left_corner_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, -PI / 2));
	Sphere top_front_left_sphere = { top_front_left_corner_ref, rndBox.radius };

	ReferenceFrame top_front_right_corner_ref = rndBox.ref;
	top_front_right_corner_ref.Translate(Vector3RotateByQuaternion({ rndBox.extents.x, rndBox.extents.y, rndBox.extents.z }, top_front_right_corner_ref.q));
	Sphere top_front_right_sphere = { top_front_right_corner_ref, rndBox.radius };

	ReferenceFrame top_back_left_corner_ref = rndBox.ref;
	top_back_left_corner_ref.Translate(Vector3RotateByQuaternion({ -rndBox.extents.x, rndBox.extents.y, -rndBox.extents.z }, top_back_left_corner_ref.q));
	top_back_left_corner_ref.q = QuaternionMultiply(top_back_left_corner_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, PI));
	Sphere top_back_left_sphere = { top_back_left_corner_ref, rndBox.radius };

	ReferenceFrame top_back_right_corner_ref = rndBox.ref;
	top_back_right_corner_ref.Translate(Vector3RotateByQuaternion({ rndBox.extents.x, rndBox.extents.y, -rndBox.extents.z }, top_back_right_corner_ref.q));
	top_back_right_corner_ref.q = QuaternionMultiply(top_back_right_corner_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, PI / 2));
	top_back_right_corner_ref.q = QuaternionMultiply(top_back_right_corner_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, PI));
	Sphere top_back_right_sphere = { top_back_right_corner_ref, rndBox.radius };

	// Bottom
	ReferenceFrame bottom_front_left_corner_ref = rndBox.ref;
	bottom_front_left_corner_ref.Translate(Vector3RotateByQuaternion({ -rndBox.extents.x, -rndBox.extents.y, rndBox.extents.z }, bottom_front_left_corner_ref.q));
	bottom_front_left_corner_ref.q = QuaternionMultiply(bottom_front_left_corner_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, PI));
	bottom_front_left_corner_ref.q = QuaternionMultiply(bottom_front_left_corner_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, PI));
	Sphere bottom_front_left_sphere = { bottom_front_left_corner_ref, rndBox.radius };

	ReferenceFrame bottom_front_right_corner_ref = rndBox.ref;
	bottom_front_right_corner_ref.Translate(Vector3RotateByQuaternion({ rndBox.extents.x, -rndBox.extents.y, rndBox.extents.z }, bottom_front_right_corner_ref.q));
	bottom_front_right_corner_ref.q = QuaternionMultiply(bottom_front_right_corner_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, PI));
	bottom_front_right_corner_ref.q = QuaternionMultiply(bottom_front_right_corner_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, PI / 2));
	Sphere bottom_front_right_sphere = { bottom_front_right_corner_ref, rndBox.radius };

	ReferenceFrame bottom_back_left_corner_ref = rndBox.ref;
	bottom_back_left_corner_ref.Translate(Vector3RotateByQuaternion({ -rndBox.extents.x, -rndBox.extents.y, -rndBox.extents.z }, bottom_back_left_corner_ref.q));
	bottom_back_left_corner_ref.q = QuaternionMultiply(bottom_back_left_corner_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, PI));
	bottom_back_left_corner_ref.q = QuaternionMultiply(bottom_back_left_corner_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, -PI / 2));
	Sphere bottom_back_left_sphere = { bottom_back_left_corner_ref, rndBox.radius };

	ReferenceFrame bottom_back_right_corner_ref = rndBox.ref;
	bottom_back_right_corner_ref.Translate(Vector3RotateByQuaternion({ rndBox.extents.x, -rndBox.extents.y, -rndBox.extents.z }, bottom_back_right_corner_ref.q));
	bottom_back_right_corner_ref.q = QuaternionMultiply(bottom_back_right_corner_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, PI));
	bottom_back_right_corner_ref.q = QuaternionMultiply(bottom_back_right_corner_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, 0));
	Sphere bottom_back_right_sphere = { bottom_back_right_corner_ref, rndBox.radius };

	#pragma endregion

	#pragma region horizontal cylinders

	// Front
	ReferenceFrame front_top_horizontal_ref = rndBox.ref;
	front_top_horizontal_ref.Translate(Vector3RotateByQuaternion({ 0 ,rndBox.extents.y, rndBox.extents.z }, front_top_horizontal_ref.q));
	front_top_horizontal_ref.q = QuaternionMultiply(front_top_horizontal_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI / 2));
	Cylinder front_top_horizontal_cylinder = { front_top_horizontal_ref, rndBox.extents.x, rndBox.radius };

	ReferenceFrame front_bottom_horizontal_ref = rndBox.ref;
	front_bottom_horizontal_ref.Translate(Vector3RotateByQuaternion({ 0 ,-rndBox.extents.y, rndBox.extents.z }, front_bottom_horizontal_ref.q));
	front_bottom_horizontal_ref.q = QuaternionMultiply(front_bottom_horizontal_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI / 2));
	front_bottom_horizontal_ref.q = QuaternionMultiply(front_bottom_horizontal_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, -PI / 2));
	Cylinder front_bottom_horizontal_cylinder = { front_bottom_horizontal_ref, rndBox.extents.x, rndBox.radius };

	// Back
	ReferenceFrame back_top_horizontal_ref = rndBox.ref;
	back_top_horizontal_ref.Translate(Vector3RotateByQuaternion({ 0 ,rndBox.extents.y, -rndBox.extents.z }, back_top_horizontal_ref.q));
	back_top_horizontal_ref.q = QuaternionMultiply(back_top_horizontal_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI / 2));
	back_top_horizontal_ref.q = QuaternionMultiply(back_top_horizontal_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, PI / 2));
	Cylinder back_top_horizontal_cylinder = { back_top_horizontal_ref, rndBox.extents.x, rndBox.radius };

	ReferenceFrame back_bottom_horizontal_ref = rndBox.ref;
	back_bottom_horizontal_ref.Translate(Vector3RotateByQuaternion({ 0 ,-rndBox.extents.y, -rndBox.extents.z }, back_bottom_horizontal_ref.q));
	back_bottom_horizontal_ref.q = QuaternionMultiply(back_bottom_horizontal_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI / 2));
	back_bottom_horizontal_ref.q = QuaternionMultiply(back_bottom_horizontal_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, PI));
	Cylinder back_bottom_horizontal_cylinder = { back_bottom_horizontal_ref, rndBox.extents.x, rndBox.radius };

	// Left
	ReferenceFrame left_top_horizontal_ref = rndBox.ref;
	left_top_horizontal_ref.Translate(Vector3RotateByQuaternion({ -rndBox.extents.x ,rndBox.extents.y, 0 }, left_top_horizontal_ref.q));
	left_top_horizontal_ref.q = QuaternionMultiply(left_top_horizontal_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI / 2));
	left_top_horizontal_ref.q = QuaternionMultiply(left_top_horizontal_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, PI / 2));
	left_top_horizontal_ref.q = QuaternionMultiply(left_top_horizontal_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, PI / 2));
	Cylinder left_top_horizontal_cylinder = { left_top_horizontal_ref, rndBox.extents.z, rndBox.radius };

	ReferenceFrame left_bottom_horizontal_ref = rndBox.ref;
	left_bottom_horizontal_ref.Translate(Vector3RotateByQuaternion({ -rndBox.extents.x ,-rndBox.extents.y, 0 }, left_bottom_horizontal_ref.q));
	left_bottom_horizontal_ref.q = QuaternionMultiply(left_bottom_horizontal_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI / 2));
	left_bottom_horizontal_ref.q = QuaternionMultiply(left_bottom_horizontal_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, PI / 2));
	left_bottom_horizontal_ref.q = QuaternionMultiply(left_bottom_horizontal_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, PI));
	Cylinder left_bottom_horizontal_cylinder = { left_bottom_horizontal_ref, rndBox.extents.z, rndBox.radius };

	// Right
	ReferenceFrame right_top_horizontal_ref = rndBox.ref;
	right_top_horizontal_ref.Translate(Vector3RotateByQuaternion({ rndBox.extents.x ,rndBox.extents.y, 0 }, right_top_horizontal_ref.q));
	right_top_horizontal_ref.q = QuaternionMultiply(right_top_horizontal_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI / 2));
	right_top_horizontal_ref.q = QuaternionMultiply(right_top_horizontal_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, PI / 2));
	right_top_horizontal_ref.q = QuaternionMultiply(right_top_horizontal_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, 0));
	Cylinder right_top_horizontal_cylinder = { right_top_horizontal_ref, rndBox.extents.z, rndBox.radius };

	ReferenceFrame right_bottom_horizontal_ref = rndBox.ref;
	right_bottom_horizontal_ref.Translate(Vector3RotateByQuaternion({ rndBox.extents.x ,-rndBox.extents.y, 0 }, right_bottom_horizontal_ref.q));
	right_bottom_horizontal_ref.q = QuaternionMultiply(right_bottom_horizontal_ref.q, QuaternionFromAxisAngle({ 0,0,1 }, PI / 2));
	right_bottom_horizontal_ref.q = QuaternionMultiply(right_bottom_horizontal_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, PI / 2));
	right_bottom_horizontal_ref.q = QuaternionMultiply(right_bottom_horizontal_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, -PI / 2));
	Cylinder right_bottom_horizontal_cylinder = { right_bottom_horizontal_ref, rndBox.extents.z, rndBox.radius };

	#pragma endregion

	Quad faces[6] =
	{
		front_quad,
		back_quad,
		left_quad,
		right_quad,
		top_quad,
		bottom_quad
	};

	Cylinder cylinders[12] =
	{
		left_top_horizontal_cylinder,
		left_bottom_horizontal_cylinder,
		right_top_horizontal_cylinder,
		right_bottom_horizontal_cylinder,
		front_top_horizontal_cylinder,
		front_bottom_horizontal_cylinder,
		back_top_horizontal_cylinder,
		back_bottom_horizontal_cylinder,

		left_front_vertical_cylinder,
		left_back_vertical_cylinder,
		right_front_vertical_cylinder,
		right_back_vertical_cylinder
	};

	Sphere spheres[8] = {
		top_front_left_sphere,
		top_front_right_sphere,
		bottom_front_left_sphere,
		bottom_front_right_sphere,
		top_back_left_sphere,
		top_back_right_sphere,
		bottom_back_left_sphere,
		bottom_back_right_sphere,
	};

	Vector3 tempInterPt;
	Vector3 tempInterNormal;
	t = INFINITY;
	float tempT;
	bool intersect = false;

	// Vérification sur les faces
	for (int i = 0; i < 6; i++)
	{
		if (IntersectSegmentQuad(seg, faces[i], tempT, tempInterPt, tempInterNormal))
		{
			intersect = true;
			if (tempT < t)
			{
				t = tempT;
				interPt = tempInterPt;
				interNormal = tempInterNormal;
			}
		}
	}

	// Vérification sur les arrêtes
	for (int i = 0; i < 12; i++)
	{
		if (IntersectSegmentCylinder(seg, cylinders[i], tempT, tempInterPt, tempInterNormal))
		{
			intersect = true;
			if (tempT < t)
			{
				t = tempT;
				interPt = tempInterPt;
				interNormal = tempInterNormal;
			}
		}
	}

	// Vérification sur les coins
	for (int i = 0; i < 8; i++)
	{
		if (IntersectSegmentSphere(seg, spheres[i], tempT, tempInterPt, tempInterNormal))
		{
			intersect = true;
			if (tempT < t)
			{
				t = tempT;
				interPt = tempInterPt;
				interNormal = tempInterNormal;
			}
		}
	}

	return intersect;
}
#pragma endregion

#pragma region Collisions

bool IsMax(float x, float tab[], int n)
{
	for (int i = 0; i < n; i++)
	{
		if (tab[i] > x) return false;
	}
	return true;
}

bool IsMin(float x, float tab[], int n)
{
	for (int i = 0; i < n; i++)
	{
		if (tab[i] < x) return false;
	}
	return true;
}

char GetNearestFace(Vector3 pos, Box box)
{
	Vector3 localPos = GlobalToLocalPos(pos, box.ref);

	float coefs[6] =
	{
		localPos.x / box.extents.x,  // right
		localPos.x / -box.extents.x, // left
		localPos.y / box.extents.y,	 // top
		localPos.y / -box.extents.y, // bottom
		localPos.z / box.extents.z,  // front
		localPos.z / -box.extents.z  // back
	};

	if (IsMax(coefs[0], coefs, 6)) return 'R';
	else if (IsMax(coefs[1], coefs, 6)) return 'L';
	else if (IsMax(coefs[2], coefs, 6)) return 'T';
	else if (IsMax(coefs[3], coefs, 6)) return 'B';
	else if (IsMax(coefs[4], coefs, 6)) return 'f';
	else if (IsMax(coefs[5], coefs, 6)) return 'b';

	return '0';
}

bool GetSphereNewPositionAndVelocityIfCollidingWithRoundedBox(Sphere sphere, RoundedBox rndBox, Vector3 velocity, float deltaTime, float& colT, Vector3& colSpherePos, Vector3& colNormal, Vector3& newPosition, Vector3& newVelocity)
{
	RoundedBox minkowski = { rndBox.ref, rndBox.extents, rndBox.radius + sphere.radius };

	Vector3 A = sphere.ref.origin;
	Vector3 B = Vector3Add(sphere.ref.origin, Vector3Scale(velocity, deltaTime));
	Segment AB = { A, B };

	Box OBB = { minkowski.ref, {minkowski.extents.x + minkowski.radius, minkowski.extents.y + minkowski.radius, minkowski.extents.z + minkowski.radius} };

	// Check faces
	if (IsPointInsideBox(OBB, B) && IntersectSegmentBox(AB, OBB, colT, colSpherePos, colNormal))
	{
		Vector3 localColPt = GlobalToLocalPos(colSpherePos, OBB.ref);
		char nearestFace = GetNearestFace(colSpherePos, OBB);

		bool inX = false;
		bool inY = false;
		bool inZ = false;

		if (nearestFace == 'L' || nearestFace == 'R')
		{
			inX = fabs(localColPt.x) <= rndBox.extents.x + rndBox.radius + sphere.radius + EPSILON;
			inY = fabs(localColPt.y) <= rndBox.extents.y + rndBox.radius + EPSILON;
			inZ = fabs(localColPt.z) <= rndBox.extents.z + rndBox.radius + EPSILON;
		}
		else if (nearestFace == 'T' || nearestFace == 'B')
		{
			inX = fabs(localColPt.x) <= rndBox.extents.x + rndBox.radius + EPSILON;
			inY = fabs(localColPt.y) <= rndBox.extents.y + rndBox.radius + sphere.radius + EPSILON;
			inZ = fabs(localColPt.z) <= rndBox.extents.z + rndBox.radius + EPSILON;
		}
		else if (nearestFace == 'f' || nearestFace == 'b')
		{
			inX = fabs(localColPt.x) <= rndBox.extents.x + rndBox.radius + EPSILON;
			inY = fabs(localColPt.y) <= rndBox.extents.y + rndBox.radius + EPSILON;
			inZ = fabs(localColPt.z) <= rndBox.extents.z + rndBox.radius + sphere.radius + EPSILON;
		}

		if (inX && inY && inZ) {
			colSpherePos = Vector3Add(colSpherePos, Vector3Scale(Vector3Normalize(colNormal), -sphere.radius));
			newVelocity = Vector3Reflect(velocity, colNormal); // Changement de direction
			newPosition = Vector3Add(sphere.ref.origin, Vector3Scale(newVelocity, deltaTime * (1 - colT)));  // Ajout de la distance restante dans la bonne direction
			return true;
		}
	}
	// Vérification du reste seulement si nécéssaire
	if (IsPointInsideBox(OBB, B) && IntersectSegmentRoundedBox(AB, minkowski, colT, colSpherePos, colNormal))
	{
		colSpherePos = Vector3Add(colSpherePos, Vector3Scale(Vector3Normalize(colNormal), -sphere.radius));
		newVelocity = Vector3Reflect(velocity, colNormal); // Changement de direction
		newPosition = Vector3Add(sphere.ref.origin, Vector3Scale(newVelocity, deltaTime * (1 - colT)));  // Ajout de la distance restante dans la bonne direction
		return true;
	}
	return false;
}

bool GetSphereNewPositionAndVelocityIfCollidingWithRoundedBoxes(Sphere sphere, const std::vector<RoundedBox>& rndBoxes, Vector3 velocity, float deltaTime, float& colT, Vector3& colSpherePos, Vector3& colNormal, Vector3& newPosition, Vector3& newVelocity, Color& newColor)
{
	colT = FLT_MAX;
	bool collided = false;
	for (const RoundedBox& rndBox : rndBoxes)
	{
		float t;
		Vector3 colPt, normal, newPos, newVel;
		// Collisionne avec toutes les boxs et retient celle la plus proche de la sphère
		if (GetSphereNewPositionAndVelocityIfCollidingWithRoundedBox(sphere, rndBox, velocity, deltaTime, t, colPt, normal, newPos, newVel))
		{
			if (t < colT)
			{
				colT = t;
				colSpherePos = colPt;
				colNormal = normal;
				newPosition = newPos;
				newVelocity = newVel;
				newColor = rndBox.color;
				collided = true;
			}
		}
	}
	return collided;
}

bool GetSphereNewPositionAndVelocityIfMultiCollidingWithRoundedBoxes(Sphere sphere, const std::vector<RoundedBox>& rndBoxes, Vector3 velocity, float rotInertia, Vector3 angularMomentum, float deltaTime, int nMaxSuccessiveCollisions, Vector3& newPosition, Vector3& newVelocity, Vector3 newAngularMomentum, Vector3& colNormal, Vector3& colSpherePos, Color& newColor)
{
	float colT;
	bool collide = GetSphereNewPositionAndVelocityIfCollidingWithRoundedBoxes(sphere, rndBoxes, velocity, deltaTime, colT, colSpherePos, colNormal, newPosition, newVelocity, newColor);
	int collideCount = 0;
	Sphere tempSphere = sphere;
	// Appel successif de la fonction avec la nouvelle position / velocité tant qu'il y a des intersections à gérer
	while (collideCount < nMaxSuccessiveCollisions&& collide)
	{
		collideCount++;
		tempSphere.ref.origin = newPosition;
		collide = GetSphereNewPositionAndVelocityIfCollidingWithRoundedBoxes(tempSphere, rndBoxes, newVelocity, deltaTime, colT, colSpherePos, colNormal, newPosition, newVelocity, newColor);
	}
	return collideCount > 0;
}

#pragma endregion

/// <summary>
/// Génère une grille aléatoire configurable d'obstacles
/// </summary>
/// <param name="w">Largeur</param>
/// <param name="h">Hauteur</param>
/// <param name="d">Densité d'obstacles par ligne ou colonne</param>
/// <param name="n">Nombre d'obstacles</param>
std::vector<RoundedBox> GenerateRandomObstacleGrid(float width, float depth, int d, int n = 10000, float y = 3, bool random = true)
{
	// Ajout des murs de base
	std::vector<RoundedBox> obstacles =
	{
		{ { { 0,-1,0 }, QuaternionIdentity()}, { 20, 0.2, 20 }, 0, LIGHTGRAY },
		{ { { 0,20,0 }, QuaternionIdentity()}, { 20, 0.2, 20 }, 0, BLANK},

		{ { { -20,8,0 }, QuaternionFromAxisAngle({0,1,0}, PI / 2)}, {20, 10, 0.2}, 0, GRAY},
		{ { { 20,8,0 }, QuaternionFromAxisAngle({0,1,0}, -PI / 2)}, {20, 10, 0.2}, 0, GRAY},
		{ { { 0,8,20 }, QuaternionFromAxisAngle({0,0,0}, 0)}, {20, 10, 0.2}, 0, GRAY},
		{ { { 0,8,-20 }, QuaternionFromAxisAngle({0,0,0}, 0)}, {20, 10, 0.2}, 0, GRAY},
	};

	std::vector<Vector3> coords = {};
	ReferenceFrame ref;
	Vector3 rotAxis, coord;
	Color color;


	float val, rot, radius, x, z;
	int index;
	float boxSize = (min(width, depth) / d - 1) / 2;
	if (boxSize <= 0) return obstacles;
	if (n > d * d) n = d * d;

	// Génère les coordonnées possibles
	for (int i = 0; i < d; i++)
	{
		x = (-width / 2) + (width / d) / 2 + (width / d) * i;
		for (int j = 0; j < d; j++)
		{
			z = (-depth / 2) + (depth / d) / 2 + (depth / d) * j;
			coords.push_back({ x, y, z });
		}
	}

	// Choix aléatoire (ou non) de coordonnées et création des obstacles
	// Une coordonnée n'est choisie qu'une seule fois
	for (int i = 0; i < n; i++)
	{
		if (random)
		{
			index = rand() % coords.size();
			coord = coords[index];
			coords.erase(coords.begin() + index);

			rot = rand();
			rotAxis = { (float)(rand() % 2), (float)(rand() % 2), (float)(rand() % 2) };
			ref = { coord, QuaternionFromAxisAngle(rotAxis, rot) };

			radius = boxSize * (rand() % 100) / 100.f;
			color = { (unsigned char)(rand() % 255),(unsigned char)(rand() % 255),(unsigned char)(rand() % 255), 255 };
			obstacles.push_back({ ref, {boxSize - radius, boxSize - radius, boxSize - radius}, radius, color });
		}
		// Possibilité d'avoir une grille prédéfinie pour une meilleure fiabilité des benchmarks
		else
		{
			coord = coords[i];

			rot = (i % 4) * PI / 8;
			rotAxis = { (float)(i % 3),(float)(i % 6),(float)(i % 9) };
			ref = { coord, QuaternionFromAxisAngle(rotAxis, rot) };

			radius = (float)(i % 9) / 9.f;
			color = SKYBLUE;
			obstacles.push_back({ ref, {boxSize - radius, boxSize - radius, boxSize - radius}, radius, color });
		}
	}

	return obstacles;
}

// Utilisé pour les contrôles du joueur
void Game(RoundedBox& player, Sphere ball, float deltaTime)
{
	static float moveSpeed = 10;

	// Contrôles (clavier)
	Vector3 vel = { 0, 0, 0 };
	if (IsKeyDown(KEY_LEFT)) vel.x = -1;
	if (IsKeyDown(KEY_RIGHT)) vel.x = 1;
	if (IsKeyDown(KEY_UP)) vel.z = -1;
	if (IsKeyDown(KEY_DOWN)) vel.z = 1;
	vel = Vector3Scale(Vector3Normalize(vel), deltaTime * moveSpeed);
	player.ref.origin = Vector3Add(player.ref.origin, vel);

	// Contrôles (souris)
	static Vector2 prevMousePos = { 0, 0 };
	Vector2 mousePos = GetMousePosition();
	Vector2 mouseVect = Vector2Subtract(mousePos, prevMousePos);
	prevMousePos = GetMousePosition();
	player.ref.RotateByQuaternion(QuaternionFromAxisAngle({ mouseVect.y, 0, -mouseVect.x }, PI * deltaTime));
}

// Génère une couleur sinusoidale avec le temps (va de noir à la couleur choisie)
Color GetSinusoidalColor(Color color, float time)
{
	static int colorSpeed = 5;
	float colorPercentage = 0.5 + (sin(time * colorSpeed) + 1) / 4;
	unsigned char r = colorPercentage * color.r;
	unsigned char g = colorPercentage * color.g;
	unsigned char b = colorPercentage * color.b;
	return { r, g, b, 255 };
}

int main(int argc, char* argv[])
{
	#pragma region Init

	// Initialization
	//--------------------------------------------------------------------------------------
	float screenSizeCoef = .9f;
	const int screenWidth = 1920 * screenSizeCoef;
	const int screenHeight = 1080 * screenSizeCoef;

	InitWindow(screenWidth, screenHeight, "ESIEE - E3FI - 2022 - 2023 - Maths 3D");

	SetTargetFPS(144);

	// CAMERA 
	Vector3 cameraPos = SphericalToCartesian({ 15, 0, PI / 4.f });
	Camera camera = { 0 };
	camera.position = cameraPos;
	camera.target = { 0, 0, 0 };
	camera.up = { 0, 1, 0 };
	camera.fovy = 90;
	camera.type = CAMERA_PERSPECTIVE;
	SetCameraMode(camera, CAMERA_CUSTOM);  // Set an orbital camera mode

	//--------------------------------------------------------------------------------------

	#pragma endregion

	#pragma region Setup

	Vector3 initialPos;
	Vector3 initialVel;
	Vector3 pos;
	Vector3 vel;
	Vector3 rot = { 0, 0, 0 };
	Color sphereColor = SKYBLUE;
	Color highlightColor = BLACK;

	std::vector<RoundedBox> objects;
	std::vector<RoundedBox> debugObjects;
	RoundedBox basePlayer = { { { 0,0,0 }, QuaternionIdentity() }, {3.5, 0.2, 3.5}, 1, SKYBLUE };
	RoundedBox player = basePlayer;

	int score = 0;
	bool game = true;
	if (game)
	{
		initialPos = { 0, 10, 0 };
		initialVel = { 0, 10, 0 };
		camera.position = SphericalToCartesian({ 15, 0, PI / 4.f });
		objects = { player };
	}
	else
	{
		initialPos = { 0, 15, 0 };
		initialVel = { 15, 4, 10 };
		debugObjects = GenerateRandomObstacleGrid(38, 38, 5, 10000, 3, true);
		objects = debugObjects;
	}

	// Physics
	pos = initialPos;
	vel = initialVel;
	float mass = 10;
	float radius = 1;
	float E = 0.5f * (mass * Vector3Length(initialVel) + mass * Vector3Length(GRAVITY) * initialPos.y);
	float initialE = E;
	float rotInertia = (2.f / 5.f) * (mass * powf(2, radius));
	float rotSpeed = 0.25;
	float friction = 0.95;
	Vector3 omegaC = { 0,0,0 };
	Vector3 angularMomentum = { 0,0,0 };
	int nMaxSuccessiveCollisions = 5;

	float averageFPS = 0;

	#pragma endregion

	// Main game loop
	while (!WindowShouldClose()) // Detect window close button or ESC key
	{
		// Update
		//----------------------------------------------------------------------------------
		// TODO: Update your variables here
		//----------------------------------------------------------------------------------

		float deltaTime = GetFrameTime();
		float time = (float)GetTime();
		if (!game) MyUpdateOrbitalCamera(&camera, deltaTime);

		// Draw
		//----------------------------------------------------------------------------------
		BeginDrawing();

		ClearBackground(RAYWHITE);

		BeginMode3D(camera);
		{
			#pragma region display tests

			//ReferenceFrame ref;
			//static float angle = 0;
			//angle += 0.05;
			//Vector3 axes = { 1, 1, 1 };

			////PLANE DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ camera.position.x, -10, camera.position.z },
			//	QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
			//Quad plane = { ref,{10000, 0, 10000} };
			//MyDrawQuad(plane, true, false, LIGHTGRAY);

			//// DISK DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ -20,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize(axes), angle));
			//Disk disk = { ref, 2 };
			//MyDrawDisk(disk, 15, true, true);

			////BOX DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ 0,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize(axes), angle));
			//Box box = { ref, {1, 1, 2} };
			//MyDrawBox(box, true, true, BLUE, BLACK);

			///*Vector3 sphPos = { 1.8,0,0 };
			//ref = ReferenceFrame(
			//	sphPos,
			//	QuaternionFromAxisAngle(Vector3Normalize(axes), angle));
			//Sphere sphere = { ref, .15f };
			//MyDrawSphere(sphere, 10, 10, true, true, GREEN, BLACK);
			//cout << IsPointInsideBox(box, sphPos) << "\n";*/

			////QUAD DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ -10,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize(axes), angle));
			//Quad quad = { ref,{1, 0, 2} };
			//MyDrawQuad(quad, true, true);

			//// SPHERE DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ -5,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize(axes), angle));
			//Sphere sphere = { ref, 2 };
			//MyDrawSphere(sphere, 15, 15, true, true);

			//// CYLINDER DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ 0,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize(axes), angle));
			//Cylinder cylinder = { ref, 1, 2 };
			//MyDrawCylinder(cylinder, 15, false, true, true);

			//// CYLINDER QUARTER DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ 5,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize(axes), angle));
			//CylinderQuarter cylinderQuarter = { ref, 1, 2 };
			//MyDrawCylinderQuarter(cylinderQuarter, 30, true, true, true);

			//// SPHERE CORNER DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ 10,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize(axes), angle));
			//SphereCorner sphCorner = { ref, 2 };
			//MyDrawSphereCorner(sphCorner, 10, 10, true, true);

			//// CAPSULE DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ 15,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize(axes), angle));
			//Capsule capsule = { ref, 2, 2 };
			//MyDrawCapsule(capsule, 25, 25, true, true, BLUE, WHITE);

			//// HEMISPHERE DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ 20,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize(axes), angle));
			//Hemisphere hemisphere = { ref, 2 };
			//MyDrawHemisphere(hemisphere, 15, 15, true, true);

			//// ROUNDEDBOX DISPLAY TEST
			//ref = ReferenceFrame(
			//	{ 5,8,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize(axes), angle));
			///* {ref, {Longeur Largeur Hauteur}, Radius} */
			//RoundedBox roundedBox = { ref, {5,2,3}, 2 };
			//MyDrawRoundedBox(roundedBox, 10, true, true, RED, WHITE);

#pragma endregion

			#pragma region game management

			Sphere sphere = { {pos, QuaternionFromAxisAngle(Vector3Normalize(angularMomentum), rotInertia * time * rotSpeed)}, radius, GetSinusoidalColor(sphereColor, time) };

			if (IsKeyPressed(KEY_G))
			{
				game = true;
				player = basePlayer;
				initialPos = { 0, 5, 0 };
				initialVel = { 0, 10, 0 };
				pos = initialPos;
				vel = initialVel;
				score = 0;
				camera.position = SphericalToCartesian({ 15, 0, PI / 4.f });
			}
			else if (IsKeyPressed(KEY_D))
			{
				game = false;
				debugObjects = GenerateRandomObstacleGrid(38, 38, 5, 10000, 3, true);
				objects = debugObjects;				
				initialPos = { 0, 15, 0 };
				initialVel = { 15, 4, 10 };
				pos = initialPos;
				vel = { (float)(rand() % 40 - 20), (float)(rand() % 20), (float)(rand() % 40 - 20) };
				camera.position = SphericalToCartesian({ 30,  PI / 4.f, PI / 4.f });
			}

			if (game)
			{
				objects = { player };
				Game(player, sphere, deltaTime);
				if (sphere.ref.origin.y < -20)
				{
					pos = initialPos;
					vel = { 0,10,0 };
					score = 0;
					player = basePlayer;
				}
			}
			else objects = debugObjects;

			#pragma endregion

			#pragma region Physics and drawing

			float colT;
			Vector3 colSpherePos;
			Vector3 colNormal;
			Vector3 newPos;
			Vector3 newVel;
			Vector3 colPt;
			Color newColor;

			vel = Vector3Add(vel, Vector3Scale(GRAVITY, deltaTime));
			float k = mass * Vector3Length(GRAVITY) * (sphere.ref.origin.y - initialPos.y);
			float vitesse = sqrtf(2 * (E - k) / mass);
			if (E > k) vel = Vector3Scale(Vector3Normalize(vel), vitesse);

			bool collide = GetSphereNewPositionAndVelocityIfMultiCollidingWithRoundedBoxes(sphere, objects, vel, rotInertia, angularMomentum, deltaTime, 100, newPos, newVel, angularMomentum, colNormal, colSpherePos, newColor);
			if (collide)
			{
				//E *= 0.8;
				omegaC = GlobalToLocalVect(colSpherePos, sphere.ref);
				angularMomentum = Vector3Subtract(angularMomentum, Vector3Scale(Vector3CrossProduct(omegaC, Vector3Subtract(vel, Vector3Add(Vector3Scale(colNormal, Vector3DotProduct(vel, colNormal)), Vector3CrossProduct(rot, omegaC)))), 1 * deltaTime));
				vel = newVel;
				pos = newPos;
				sphereColor = newColor;
				score += 1;
			}
			else pos = Vector3Add(pos, Vector3Scale(vel, deltaTime));

			MyDrawSphere(sphere, 15, 15, true, true, sphere.color);
			for (const RoundedBox& rndBox : objects)
			{
				if (rndBox.color.a != 0)
				{
					MyDrawRoundedBox(rndBox, 3, true, true, rndBox.color);
				}
			}

			#pragma endregion

			#pragma region Debug controls

			if (!game && IsKeyDown(KEY_SPACE))
			{
				pos = initialPos;
				vel = { (float)(rand() % 40 - 20), (float)(rand() % 20), (float)(rand() % 40 - 20) };
			}

			if (!game && IsMouseButtonDown(MOUSE_LEFT_BUTTON))
			{
				Vector3 ray = Vector3Scale(Vector3Subtract(camera.target, camera.position), 100);
				Segment shoot = { camera.position, Vector3Add(camera.position, ray) };
				for (const RoundedBox& rndBox : objects)
				{
					if (rndBox.color.a == 0) continue;
					bool cameraTest = IntersectSegmentRoundedBox(shoot, rndBox, colT, colPt, colNormal);
					if (cameraTest)
					{
						MyDrawSphere({ {colPt,QuaternionIdentity()},.1f }, 16, 8, true, false, GetSinusoidalColor(SKYBLUE, time));
						DrawLine3D(colPt, Vector3Add(Vector3Scale(colNormal, 1), colPt), GetSinusoidalColor(SKYBLUE, time));
					}
				}
			}

			#pragma endregion

			#pragma region FPS counter

			static float lastTime = time;
			static unsigned long frames = 0;
			int interval = 1;
			if (time - lastTime >= interval)
			{
				averageFPS = (float)frames / (float)interval;
				frames = 0;
				lastTime = time;
			}
			frames += 1;

			#pragma endregion
		}
		EndMode3D();

		#pragma region Text

		static float shadowOffset = 1;
		static int fontSize = 30;
		static float spacing = 3;
		static Color fontColor = GREEN;
		static Color shadowColor = BLACK;

		DrawTextEx(GetFontDefault(), FormatText("Avg FPS: %.2f", averageFPS), { 10 + shadowOffset, 10 + shadowOffset }, fontSize, spacing, shadowColor);
		DrawTextEx(GetFontDefault(), FormatText("Avg FPS: %.2f", averageFPS), { 10, 10 }, fontSize, spacing, fontColor);

		DrawTextEx(GetFontDefault(), FormatText("Velocity: {%0.2f, %0.2f, %0.2f}", vel.x, vel.y, vel.z), { 10 + shadowOffset, 10 + 50 + shadowOffset }, fontSize, spacing, shadowColor);
		DrawTextEx(GetFontDefault(), FormatText("Velocity: {%0.2f, %0.2f, %0.2f}", vel.x, vel.y, vel.z), { 10, 10 + 50 }, fontSize, spacing, fontColor);

		DrawTextEx(GetFontDefault(), FormatText("Angular momentum: {%0.2f, %0.2f, %0.2f}", angularMomentum.x, angularMomentum.y, angularMomentum.z), { 10 + shadowOffset, 10 + 100 + shadowOffset }, fontSize, spacing, shadowColor);
		DrawTextEx(GetFontDefault(), FormatText("Angular momentum: {%0.2f, %0.2f, %0.2f}", angularMomentum.x, angularMomentum.y, angularMomentum.z), { 10, 10 + 100 }, fontSize, spacing, fontColor);

		DrawTextEx(GetFontDefault(), FormatText("D: debug mode / regenerate grid"), { 10 + shadowOffset, 10 + 150 + shadowOffset }, fontSize, spacing, shadowColor);
		DrawTextEx(GetFontDefault(), FormatText("D: debug mode / regenerate grid"), { 10, 10 + 150 }, fontSize, spacing, !game ? DARKGRAY : SKYBLUE);

		DrawTextEx(GetFontDefault(), FormatText("G: game mode"), { 10 + shadowOffset, 10 + 200 + shadowOffset }, fontSize, spacing, shadowColor);
		DrawTextEx(GetFontDefault(), FormatText("G: game mode"), { 10, 10 + 200 }, fontSize, spacing, game ? DARKGRAY : SKYBLUE);

		if (game == true)
		{
			bool arrowsPressed = IsKeyDown(KEY_LEFT) || IsKeyDown(KEY_RIGHT) || IsKeyDown(KEY_UP) || IsKeyDown(KEY_DOWN);

			static Vector2 prevMousePos = { 0, 0 };
			Vector2 mousePos = GetMousePosition();
			Vector2 mouseVect = Vector2Subtract(mousePos, prevMousePos);
			bool mouseMoving = Vector2Length(mouseVect) != 0;
			prevMousePos = GetMousePosition();
			DrawTextEx(GetFontDefault(), FormatText("Arrows keys: move"), { 10 + shadowOffset, 10 + 250 + shadowOffset }, fontSize, spacing, shadowColor);
			DrawTextEx(GetFontDefault(), FormatText("Arrows keys: move"), { 10, 10 + 250 }, fontSize, spacing, arrowsPressed ? DARKGRAY : SKYBLUE);

			DrawTextEx(GetFontDefault(), FormatText("Mouse: rotate"), { 10 + shadowOffset, 10 + 300 + shadowOffset }, fontSize, spacing, shadowColor);
			DrawTextEx(GetFontDefault(), FormatText("Mouse: rotate"), { 10, 10 + 300 }, fontSize, spacing, mouseMoving ? DARKGRAY : SKYBLUE);

			DrawTextEx(GetFontDefault(), FormatText("Score: %d", score), { 750 + shadowOffset, 50 + shadowOffset }, 40, spacing, shadowColor);
			DrawTextEx(GetFontDefault(), FormatText("Score: %d", score), { 750, 50 }, 40, spacing, SKYBLUE);

		}
		else
		{
			DrawTextEx(GetFontDefault(), FormatText("Left click: shoot ray"), { 10 + shadowOffset, 10 + 250 + shadowOffset }, fontSize, spacing, shadowColor);
			DrawTextEx(GetFontDefault(), FormatText("Left click: shoot ray"), { 10, 10 + 250 }, fontSize, spacing, IsMouseButtonDown(MOUSE_LEFT_BUTTON) ? DARKGRAY : SKYBLUE);

			DrawTextEx(GetFontDefault(), FormatText("Right click: move camera"), { 10 + shadowOffset, 10 + 300 + shadowOffset }, fontSize, spacing, shadowColor);
			DrawTextEx(GetFontDefault(), FormatText("Right click: move camera"), { 10, 10 + 300 }, fontSize, spacing, IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? DARKGRAY : SKYBLUE);

			DrawTextEx(GetFontDefault(), FormatText("Space: reset ball"), { 10 + shadowOffset, 10 + 350 + shadowOffset }, fontSize, spacing, shadowColor);
			DrawTextEx(GetFontDefault(), FormatText("Space: reset ball"), { 10, 10 + 350 }, fontSize, spacing, IsKeyDown(KEY_SPACE) ? DARKGRAY : SKYBLUE);
		}

		#pragma endregion

		EndDrawing();
		//----------------------------------------------------------------------------------
	}

	// De-Initialization
	//--------------------------------------------------------------------------------------   
	CloseWindow();        // Close window and OpenGL context
	//--------------------------------------------------------------------------------------

	return 0;
}