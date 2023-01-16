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
	ReferenceFrame ref;
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
	Vector3 localVect = Vector3Subtract(globalVect, localRef.origin);

	Vector3 rotationAxis;
	float angle;
	QuaternionToAxisAngle(localRef.q, &rotationAxis, &angle);

	Quaternion invertedQ = QuaternionFromAxisAngle(Vector3Normalize(rotationAxis), -angle);
	return Vector3RotateByQuaternion(localVect, invertedQ);
}

Vector3 LocalToGlobalPos(Vector3 localPos, ReferenceFrame localRef)
{
	localPos = Vector3RotateByQuaternion(localPos, localRef.q);
	return Vector3Add(localPos, localRef.origin);
}

Vector3 GlobalToLocalPos(Vector3 globalPos, ReferenceFrame localRef)
{
	Vector3 localPos = Vector3Subtract(globalPos, localRef.origin);

	Vector3 rotationAxis;
	float angle;
	QuaternionToAxisAngle(localRef.q, &rotationAxis, &angle);

	Quaternion invertedQ = QuaternionFromAxisAngle(Vector3Normalize(rotationAxis), -angle);
	return Vector3RotateByQuaternion(localPos, invertedQ);
}

Vector3 ProjectedPointOnLine(Vector3 linePt, Vector3 lineUnitDir, Vector3 pt) {
	Vector3 v = Vector3Subtract(pt, linePt);  // Vecteur allant de linePt à pt
	float d = Vector3DotProduct(v, lineUnitDir);  // Distance de pt au point projeté
	return Vector3Add(linePt, Vector3Scale(lineUnitDir, d));
}


bool IsPointInsideBox(Box box, Vector3 globalPt)
{
	Vector3 localPt = GlobalToLocalPos(globalPt, box.ref);
	bool w = localPt.x > -box.extents.x && localPt.x < box.extents.x;
	bool h = localPt.y > -box.extents.y && localPt.y < box.extents.y;
	bool d = localPt.z > -box.extents.z && localPt.z < box.extents.z;
	return (w && h && d);
}

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

#pragma region sphereCorner

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

#pragma region hemisphere

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

void MyDrawHemisphere(Hemisphere hemisphere, int nMeridians, int nParallels, bool drawPolygon = true, bool drawWireframe = true, Color polygonColor = LIGHTGRAY, Color wireframeColor = DARKGRAY)
{
	if (drawPolygon) MyDrawPolygonHemisphere(hemisphere, nMeridians, nParallels, polygonColor);
	if (drawWireframe) MyDrawWireframeHemisphere(hemisphere, nMeridians, nParallels, wireframeColor);
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

#pragma region capsule

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

	roundedBox.extents.x /= 2;
	roundedBox.extents.y /= 2;
	roundedBox.extents.z /= 2;
	roundedBox.radius /= 2;

#pragma region Pavé

	ReferenceFrame top_ref = ReferenceFrame(
		{ 0,roundedBox.extents.z + roundedBox.radius,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	Quad top_quad = { top_ref,{roundedBox.extents.x,0,roundedBox.extents.y} };
	MyDrawPolygonQuad(top_quad, color);

	ReferenceFrame bottom_ref = ReferenceFrame(
		{ 0,-roundedBox.extents.z - roundedBox.radius,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI));
	Quad bottom_quad = { bottom_ref,{roundedBox.extents.x,0,roundedBox.extents.y} };
	MyDrawPolygonQuad(bottom_quad, color);

	ReferenceFrame front_ref = ReferenceFrame(
		{ 0,0,roundedBox.extents.y + roundedBox.radius },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	Quad front_quad = { front_ref,{roundedBox.extents.x,0,roundedBox.extents.z} };
	MyDrawPolygonQuad(front_quad, color);

	ReferenceFrame back_ref = ReferenceFrame(
		{ 0,0,-roundedBox.extents.y - roundedBox.radius },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2));
	Quad back_quad = { back_ref,{roundedBox.extents.x,0,roundedBox.extents.z} };
	MyDrawPolygonQuad(back_quad, color);

	ReferenceFrame left_ref = ReferenceFrame(
		{ -roundedBox.extents.x - roundedBox.radius,0,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2)));
	Quad left_quad = { left_ref,{roundedBox.extents.y,0,roundedBox.extents.z} };
	MyDrawPolygonQuad(left_quad, color);

	ReferenceFrame right_ref = ReferenceFrame(
		{ roundedBox.extents.x + roundedBox.radius,0,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2)));
	Quad right_quad = { right_ref,{roundedBox.extents.y,0,roundedBox.extents.z} };
	MyDrawPolygonQuad(right_quad, color);

#pragma endregion

#pragma region Quart de cylindre vertical

	ReferenceFrame first_vertical_ref(
		{ roundedBox.extents.x,0,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	CylinderQuarter firstVerticalCylinderQuarter = { first_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(firstVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_vertical_ref(
		{ -roundedBox.extents.x,0,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), 3 * PI / 2));
	CylinderQuarter secondVerticalCylinderQuarter = { second_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(secondVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_vertical_ref(
		{ -roundedBox.extents.x,0,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	CylinderQuarter thirdVerticalCylinderQuarter = { third_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(thirdVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_vertical_ref(
		{ roundedBox.extents.x,0,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2));
	CylinderQuarter fourthVerticalCylinderQuarter = { fourth_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(fourthVerticalCylinderQuarter, nSectors, false, color);

#pragma endregion

#pragma region Quart de cylindre horizontal face supérieure

	ReferenceFrame first_superior_horizontal_ref(
		{ 0,roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2));
	CylinderQuarter firstSuperiorHorizontalCylinderQuarter = { first_superior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(firstSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_superior_horizontal_ref(
		{ 0,roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2)));
	CylinderQuarter secondSuperiorHorizontalCylinderQuarter = { second_superior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(secondSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_superior_horizontal_ref(
		{ roundedBox.extents.x,roundedBox.extents.z,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2));
	CylinderQuarter thirdSuperiorHorizontalCylinderQuarter = { third_superior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(thirdSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_superior_horizontal_ref(
		{ -roundedBox.extents.x,roundedBox.extents.z,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2)));
	CylinderQuarter fourthSuperiorHorizontalCylinderQuarter = { fourth_superior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(fourthSuperiorHorizontalCylinderQuarter, nSectors, false, color);

#pragma endregion

#pragma region Quart de cylindre horizontal face inférieur

	ReferenceFrame first_inferior_horizontal_ref(
		{ 0,-roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2));
	CylinderQuarter firstInferiorHorizontalCylinderQuarter = { first_inferior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(firstInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_inferior_horizontal_ref(
		{ 0,-roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2)));
	CylinderQuarter secondInferiorHorizontalCylinderQuarter = { second_inferior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(secondInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_inferior_horizontal_ref(
		{ roundedBox.extents.x,-roundedBox.extents.z,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	CylinderQuarter thirdInferiorHorizontalCylinderQuarter = { third_inferior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(thirdInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_inferior_horizontal_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.z,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2)));
	CylinderQuarter fourthInferiorHorizontalCylinderQuarter = { fourth_inferior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawPolygonCylinderQuarter(fourthInferiorHorizontalCylinderQuarter, nSectors, false, color);

#pragma endregion

#pragma region Coin face supérieur

	ReferenceFrame first_superior_corner_ref(
		{ roundedBox.extents.x,roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	SphereCorner firstSuperiorCorner = { first_superior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(firstSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame second_superior_corner_ref(
		{ roundedBox.extents.x,roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2));
	SphereCorner secondSuperiorCorner = { second_superior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(secondSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame third_superior_corner_ref(
		{ -roundedBox.extents.x,roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2));
	SphereCorner thirdSuperiorCorner = { third_superior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(thirdSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame fourth_superior_corner_ref(
		{ -roundedBox.extents.x,roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	SphereCorner fourthSuperiorCorner = { fourth_superior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(fourthSuperiorCorner, nSectors, nSectors, color);

#pragma endregion

#pragma region Coin face inférieur

	ReferenceFrame first_inferior_corner_ref(
		{ roundedBox.extents.x,-roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	SphereCorner firstInferiorCorner = { first_inferior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(firstInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame second_inferior_corner_ref(
		{ roundedBox.extents.x,-roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2)));
	SphereCorner secondInferiorCorner = { second_inferior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(secondInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame third_inferior_corner_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.z,roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI)));
	SphereCorner thirdInferiorCorner = { third_inferior_corner_ref, roundedBox.radius };
	MyDrawPolygonSphereCorner(thirdInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame fourth_inferior_corner_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.z,-roundedBox.extents.y },
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

	roundedBox.extents.x /= 2;
	roundedBox.extents.y /= 2;
	roundedBox.extents.z /= 2;
	roundedBox.radius /= 2;

#pragma region Pavé

	ReferenceFrame top_ref = ReferenceFrame(
		{ 0,roundedBox.extents.z + roundedBox.radius,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	Quad top_quad = { top_ref,{roundedBox.extents.x,0,roundedBox.extents.y} };
	MyDrawWireframeQuad(top_quad, color);

	ReferenceFrame bottom_ref = ReferenceFrame(
		{ 0,-roundedBox.extents.z - roundedBox.radius,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI));
	Quad bottom_quad = { bottom_ref,{roundedBox.extents.x,0,roundedBox.extents.y} };
	MyDrawWireframeQuad(bottom_quad, color);

	ReferenceFrame front_ref = ReferenceFrame(
		{ 0,0,roundedBox.extents.y + roundedBox.radius },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	Quad front_quad = { front_ref,{roundedBox.extents.x,0,roundedBox.extents.z} };
	MyDrawWireframeQuad(front_quad, color);

	ReferenceFrame back_ref = ReferenceFrame(
		{ 0,0,-roundedBox.extents.y - roundedBox.radius },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2));
	Quad back_quad = { back_ref,{roundedBox.extents.x,0,roundedBox.extents.z} };
	MyDrawWireframeQuad(back_quad, color);

	ReferenceFrame left_ref = ReferenceFrame(
		{ -roundedBox.extents.x - roundedBox.radius,0,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2)));
	Quad left_quad = { left_ref,{roundedBox.extents.y,0,roundedBox.extents.z} };
	MyDrawWireframeQuad(left_quad, color);

	ReferenceFrame right_ref = ReferenceFrame(
		{ roundedBox.extents.x + roundedBox.radius,0,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2)));
	Quad right_quad = { right_ref,{roundedBox.extents.y,0,roundedBox.extents.z} };
	MyDrawWireframeQuad(right_quad, color);

#pragma endregion

#pragma region Quart de cylindre vertical

	ReferenceFrame first_vertical_ref(
		{ roundedBox.extents.x,0,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	CylinderQuarter firstVerticalCylinderQuarter = { first_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(firstVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_vertical_ref(
		{ -roundedBox.extents.x,0,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), 3 * PI / 2));
	CylinderQuarter secondVerticalCylinderQuarter = { second_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(secondVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_vertical_ref(
		{ -roundedBox.extents.x,0,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	CylinderQuarter thirdVerticalCylinderQuarter = { third_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(thirdVerticalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_vertical_ref(
		{ roundedBox.extents.x,0,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2));
	CylinderQuarter fourthVerticalCylinderQuarter = { fourth_vertical_ref, roundedBox.extents.z, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(fourthVerticalCylinderQuarter, nSectors, false, color);

#pragma endregion

#pragma region Quart de cylindre horizontal face supérieure

	ReferenceFrame first_superior_horizontal_ref(
		{ 0,roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2));
	CylinderQuarter firstSuperiorHorizontalCylinderQuarter = { first_superior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(firstSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_superior_horizontal_ref(
		{ 0,roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2)));
	CylinderQuarter secondSuperiorHorizontalCylinderQuarter = { second_superior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(secondSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_superior_horizontal_ref(
		{ roundedBox.extents.x,roundedBox.extents.z,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2));
	CylinderQuarter thirdSuperiorHorizontalCylinderQuarter = { third_superior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(thirdSuperiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_superior_horizontal_ref(
		{ -roundedBox.extents.x,roundedBox.extents.z,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2)));
	CylinderQuarter fourthSuperiorHorizontalCylinderQuarter = { fourth_superior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(fourthSuperiorHorizontalCylinderQuarter, nSectors, false, color);

#pragma endregion

#pragma region Quart de cylindre horizontal face inférieur

	ReferenceFrame first_inferior_horizontal_ref(
		{ 0,-roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2));
	CylinderQuarter firstInferiorHorizontalCylinderQuarter = { first_inferior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(firstInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame second_inferior_horizontal_ref(
		{ 0,-roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), -PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2)));
	CylinderQuarter secondInferiorHorizontalCylinderQuarter = { second_inferior_horizontal_ref, roundedBox.extents.x, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(secondInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame third_inferior_horizontal_ref(
		{ roundedBox.extents.x,-roundedBox.extents.z,0 },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	CylinderQuarter thirdInferiorHorizontalCylinderQuarter = { third_inferior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(thirdInferiorHorizontalCylinderQuarter, nSectors, false, color);

	ReferenceFrame fourth_inferior_horizontal_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.z,0 },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2)));
	CylinderQuarter fourthInferiorHorizontalCylinderQuarter = { fourth_inferior_horizontal_ref, roundedBox.extents.y, roundedBox.radius };
	MyDrawWireframeCylinderQuarter(fourthInferiorHorizontalCylinderQuarter, nSectors, false, color);

#pragma endregion

#pragma region Coin face supérieur

	ReferenceFrame first_superior_corner_ref(
		{ roundedBox.extents.x,roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));
	SphereCorner firstSuperiorCorner = { first_superior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(firstSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame second_superior_corner_ref(
		{ roundedBox.extents.x,roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2));
	SphereCorner secondSuperiorCorner = { second_superior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(secondSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame third_superior_corner_ref(
		{ -roundedBox.extents.x,roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), -PI / 2));
	SphereCorner thirdSuperiorCorner = { third_superior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(thirdSuperiorCorner, nSectors, nSectors, color);

	ReferenceFrame fourth_superior_corner_ref(
		{ -roundedBox.extents.x,roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI));
	SphereCorner fourthSuperiorCorner = { fourth_superior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(fourthSuperiorCorner, nSectors, nSectors, color);

#pragma endregion

#pragma region Coin face inférieur

	ReferenceFrame first_inferior_corner_ref(
		{ roundedBox.extents.x,-roundedBox.extents.z,roundedBox.extents.y },
		QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2));
	SphereCorner firstInferiorCorner = { first_inferior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(firstInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame second_inferior_corner_ref(
		{ roundedBox.extents.x,-roundedBox.extents.z,-roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI / 2),
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 2)));
	SphereCorner secondInferiorCorner = { second_inferior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(secondInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame third_inferior_corner_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.z,roundedBox.extents.y },
		QuaternionMultiply(
			QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI),
			QuaternionFromAxisAngle(Vector3Normalize({ 0,1,0 }), PI)));
	SphereCorner thirdInferiorCorner = { third_inferior_corner_ref, roundedBox.radius };
	MyDrawWireframeSphereCorner(thirdInferiorCorner, nSectors, nSectors, color);

	ReferenceFrame fourth_inferior_corner_ref(
		{ -roundedBox.extents.x,-roundedBox.extents.z,-roundedBox.extents.y },
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
	Vector3 a = LocalToGlobalPos(segment.a, segment.ref);
	Vector3 b = LocalToGlobalPos(segment.b, segment.ref);

	DrawLine3D(a, b, color);
	MyDrawPolygonSphere({ {a,QuaternionIdentity()},.15f }, 16, 8, colorA);
	MyDrawPolygonSphere({ {b,QuaternionIdentity()},.15f }, 16, 8, colorB);

}

#pragma endregion

#pragma endregion

#pragma region Intersections

#pragma region IntersectLinePlane

bool IntersectLinePlane(Line line, Plane plane, float& t, Vector3& interPt, Vector3& interNormal)
{
	// no intersection if line is parallel to the plane
	float dotProd = Vector3DotProduct(plane.normal, line.dir);
	if (fabsf(dotProd) < EPSILON) return false;

	// intersection: t, interPt & interNormal
	t = (plane.d - Vector3DotProduct(plane.normal, line.pt)) / dotProd;
	interPt = Vector3Add(line.pt, Vector3Scale(line.dir, t)); // OM = OA+tAB
	interNormal = Vector3Scale(plane.normal,Vector3DotProduct(Vector3Subtract(line.pt, interPt), plane.normal) < 0 ? -1.f : 1.f);
	return true;
}

bool IntersectSegmentPlane(Segment segment, Plane plane, float& t, Vector3& interPt, Vector3& interNormal)
{
	// create a line from the two points of the segment
	Line line;
	line.pt = segment.a;
	line.dir = Vector3Subtract(segment.b, segment.a);

	// check for intersection between the line and plane
	if (!IntersectLinePlane(line, plane, t, interPt, interNormal)) return false;

	// check if intersection point is between the two points of the segment
	if (t < 0 || t > 1) return false;

	return true;
}

bool IntersectSegmentQuad(Segment segment, Quad quad, float& t, Vector3& interPt, Vector3& interNormal)
{
	// create a plane from the quad's reference frame
	Plane plane;
	plane.normal = quad.ref.origin;
	plane.d = Vector3DotProduct(quad.ref.origin, plane.normal);

	// check for intersection between the segment and the plane
	if (!IntersectSegmentPlane(segment, plane, t, interPt, interNormal)) return false;

	// check if the intersection point is inside the quad
	if (fabsf(interPt.x) > fabsf(quad.extents.x) || fabsf(interPt.z) > fabsf(quad.extents.z)) return false;

	return true;
}

bool IntersectSegmentSphere(Segment seg, Sphere s, float& t, Vector3& interPt, Vector3& interNormal)
{
	// revérifier tous les référentiels
	Vector3 OA = LocalToGlobalPos(seg.a, seg.ref);
	Vector3 AB = Vector3Subtract(seg.b, seg.a);

	Vector3 A_sphereRef = GlobalToLocalPos(LocalToGlobalPos(seg.a, seg.ref), s.ref);
	float AB_dot = Vector3DotProduct(AB, AB);
	float A_sphereRef_Dot = Vector3DotProduct(AB, A_sphereRef);
	t = A_sphereRef_Dot / AB_dot;
	cout << t << "\n";
	if (t < 0 || t > 1) return false;

	Vector3 OM = Vector3Add(OA, Vector3Scale(AB, t));
	interPt = LocalToGlobalPos(OM, s.ref);

	cout << "a: " << Vector3DotProduct(, OM) << ", b: " << (s.radius * s.radius) << "\n";
	if (Vector3DotProduct(OM, OM) == (s.radius * s.radius))
	{
		return true;
	}

	// 
	//if (t < 0 || t > 1) return false;
	//interPt = LocalToGlobalPos(Vector3Add(seg.a, Vector3Scale(segDir, t)), s.ref);
	//Vector3 spherePos = Vector3RotateByQuaternion(s.ref.origin, s.ref.q);
	//float dist = Vector3Distance(interPt, spherePos);
	//if (dist > s.radius) return false;

	//interNormal = Vector3Normalize(Vector3Subtract(interPt, spherePos));
	return true;
}

//bool IntersectSegmentQuad(Segment seg, Quad quad, float& t, Vector3& interPt, Vector3& interNormal)
//{
//	// Check for intersection between the segment and the plane of the quad
//	if (!IntersectSegmentPlane(seg, quad.plane, t, interPt, interNormal)) return false;
//
//	// check if intersection point is inside the quad
//	Vector3 quadPt = Vector3TransformCoord(interPt, quad.ref);
//	if (quadPt.x < -quad.extents.x || quadPt.x > quad.extents.x ||
//		quadPt.y < -quad.extents.y || quadPt.y > quad.extents.y)
//	{
//		return false;
//	}
//	return true;
//}

#pragma endregion

#pragma endregion


int main(int argc, char* argv[])
{
	// ICI LES TEST

	#pragma region Init

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

	#pragma endregion

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
			DrawGrid(20, 1); // Draw a grid
			DrawLine3D({ 0 }, { 0,10,0 }, DARKGRAY);
			DrawSphere({ 10,0,0 }, .2f, RED);
			DrawSphere({ 0,10,0 }, .2f, GREEN);
			DrawSphere({ 0,0,10 }, .2f, BLUE);

			ReferenceFrame ref;
			static float angle = 0;
			angle += 0.05;
			Vector3 axes = { 1, 1, 1 };

			#pragma region display tests

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

			#pragma region methods testing

			//Vector3 initialPos = { 4,0,0 };
			//ref = ReferenceFrame(
			//	initialPos,
			//	QuaternionFromAxisAngle(Vector3Normalize({1,1,1}), PI/2));
			//Vector3 global = LocalToGlobalVect({ 2,1,6 }, ref);
			//Vector3 local = GlobalToLocalVect(global, ref);

			//cout << "global: {" << global.x << ", " << global.y << ", " << global.z << "}\n";
			//cout << "local: {" << local.x << ", " << local.y << ", " << local.z << "}\n";

			//global = LocalToGlobalVect(local, ref);
			//cout << "newGlobal: {" << global.x << ", " << global.y << ", " << global.z << "}\n";

			//MyDrawPolygonSphere({ {initialPos, QuaternionIdentity()},.15f }, 16, 8, BLUE);
			//MyDrawPolygonSphere({ {global, QuaternionIdentity()},.15f }, 16, 8, RED);

			//Vector3 unitVect = { 4, 7, 1 };
			//Vector3 lineOrigin = { 0, 0, 0 };
			//Segment segment = { ref, lineOrigin, unitVect };

			//Vector3 point = { 2,4,5 };
			//Vector3 proj = ProjectedPointOnLine(lineOrigin, Vector3Normalize(unitVect), point);
			////cout << "proj: {" << proj.x << ", " << proj.y << ", " << proj.z << "}\n";
			//Segment projSegment = { ref, point, proj };

			//MyDrawSegment(segment);
			//MyDrawSegment(projSegment);

			#pragma endregion

			#pragma region Intersections

			//TESTS INTERSECTIONS
			Vector3 interPt;
			Vector3 interNormal;
			float t;

			ref = ReferenceFrame(
				{ 5,5,5 },
				QuaternionFromAxisAngle(Vector3Normalize({0,0,0}), 0));

			// TEST LINE PLANE INTERSECTION
			/*Segment segment = { ref, {-5,8,0},{5,-8,3} };
			MyDrawSegment(segment);
			Plane plane = { Vector3RotateByQuaternion({0,1,0}, QuaternionFromAxisAngle({1,0,0},time * .5f)), 2 };
			ReferenceFrame refQuad = { Vector3Scale(plane.normal, plane.d),
									   QuaternionFromVector3ToVector3({0,1,0},plane.normal) };
			Quad quad = { refQuad,{10,1,10} };
			MyDrawQuad(quad);
			Line line = { segment.a,Vector3Subtract(segment.b,segment.a) };
			if (IntersectLinePlane(line, plane, t, interPt, interNormal))
			{
				cout << "newGlobal: {" << "PAS D INTERSECTION" << "}\n";
				MyDrawPolygonSphere({ {interPt,QuaternionIdentity()},.1f }, 16, 8, BLACK);
				DrawLine3D(interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt), RED);
			}*/

			// TEST SEGMENT PLANE INTERSECTION
			//Segment segment = { ref, {-5,8,0},{5,-8,3} };
			//MyDrawSegment(segment);
			//Plane plane = { Vector3RotateByQuaternion({0,1,0}, QuaternionFromAxisAngle({1,0,0},time * .5f)), 2 };
			//ReferenceFrame refQuad = { Vector3Scale(plane.normal, plane.d),
			//						   QuaternionFromVector3ToVector3({0,1,0},plane.normal) };
			//Quad quad = { refQuad,{10,1,10} };
			//MyDrawQuad(quad);
			//if (IntersectSegmentPlane(segment, plane, t, interPt, interNormal))
			//{
			//	if (t) printf("True\n");
			//	MyDrawPolygonSphere({ {interPt,QuaternionIdentity()},.1f }, 16, 8, RED);
			//	DrawLine3D(interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt), RED);
			//}

			// TEST SEGMENT SPHERE INTERSECTION
			Sphere sphere = { ref, 2 };
			MyDrawSphere(sphere, 50, 50);
			Segment seg = { ref, {1.5,-1,0}, {1.5,5,0} };
			MyDrawSegment(seg);
			if (IntersectSegmentSphere(seg, sphere, t, interPt, interNormal)) {
				cout << "interPt: {" << interPt.x << ", " << interPt.y << ", " << interPt.z << "}\n";
				MyDrawPolygonSphere({ {interPt,QuaternionIdentity()},.1f }, 16, 8, BLUE);
				if (t) printf("True\n");
			}

			//DrawLine3D(interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt), RED);
			//Plane plane = { Vector3RotateByQuaternion({0,1,0}, QuaternionFromAxisAngle({1,0,0},time * .5f)), 2 };
			//ReferenceFrame refQuad = { Vector3Scale(plane.normal, plane.d),
			//						   QuaternionFromVector3ToVector3({0,1,0},plane.normal) };
			//Quad quad = { refQuad,{10,1,10} };
			//MyDrawQuad(quad);
			//if (IntersectSegmentPlane(segment, plane, t, interPt, interNormal))
			//{
			//	if (t) printf("True\n");
			//	MyDrawPolygonSphere({ {interPt,QuaternionIdentity()},.1f }, 16, 8, RED);
			//	DrawLine3D(interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt), RED);
			//}

			#pragma endregion

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