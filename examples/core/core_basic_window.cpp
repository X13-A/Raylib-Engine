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
	bool w = localPt.x >= -box.extents.x && localPt.x <= box.extents.x;
	bool h = localPt.y >= -box.extents.y && localPt.y <= box.extents.y;
	bool d = localPt.z >= -box.extents.z && localPt.z <= box.extents.z;
	return (w && h && d);
}

Vector3 NormalFromRef(ReferenceFrame ref)
{
	// Vérifier le sens de la normale (derrière le quad par ex)
	return LocalToGlobalPos({0,1,0}, ref);
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

	//roundedBox.extents.x /= 2;
	//roundedBox.extents.y /= 2;
	//roundedBox.extents.z /= 2;
	//roundedBox.radius /= 2;

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

	//roundedBox.extents.x /= 2;
	//roundedBox.extents.y /= 2;
	//roundedBox.extents.z /= 2;
	//roundedBox.radius /= 2;

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
	Vector3 a = GlobalToLocalPos(segment.a, quad.ref);
	Vector3 b = GlobalToLocalPos(segment.b, quad.ref);
	if (!((a.y >= 0 && b.y <= 0) || (b.y >= 0 && a.y <= 0)))
	{
		return false;
	}
	if ((b.y - a.y) == 0) return false;
	t = -a.y / (b.y - a.y);
	
	interPt = Vector3Add(a, Vector3Scale(Vector3Subtract(b, a), t));

	if (fabs(interPt.x) > fabs(quad.extents.x) || fabs(interPt.z) > fabs(quad.extents.z)) return false;

	interPt = LocalToGlobalPos(interPt ,quad.ref);

	if (a.y < 0) interNormal = Vector3RotateByQuaternion({ 0,-1,0 }, quad.ref.q);	
	else interNormal = Vector3RotateByQuaternion({ 0,1,0 }, quad.ref.q);

	return true;
}

bool IntersectSegmentDisk(Segment segment, Disk disk, float& t, Vector3& interPt, Vector3& interNormal)
{
	Vector3 a = GlobalToLocalPos(segment.a, disk.ref);
	Vector3 b = GlobalToLocalPos(segment.b, disk.ref);
	if (!((a.y >= 0 && b.y <= 0) || (b.y >= 0 && a.y <= 0)))
	{
		return false;
	}
	t = -a.y / (b.y - a.y);

	interPt = Vector3Add(a, Vector3Scale(Vector3Subtract(b, a), t));

	if (Vector3Distance(interPt, {0,0,0}) > disk.radius) return false;

	interPt = LocalToGlobalPos(interPt, disk.ref);

	if (a.y < 0) interNormal = Vector3RotateByQuaternion({ 0,-1,0 }, disk.ref.q);
	else interNormal = Vector3RotateByQuaternion({ 0,1,0 }, disk.ref.q);

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

	// Calcul des ratios menant aux points d'intersection
	float t0 = Vector3Distance(Ro, tp) / Vector3Distance(seg.a, seg.b);
	float t1 = t0 - x / Vector3Distance(seg.a, seg.b);
	float t2 = t0 + x / Vector3Distance(seg.a, seg.b);

	// Calcul des points d'intersection
	Vector3 tp1 = Vector3Add(Ro, Vector3Scale(Rd, t1));
	Vector3 tp2 = Vector3Add(Ro, Vector3Scale(Rd, t2));

	// Calcul des distances entre les points et l'origine du segment
	float tp1Dist = Vector3Distance(tp1, Ro);
	float tp2Dist = Vector3Distance(tp2, Ro);

	// Calcul du point le plus proche
	if (tp1Dist < tp2Dist)
	{
		t = t1;
		interPt = tp1;
	}
	else
	{
		t = t2;
		interPt = tp2;
	}
	
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

		t = min(t1, t2);

		if (t > 1 || t < 0) return false;

		interPt = Vector3Add(segment.a, Vector3Scale(Vector3Subtract(segment.b, segment.a), t));
		Vector3 interPtLocal = GlobalToLocalPos(interPt, cylinder.ref);
		interNormal = Vector3Subtract(interPtLocal, { 0, interPtLocal.y, 0 });
		interNormal = LocalToGlobalPos(Vector3Add(interPtLocal, interNormal), cylinder.ref);
		interNormal = Vector3Normalize(Vector3Subtract(interNormal, interPt));
		return true;
	}

	if (fabs(delta) <= EPSILON)
	{
		t = -b / (2 * a);
		return true;
	}

	return false;
}

bool IntersectSegmentCylinder(Segment segment, Cylinder cylinder, float& t, Vector3& interPt, Vector3& interNormal)
{
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

	return intersect;
}

bool IntersectSegmentCapsule(Segment segment, Capsule capsule, float& t, Vector3& interPt, Vector3& interNormal)
{
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
	left_front_vertical_ref.Translate(Vector3RotateByQuaternion({ -rndBox.extents.x ,0,rndBox.extents.z  }, left_front_vertical_ref.q));
	left_front_vertical_ref.q = QuaternionMultiply(left_front_vertical_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, 3*PI/2));
	Cylinder left_front_vertical_cylinder = { left_front_vertical_ref, rndBox.extents.y, rndBox.radius };

	ReferenceFrame right_front_vertical_ref = rndBox.ref;
	right_front_vertical_ref.Translate(Vector3RotateByQuaternion({ rndBox.extents.x ,0,rndBox.extents.z }, right_front_vertical_ref.q));
	Cylinder right_front_vertical_cylinder = { right_front_vertical_ref, rndBox.extents.y, rndBox.radius };

	ReferenceFrame right_back_vertical_ref = rndBox.ref;
	right_back_vertical_ref.Translate(Vector3RotateByQuaternion({ rndBox.extents.x ,0,-rndBox.extents.z }, right_back_vertical_ref.q));
	right_back_vertical_ref.q = QuaternionMultiply(right_back_vertical_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, PI/2));
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
	top_front_left_corner_ref.q = QuaternionMultiply(top_front_left_corner_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, -PI/2));
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
	bottom_front_right_corner_ref.q = QuaternionMultiply(bottom_front_right_corner_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, PI/2));
	Sphere bottom_front_right_sphere = { bottom_front_right_corner_ref, rndBox.radius };

	ReferenceFrame bottom_back_left_corner_ref = rndBox.ref;
	bottom_back_left_corner_ref.Translate(Vector3RotateByQuaternion({ -rndBox.extents.x, -rndBox.extents.y, -rndBox.extents.z }, bottom_back_left_corner_ref.q));
	bottom_back_left_corner_ref.q = QuaternionMultiply(bottom_back_left_corner_ref.q, QuaternionFromAxisAngle({ 1,0,0 }, PI));
	bottom_back_left_corner_ref.q = QuaternionMultiply(bottom_back_left_corner_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, -PI/2));
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
	front_bottom_horizontal_ref.q = QuaternionMultiply(front_bottom_horizontal_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, -PI/2));
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
	right_bottom_horizontal_ref.q = QuaternionMultiply(right_bottom_horizontal_ref.q, QuaternionFromAxisAngle({ 0,1,0 }, -PI/2));
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

	float tempT;
	Vector3 tempInterPt;
	Vector3 tempInterNormal;
	float minDist = INFINITY;
	float dist;
	bool intersect = false;


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

	for (int i = 0; i < 12; i++) {
		if (IntersectSegmentCylinder(seg, cylinders[i], tempT, tempInterPt, tempInterNormal))
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

	for (int i = 0; i < 8; i++) {
		if (IntersectSegmentSphere(seg, spheres[i], tempT, tempInterPt, tempInterNormal))
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

bool GettingCloseToDest(Vector3 pos, Vector3 acc, float deltaTime, Vector3 dest)
{
	float actualNorm = Vector3Distance(pos, dest);
	float nextNorm = Vector3Distance(Vector3Add(pos, Vector3Scale(acc, deltaTime)), dest);

	if (nextNorm < actualNorm) return true;
	return false;
}

bool WillCollideWithBox(Vector3 pos, Vector3 acc, float deltaTime, Box box, string &face)
{
	Vector3 dest = Vector3Add(pos, Vector3Scale(acc, deltaTime * 10));
	Segment dir = { pos, dest };


	float t;
	Vector3 interPt;
	Vector3 interNormal;

	if (!IntersectSegmentBox(dir, box, t, interPt, interNormal)) return false;

	MyDrawSegment({ interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt) });

	return true;
}

bool GetSphereNewPositionAndVelocityIfCollidingWithRoundedBox(Sphere sphere, RoundedBox rndBox, Vector3 velocity, float deltaTime, float& colT, Vector3& colSpherePos, Vector3& colNormal, Vector3& newPosition, Vector3& newVelocity)
{
	RoundedBox minkowski = { rndBox.ref, {rndBox.extents.x + rndBox.radius, rndBox.extents.y + rndBox.radius, rndBox.extents.z + rndBox.radius}, sphere.radius};
	Vector3 A = sphere.ref.origin;
	Vector3 B = Vector3Add(sphere.ref.origin, Vector3Scale(velocity, deltaTime));
	Segment AB = { A, B };
	Vector3 colPt;

	if (IntersectSegmentRoundedBox(AB, minkowski, colT, colPt, colNormal))
	{
		newVelocity = Vector3Reflect(velocity, colNormal); // Changement de direction
		newPosition = Vector3Add(colPt, Vector3Scale(newVelocity, deltaTime * (1 - colT)));  // Ajout de la distance restante dans la bonne direction
		return true;
		return true;
	}
	//Box OBB = { minkowski.ref, {minkowski.extents.x + minkowski.radius, minkowski.extents.y + minkowski.radius, minkowski.extents.z + minkowski.radius} };
	//
	//Vector3 A = sphere.ref.origin;
	//Vector3 B = Vector3Add(sphere.ref.origin, Vector3Scale(velocity, deltaTime));
	//Segment AB = { A, B };
	//
	////MyDrawWireframeBox(OBB, BLUE);
	////MyDrawWireframeRoundedBox(minkowski, 10, RED);
	//
	//Vector3 colPt;
	//if (IntersectSegmentBox(AB, OBB, colT, colPt, colNormal))
	//{
	//	Vector3 localColPt = GlobalToLocalPos(colPt, OBB.ref);
	//	char nearestFace = GetNearestFace(colPt, OBB);
	//
	//	bool inX = false;
	//	bool inY = false;
	//	bool inZ = false;
	//
	//	if (nearestFace == 'L' || nearestFace == 'R')
	//	{
	//		inX = fabs(localColPt.x) <= rndBox.extents.x + rndBox.radius + sphere.radius + EPSILON;
	//		inY = fabs(localColPt.y) <= rndBox.extents.y + rndBox.radius + EPSILON;
	//		inZ = fabs(localColPt.z) <= rndBox.extents.z + rndBox.radius + EPSILON;
	//	}
	//	else if (nearestFace == 'T' || nearestFace == 'B')
	//	{
	//		inX = fabs(localColPt.x) <= rndBox.extents.x + rndBox.radius + EPSILON;
	//		inY = fabs(localColPt.y) <= rndBox.extents.y + rndBox.radius + sphere.radius + EPSILON;
	//		inZ = fabs(localColPt.z) <= rndBox.extents.z + rndBox.radius + EPSILON;
	//	}
	//	else if (nearestFace == 'f' || nearestFace == 'b')
	//	{
	//		inX = fabs(localColPt.x) <= rndBox.extents.x + rndBox.radius + EPSILON;
	//		inY = fabs(localColPt.y) <= rndBox.extents.y + rndBox.radius + EPSILON;
	//		inZ = fabs(localColPt.z) <= rndBox.extents.z + rndBox.radius + sphere.radius + EPSILON;
	//	}
	//
	//	if (!(inX && inY && inZ)) return false;
	//	newVelocity = Vector3Reflect(velocity, colNormal); // Changement de direction
	//	newPosition = Vector3Add(colPt, Vector3Scale(newVelocity, deltaTime * (1-colT)));  // Ajout de la distance restante dans la bonne direction
	//	return true;
	//}
	//else if (IsPointInsideBox(OBB, B))
	//{
	//	if (IntersectSegmentRoundedBox(AB, minkowski, colT, colPt, colNormal))
	//	{
	//		newVelocity = Vector3Reflect(velocity, colNormal); // Changement de direction
	//		newPosition = Vector3Add(colPt, Vector3Scale(newVelocity, deltaTime * (1 - colT)));  // Ajout de la distance restante dans la bonne direction
	//		return true;
	//	}
	//}
	return false;
}

bool GetSphereNewPositionAndVelocityIfCollidingWithRoundedBoxes(Sphere sphere, const std::vector<RoundedBox>& rndBoxes, Vector3 velocity, float deltaTime, float& colT, Vector3& colSpherePos, Vector3& colNormal, Vector3& newPosition, Vector3& newVelocity)
{
	colT = FLT_MAX;
	bool collided = false;
	for (const RoundedBox &rndBox : rndBoxes)
	{
		float t;
		Vector3 colPt, normal, newPos, newVel;
		if (GetSphereNewPositionAndVelocityIfCollidingWithRoundedBox(sphere, rndBox, velocity, deltaTime, t, colPt, normal, newPos, newVel))
		{
			if (t < colT)
			{
				colT = t;
				colSpherePos = colPt;
				colNormal = normal;
				newPosition = newPos;
				newVelocity = newVel;
				collided = true;
			}
		}
	}
	return collided;
}

#pragma endregion

Vector3 GetAcc(Vector3 vel)
{
	Vector3 G = { 0, -0, 0 };
	return Vector3Add(G, vel);
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

	SetTargetFPS(240);

	// CAMERA 
	Vector3 cameraPos = { 8, 15, 14 };
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
	
	Vector3 initialPos = { 0, 2, 0 };
	Vector3 initialVel = { 16, 8, 0 };
	
	Vector3 pos = initialPos;
	Vector3 vel = initialVel;
	Vector3 rot = { 0, 0, 0 };
	float mass = 10;
	float radius = 1;

	#pragma endregion

	// Main game loopx
	while (!WindowShouldClose()) // Detect window close button or ESC key
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
			if (!IsMouseButtonDown(MOUSE_LEFT_BUTTON))
			{
				//3D REFERENTIAL
				DrawGrid(39, 1); // Draw a grid
				DrawLine3D({ 0 }, { 0,18.5,0 }, DARKGRAY);
				DrawSphere({ 18.5,0,0 }, .2f, RED);
				DrawSphere({ 0,18.5,0 }, .2f, GREEN);
				DrawSphere({ 0,0,18.5 }, .2f, BLUE);
			}

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

			#pragma region intersections tests

			//TESTS INTERSECTIONS
			//Vector3 interPt;
			//Vector3 interNormal;
			//float t;

			//ref = ReferenceFrame(
			//	{ 0,0,0 },
			//	QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), 0));

			// TEST LINE PLANE INTERSECTION
			/*Segment segment = { {-5,8,0},{5,-8,3} };
			MyDrawSegment(segment);
			Plane plane = { Vector3RotateByQuaternion({10,1,0}, QuaternionFromAxisAngle({1,0,0},time * .5f)), 2 };
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
			/*Segment segment = { {-5,8,0},{5,-8,3} };
			MyDrawSegment(segment);
			Plane plane = { Vector3RotateByQuaternion({2,4,0}, QuaternionFromAxisAngle({1,1,1},time * .5f)), 2 };
			ReferenceFrame refQuad = { Vector3Scale(plane.normal, plane.d),
									   QuaternionFromVector3ToVector3({0,1,0},plane.normal) };
			Quad quad = { refQuad,{10,1,10} };
			MyDrawQuad(quad);
			if (IntersectSegmentPlane(segment, plane, t, interPt, interNormal))
			{
				if (t) printf("True\n");
				MyDrawPolygonSphere({ {interPt,QuaternionIdentity()},.1f }, 16, 8, RED);
				DrawLine3D(interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt), RED);
			}*/

			// TEST SEGMENT QUAD INTERSECTION
			/*Segment segment = { {-5,8,0},{5,-8,3} };
			MyDrawSegment(segment);
			ReferenceFrame refQuad = ReferenceFrame({ 0,2,0 }, QuaternionFromAxisAngle({ 0,0,1 }, time));
			Quad quad = { refQuad,{4,1,4} };
			MyDrawQuad(quad);
			bool test = IntersectSegmentQuad(segment, quad, t, interPt, interNormal);
			if (test)
			{
				MyDrawPolygonSphere({ {interPt,QuaternionIdentity()},.1f }, 16, 8, RED);
				DrawLine3D(interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt), RED);
			}*/

			//MyDrawPolygonSphere({ refQuad, 0.05f }, 16, 16, BLACK);
			/*Vector3 i = LocalToGlobalPos(refQuad.i, refQuad);
			Vector3 j = NormalFromRef(refQuad);
			Vector3 k = LocalToGlobalPos(refQuad.k, refQuad);
			DrawLine3D(refQuad.origin, i, RED);
			DrawLine3D(refQuad.origin, j, GREEN);
			DrawLine3D(refQuad.origin, k, BLUE);
			MyDrawPolygonSphere({ ReferenceFrame(i, QuaternionIdentity()), 0.05f }, 16, 16, RED);
			MyDrawPolygonSphere({ ReferenceFrame(j, QuaternionIdentity()), 0.05f }, 16, 16, GREEN);
			MyDrawPolygonSphere({ ReferenceFrame(k, QuaternionIdentity()), 0.05f }, 16, 16, BLUE);

			 TEST SEGMENT SPHERE INTERSECTION SIMPLE
			Segment segment = { ref, {-5,0,0}, {5,0,0} };
			MyDrawSegment(segment);
			ReferenceFrame refSphere = { {0,0,0}, QuaternionIdentity() };
			Sphere sphere = { refSphere, 1 };
			MyDrawSphere(sphere, 15, 15);
			bool test = IntersectSegmentSphere(segment, sphere, t, interPt, interNormal);
			if (test)
			{
				MyDrawPolygonSphere({ {interPt,QuaternionIdentity()},.1f }, 16, 8, BLUE);
				DrawLine3D(interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt), RED);
			}
			cout << "res: " << test << "\n";*/

			// TEST SEGMENT SPHERE INTERSECTION ROTATE
			/*Segment segment = { ref,{5,-5,10}, {1,6,0} };
			MyDrawSegment(segment);
			Plane plane = { Vector3RotateByQuaternion({0,1,0}, QuaternionFromAxisAngle({1,0,0},time * .5f)), 2 };
			ReferenceFrame refSphere = { Vector3Scale(plane.normal, plane.d),
									   QuaternionFromVector3ToVector3({0,1,0},plane.normal) };
			Sphere sphere = { refSphere, 5 };
			MyDrawSphere(sphere, 15, 15);
			bool test = IntersectSegmentSphere(segment, sphere, t, interPt, interNormal);
			if (test)
			{
				MyDrawPolygonSphere({ {interPt,QuaternionIdentity()},.1f }, 16, 8, BLUE);
				DrawLine3D(interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt), RED);
			}
			cout << "intersection: " << test << "\n";*/


			// TEST SEGMENT DISK INTERSECTION
			/*Segment segment = { {-5,8,0},{5,-8,3} };
			MyDrawSegment(segment);
			Plane plane = { Vector3RotateByQuaternion({0,1,0}, QuaternionFromAxisAngle({1,0,0},time * .5f)), 2 };
			ReferenceFrame refDisk = { Vector3Scale(plane.normal, plane.d),
									   QuaternionFromVector3ToVector3({0,1,0},plane.normal) };
			Disk disk = { refDisk, 5 };
			MyDrawDisk(disk, 15);
			bool test = IntersectSegmentDisk(segment, disk, t, interPt, interNormal);
			if (test)
			{
				MyDrawPolygonSphere({ {interPt,QuaternionIdentity()},.1f }, 16, 8, RED);
				DrawLine3D(interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt), RED);
			}
			cout << "res: " << test << "\n";*/

			// TEST SEGMENT BOX INTERSECTION
			/*Segment segment = { {5,-5,10}, {1,6,-5} };
			MyDrawSegment(segment);
			Plane plane = { Vector3RotateByQuaternion({0,1,0}, QuaternionFromAxisAngle({1,0,0},time * .5f)), 2 };
			ReferenceFrame refBox = { Vector3Scale(plane.normal, plane.d),
									   QuaternionFromVector3ToVector3({0,1,0},plane.normal) };
			Box box = { refBox, {5,2,3} };
			MyDrawBox(box, false, true);
			bool test = IntersectSegmentBox(segment, box, t, interPt, interNormal);
			if (test)
			{
				MyDrawPolygonSphere({ {interPt,QuaternionIdentity()},.1f }, 16, 8, BLUE);
				DrawLine3D(interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt), RED);
			}
			cout << "res: " << test << "\n";*/

			// TEST SEGMENT CYLINDER INTERSECTION
			//ReferenceFrame refCylindre = { {0,0,1}, QuaternionFromAxisAngle({1,0,0}, time)};
			//Cylinder cylindre = { refCylindre, 4, 1 };
			//MyDrawCylinder(cylindre, 15);
			//Segment segment = { {-5,8,0},{1, 1, 1}};
			//MyDrawSegment(segment);
			//bool test = IntersectSegmentCylinder(segment, cylindre, t, interPt, interNormal);
			//if (test)
			//{
			//	MyDrawPolygonSphere({ {interPt,QuaternionIdentity()},.1f }, 16, 8, RED);
			//	DrawLine3D(interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt), RED);
			//}

			// TEST SEGMENT CAPSULE INTERSECTION
			//ReferenceFrame refCapsule = { {0,0,1}, QuaternionFromAxisAngle({0,0,1}, time)};
			//Capsule capsule = { refCapsule, 4, 1 };
			//MyDrawCapsule(capsule, 15, 15);
			//Segment segment = { {0,8,1},{0,-8,1} };
			//bool test = IntersectSegmentCapsule(segment, capsule, t, interPt, interNormal);
			//if (test)
			//{
			//	MyDrawPolygonSphere({ {interPt,QuaternionIdentity()},.1f }, 16, 8, RED);
			//	DrawLine3D(interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt), RED);
			//}
			//cout << "res: " << test << "\n";
			

			// TEST SEGMENT ROUNDED BOX INTERSECTION
			//ReferenceFrame refRndBox = { {3,4,1}, QuaternionFromAxisAngle({1,1,1}, time) };
			//RoundedBox rndBox = { refRndBox, {1,2,3}, 1};
			//MyDrawRoundedBox(rndBox, 5, true, true);
			//Segment segment = { {0,8,1},{0,-8,1} };
			//MyDrawSegment(segment);
			//bool test = IntersectSegmentRoundedBox(segment, rndBox, t, interPt, interNormal);
			//if (test)
			//{
			//	MyDrawPolygonSphere({ {interPt,QuaternionIdentity()},.1f }, 16, 8, RED);
			//	DrawLine3D(interPt, Vector3Add(Vector3Scale(interNormal, 1), interPt), RED);
			//}
			//cout << "res: " << test << "\n";

			#pragma endregion
			
			#pragma region collision testing
			int colorSpeed = 5;
			float colorPercentage = (sin(time * colorSpeed) + 1) / 2;
			unsigned char r = colorPercentage * 153;
			unsigned char g = colorPercentage * 204;
			unsigned char b = colorPercentage * 255;
			Color highlightColor = { r, g, b, 255 };
			//printf("%d, %d, %d\n", r, g, b);
			Sphere sphere = { {pos, QuaternionIdentity() }, radius };

			std::vector<RoundedBox> boxes =
			{
				{ { { 0,-1,0 }, QuaternionIdentity()}, { 20, 0.2, 20 }, 0, LIGHTGRAY }, // Floor
				{ { { 0,9,0 }, QuaternionIdentity()}, { 20, 0.2, 20 }, 0, BLANK}, // Ceiling
				{ { { -20,4,0 }, QuaternionFromAxisAngle({0,1,0}, PI / 2)}, {20, 5, 0.2}, 0, GRAY}, // Wall
				{ { { 20,4,0 }, QuaternionFromAxisAngle({0,1,0}, -PI / 2)}, {20, 5, 0.2}, 0, GRAY}, // Wall
				{ { { 0,4,20 }, QuaternionFromAxisAngle({0,0,0}, 0)}, {20, 5, 0.2}, 0, GRAY}, // Wall
				{ { { 0,4,-20 }, QuaternionFromAxisAngle({0,0,0}, 0)}, {20, 5, 0.2}, 0, GRAY}, // Wall
				
				{ { { 8, 3, 0}, QuaternionIdentity() }, {1,1,1}, 1.f, RED },
				{ { { -8, 3, 0}, QuaternionIdentity() }, {1,1,1}, 1.f, RED },
				{ { { 0, 3, 8}, QuaternionIdentity() }, {1,1,1}, 1.f, RED },
				{ { { 0, 3, -8}, QuaternionIdentity() }, {1,1,1}, 1.f, RED }
			};

			float colT;
			Vector3 colSpherePos;
			Vector3 colNormal;
			Vector3 newPos;
			Vector3 newVel;
			Vector3 colPt;

			bool collide = GetSphereNewPositionAndVelocityIfCollidingWithRoundedBoxes(sphere, boxes, vel, deltaTime, colT, colSpherePos, colNormal, newPos, newVel);
			if (collide)
			{
				vel = newVel;
				pos = newPos;
			}
			else
			{
				pos = Vector3Add(pos, Vector3Scale(vel, deltaTime));
			}

			MyDrawSphere(sphere, 30, 30, true, false, highlightColor);

			for (const RoundedBox& rndBox : boxes)
			{
				if (rndBox.color.a != 0)
				{
					MyDrawRoundedBox(rndBox, 10, true, false, rndBox.color);
				}
			}

			if (IsMouseButtonDown(MOUSE_MIDDLE_BUTTON))
			{
				pos = initialPos;
				vel = initialVel;
			}

#pragma region shootSegment

			if (IsMouseButtonDown(MOUSE_LEFT_BUTTON))
			{
				Vector3 ray = Vector3Scale(Vector3Subtract(camera.target, camera.position), 100);
				Segment shoot = { camera.position, Vector3Add(camera.position, ray) };

				for (const RoundedBox& rndBox : boxes)
				{
					if (rndBox.color.a == 0) continue;
					bool cameraTest = IntersectSegmentRoundedBox(shoot, rndBox, colT, colPt, colNormal);
					if (cameraTest)
					{
						MyDrawSphere({ {colPt,QuaternionIdentity()},.1f }, 16, 8, true, false, highlightColor);
						DrawLine3D(colPt, Vector3Add(Vector3Scale(colNormal, 1), colPt), highlightColor);
					}
				}
			}

			#pragma endregion

			printf("%f FPS\n", 1/deltaTime);
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